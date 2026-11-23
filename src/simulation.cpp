// The simulation itself. API is in occupancy/simulation.hpp.
//
// Pipeline: load morphologies -> project to 2D -> seed N cells with a random
// arbor rotation -> walk each soma over DIV, dragging its arbor, while the
// arbor reveals on the DIV schedule -> collect centrelines -> dilate once by
// the territory radius -> occupancy = covered / field area.
//
// WATCH THE RNG DRAW ORDER: reorder a draw and every logged number moves, pins
// included. Headings and seek decisions are drawn only when migrating, so
// mig_speed=0 leaves the stream alone, and the cue machinery draws nothing, so
// every coupling at weight 0 is bit-identical to no coupling.

#include "occupancy/simulation.hpp"
#include "occupancy/compartments.hpp"
#include "occupancy/distance_transform.hpp"
#include "occupancy/schedule.hpp"
#include "occupancy/segment.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace occupancy
{

namespace
{

/**
 * @brief Exact EDT over the scaffold geometry grid, in um.
 */
std::vector<double> geometry_edt(const std::vector<std::uint8_t> &cells,
                                 int gnx, int gny, double field_um, bool inside)
{
  const double comp_w_um = field_um / gnx;
  const double comp_h_um = field_um / gny;
  const double INF = 1e20;
  std::vector<double> out_um((std::size_t)gnx * gny, 0.0);
  // Pass 1 columns, pass 2 rows. The transform works in index units, so on a
  // non-square grid the two axes would come out in different um unless the
  // column pass is rescaled. Convert the y distances into column-width units
  // first, then one multiply by that width at the end gives um. thank you
  // google and claude for this one
  const int longest_side = std::max(gnx, gny);
  const double y_to_x_units = (comp_h_um / comp_w_um) * (comp_h_um / comp_w_um);
  std::vector<double> col_sq((std::size_t)gnx * gny);
  std::vector<double> in((std::size_t)longest_side),
      out((std::size_t)longest_side);
  std::vector<int> sv((std::size_t)longest_side);
  std::vector<double> sz((std::size_t)longest_side + 1);
  for (int x = 0; x < gnx; ++x)
  {
    for (int y = 0; y < gny; ++y)
    {
      const bool is_scaffold = cells[(std::size_t)y * gnx + x] != 0;
      in[(std::size_t)y] = (is_scaffold != inside) ? 0.0 : INF;
    }
    distance_transform_1d(in.data(), out.data(), gny, sv.data(), sz.data());
    for (int y = 0; y < gny; ++y)
      col_sq[(std::size_t)y * gnx + x] =
          std::min(INF, out[(std::size_t)y] * y_to_x_units);
  }
  for (int y = 0; y < gny; ++y)
  {
    for (int x = 0; x < gnx; ++x)
      in[(std::size_t)x] = col_sq[(std::size_t)y * gnx + x];
    distance_transform_1d(in.data(), out.data(), gnx, sv.data(), sz.data());
    for (int x = 0; x < gnx; ++x)
      out_um[(std::size_t)y * gnx + x] =
          std::sqrt(out[(std::size_t)x]) * comp_w_um;
  }
  return out_um;
}

} // namespace

OccupancyMode parse_occupancy_mode(const std::string &s)
{
  if (s == "inst")
    return OccupancyMode::Instant;
  if (s == "cum")
    return OccupancyMode::Cumulative;
  throw std::invalid_argument("unknown occupancy mode '" + s +
                              "' (expected cum or inst)");
}

const char *occupancy_mode_name(OccupancyMode m)
{
  return (m == OccupancyMode::Instant) ? "inst" : "cum";
}

std::string validate(const Params &p)
{
  std::ostringstream e;
  if (p.morphology_paths.empty())
    return "no morphology supplied";
  for (const auto &m : p.morphology_paths)
    if (m.empty())
      return "an empty morphology path was supplied";
  if (p.n_cells < 1)
    e << "n_cells must be >= 1 (got " << p.n_cells << "); ";
  if (p.scaffold_um <= 0)
    e << "scaffold_um must be > 0 (got " << p.scaffold_um << "); ";
  if (p.territory_um < 0)
    e << "territory_um must be >= 0 (got " << p.territory_um << "); ";
  if (p.resolution_um <= 0)
    e << "resolution_um must be > 0 (got " << p.resolution_um << "); ";
  if (p.div_days < 0)
    e << "div_days must be >= 0 (got " << p.div_days << "); ";
  if (p.resolution_um > 0 && p.scaffold_um > 0)
  {
    const long long grid_px = std::llround(p.scaffold_um / p.resolution_um);
    if (grid_px < 1)
      e << "a " << p.scaffold_um << " um field at --res " << p.resolution_um
        << " um/px rounds to a " << grid_px
        << "-pixel raster (occupancy 0/0) - --res is micrometres PER PIXEL; ";
  }
  if (p.seed_fraction <= 0 || p.seed_fraction > 1.0)
    e << "seed_fraction must be in (0,1] - it is a FRACTION of the scaffold "
         "side, not a radius in um (got "
      << p.seed_fraction << "); ";
  if (p.mig_speed_um_per_h < 0)
    e << "mig_speed must be >= 0 (got " << p.mig_speed_um_per_h << "); ";
  if (p.mig_persist_um <= 0)
    e << "mig_persist must be > 0 (got " << p.mig_persist_um << "); ";
  if (p.dt_hours <= 0)
    e << "dt_hours must be > 0 (got " << p.dt_hours << "); ";
  if (p.seek_weight < 0 || p.seek_weight > 1)
    e << "seek_weight must be in [0,1] (got " << p.seek_weight << "); ";
  if (p.p_seek < 0 || p.p_seek > 1)
    e << "p_seek must be in [0,1] (got " << p.p_seek << "); ";
  if (p.cohere_weight < 0 || p.cohere_weight > 1)
    e << "cohere_weight must be in [0,1] (got " << p.cohere_weight << "); ";
  // As with the stiffness couplings: a weight that could never act is refused
  // rather than returning a null that reads like biology.
  if (p.cohere_weight > 0 && p.mig_speed_um_per_h <= 0)
    e << "cohere_weight with migration off can never act - set mig_speed or "
         "drop the weight; ";
  if (p.cohere_weight > 0 && p.n_cells < 2)
    e << "cohere_weight needs at least 2 cells to have a neighbour (got "
      << p.n_cells << "); ";
  if (p.cohere_gate == CohereGate::Scaffold && p.scaffold_cells.empty())
    e << "cohere-gate=scaffold needs a scaffold geometry to gate on; ";
  if (p.cohere_gate == CohereGate::Scaffold && p.cohere_weight <= 0)
    e << "cohere-gate=scaffold with no cohere_weight can never act - set the "
         "weight or drop the gate; ";
  if (!p.scaffold_cells.empty())
  {
    if (p.scaffold_nx < 1 || p.scaffold_ny < 1 ||
        p.scaffold_cells.size() !=
            (std::size_t)p.scaffold_nx * (std::size_t)p.scaffold_ny)
      e << "scaffold_cells size does not match scaffold_nx*scaffold_ny; ";
    bool any_scaffold = false;
    for (auto v : p.scaffold_cells)
      any_scaffold = any_scaffold || v;
    if (!any_scaffold)
      e << "scaffold_cells has no scaffold compartment; ";
  }
  if (p.plate_scaffold && p.scaffold_cells.empty())
    e << "plate_scaffold is on but no scaffold geometry was supplied; ";
  if (p.seek_mode == SeekMode::Scaffold && p.scaffold_cells.empty())
    e << "seek=scaffold needs a scaffold geometry to steer toward; ";
  if (p.seek_mode == SeekMode::Scaffold && p.seek_weight > 0 &&
      p.mig_speed_um_per_h <= 0)
    e << "seek=scaffold with migration off can never steer - set mig_speed or "
         "drop the seek; ";
  if (p.stiff_weight < -1 || p.stiff_weight > 1)
    e << "stiff_weight must be in [-1,1] (got " << p.stiff_weight << "); ";
  if (p.stiff_speed_weight < 0 || p.stiff_speed_weight > 1)
    e << "stiff_speed_weight must be in [0,1] (got " << p.stiff_speed_weight
      << "); ";
  // nonzero, not > 0: a negative stiff_weight is up-gradient and needs maps too
  const bool stiff_coupled = p.stiff_weight != 0 || p.stiff_speed_weight > 0;
  if (stiff_coupled)
  {
    if (p.stiff_maps.empty())
      e << "a stiffness coupling weight is nonzero but no stiffness maps were "
           "supplied; ";
    else if (p.stiff_divs.size() != p.stiff_maps.size())
      e << "stiff_divs has " << p.stiff_divs.size() << " entries but there are "
        << p.stiff_maps.size() << " stiffness maps - one DIV per map; ";
    else
    {
      for (const auto &m : p.stiff_maps)
        if (m.size() != (std::size_t)p.stiff_nx * (std::size_t)p.stiff_ny)
        {
          e << "a stiffness map does not match stiff_nx*stiff_ny; ";
          break;
        }
      // run() takes the last div <= now, which is only the largest one if the
      // list ascends. The loader always ascends, so this catches hand-built
      // Params.
      for (std::size_t i = 1; i < p.stiff_divs.size(); ++i)
        if (p.stiff_divs[i] <= p.stiff_divs[i - 1])
        {
          e << "stiff_divs must be strictly ascending (entry " << i << " is "
            << p.stiff_divs[i] << " after " << p.stiff_divs[i - 1]
            << ") - the map in effect is selected by scanning them in order; ";
          break;
        }
      // A coupled run whose maps all start after the run ends would apply
      // nothing at any step and still get reported as coupled; that null
      // looked exactly like "stiffness does nothing", so reject it instead.
      if (!p.stiff_divs.empty() && p.stiff_divs.front() > p.div_days)
        e << "stiffness coupling requested but the earliest map (DIV "
          << p.stiff_divs.front() << ") is after the whole run (DIV "
          << p.div_days << ") - no map would ever be in effect; ";
      if (p.mig_speed_um_per_h <= 0)
        e << "stiffness coupling acts on migration, which is off (mig_speed "
             "0) - the coupling could never do anything; ";
    }
  }
  std::string s = e.str();
  if (!s.empty() && s.size() >= 2)
    s.erase(s.size() - 2);
  return s;
}

namespace
{
inline void rotate(double &x, double &y, double c, double s)
{
  double nx = c * x - s * y, ny = s * x + c * y;
  x = nx;
  y = ny;
}

/// Percentile of a sorted vector, linear interpolation, q in [0,1].
double percentile_sorted(const std::vector<double> &sorted, double q)
{
  if (sorted.empty())
    return 0.0;
  const double pos = q * (double)(sorted.size() - 1);
  const std::size_t lo = (std::size_t)pos;
  const std::size_t hi = std::min(lo + 1, sorted.size() - 1);
  const double frac = pos - (double)lo;
  return sorted[lo] + frac * (sorted[hi] - sorted[lo]);
}
} // namespace

Result run(const Params &p)
{
  std::string bad = validate(p);
  if (!bad.empty())
    throw std::invalid_argument("invalid parameters: " + bad);

  // Open every path so a typo gets caught, but only parse the templates we
  // actually sample: 3 cells from a pool of 50 is 3 parses.
  for (const auto &path : p.morphology_paths)
  {
    std::ifstream probe(path);
    if (!probe)
      throw std::runtime_error("cannot open morphology " + path);
  }
  const int n_templates = (int)p.morphology_paths.size();
  std::vector<MorphologyTemplate> pool((std::size_t)n_templates);
  std::vector<std::uint8_t> loaded((std::size_t)n_templates, 0);
  auto ensure_loaded = [&](int i)
  {
    if (loaded[(std::size_t)i])
      return;
    pool[(std::size_t)i] = load_template(p.morphology_paths[(std::size_t)i]);
    if (pool[(std::size_t)i].segments.empty())
      throw std::runtime_error("no segments loaded from " +
                               p.morphology_paths[(std::size_t)i]);
    loaded[(std::size_t)i] = 1;
  };

  const int n_cells = p.n_cells;
  const double field_um = p.scaffold_um;
  const double um_per_px = p.resolution_um;
  const int grid_px =
      (int)std::llround(field_um / um_per_px); // grid side, pixels
  const double territory_px = p.territory_um / um_per_px;
  // margin so arbors hanging off the field still dilate coverage back into it
  const int margin_px = (int)std::ceil(territory_px);
  const int padded_px = grid_px + 2 * margin_px;

  std::vector<std::uint8_t> centreline_px((std::size_t)padded_px * padded_px,
                                          0);

  std::mt19937 rng(p.seed);
  std::uniform_real_distribution<double> random_angle(0.0, 2 * M_PI);
  const double seed_half_um = 0.5 * p.seed_fraction * field_um;
  std::uniform_real_distribution<double> random_pos(-seed_half_um,
                                                    seed_half_um);
  std::normal_distribution<double> random_normal(0.0, 1.0);

  std::ofstream segments_file;
  const bool emit_segments = !p.segments_path.empty();
  if (emit_segments)
  {
    segments_file.open(p.segments_path);
    segments_file << "cell,x0,y0,x1,y1\n";
  }

  // Draw one cell's revealed arbor at (origin_x, origin_y), rotated by
  // (cos_rot, sin_rot). write_pixels sets centreline pixels, write_csv appends
  // to the segment CSV.
  auto rasterize = [&](int template_index, double cos_rot, double sin_rot,
                       double origin_x, double origin_y, double reveal_um,
                       bool write_pixels, bool write_csv, int cell_index)
  {
    // Sorted by d0, so once the reveal front is behind one segment it is behind
    // the rest too. Breaking here makes this O(revealed) instead of O(whole
    // arbor), and it runs once per cell per timestep.
    for (const auto &seg : pool[(std::size_t)template_index].segments)
    {
      if (seg.d0 > reveal_um)
        break;
      double frac_visible = 1.0;
      if (seg.d1 > reveal_um)
      {
        frac_visible = (reveal_um - seg.d0) / std::max(1e-9, (seg.d1 - seg.d0));
        if (frac_visible <= 0)
          continue;
      }
      double x0 = seg.x0, y0 = seg.y0;
      double x1 = seg.x0 + frac_visible * (seg.x1 - seg.x0);
      double y1 = seg.y0 + frac_visible * (seg.y1 - seg.y0);
      rotate(x0, y0, cos_rot, sin_rot);
      rotate(x1, y1, cos_rot, sin_rot);
      x0 += origin_x;
      y0 += origin_y;
      x1 += origin_x;
      y1 += origin_y;
      if (write_csv && emit_segments)
        segments_file << cell_index << "," << x0 << "," << y0 << "," << x1
                      << "," << y1 << "\n";
      if (!write_pixels)
        continue;
      double px0 = (x0 + field_um / 2) / um_per_px + margin_px,
             py0 = (y0 + field_um / 2) / um_per_px + margin_px;
      double px1 = (x1 + field_um / 2) / um_per_px + margin_px,
             py1 = (y1 + field_um / 2) / um_per_px + margin_px;
      double dx_px = px1 - px0, dy_px = py1 - py0;
      int n_samples =
          (int)std::ceil(std::max(std::fabs(dx_px), std::fabs(dy_px))) + 1;
      // Hottest loop in the engine, so hoist the reciprocal out. IEEE754 does
      // not promise sample*inv == sample/n, so the masks were diffed byte for
      // byte against the regression conditions before this went in.
      const double inv_samples = 1.0 / n_samples;
      for (int sample = 0; sample <= n_samples; ++sample)
      {
        double along = (double)sample * inv_samples;
        int px = (int)std::llround(px0 + along * dx_px);
        int py = (int)std::llround(py0 + along * dy_px);
        if (px >= 0 && px < padded_px && py >= 0 && py < padded_px)
          centreline_px[(std::size_t)py * padded_px + px] = 1;
      }
    }
  };

  // --- per-cell state. Draw order matters, see the file header. ---
  const bool confined = p.plate_scaffold;
  auto on_scaffold = [&](double x, double y) -> bool
  {
    int comp_x, comp_y;
    compartment_of(x, y, field_um, p.scaffold_nx, p.scaffold_ny, comp_x,
                   comp_y);
    return p.scaffold_cells[(std::size_t)comp_y * p.scaffold_nx + comp_x] != 0;
  };

  std::vector<double> soma_x(n_cells), soma_y(n_cells), arbor_cos(n_cells),
      arbor_sin(n_cells), chead(n_cells, 0.0);
  std::vector<std::uint8_t> seeker(n_cells, 0);
  std::vector<int> ctmpl(n_cells, 0);
  for (int n = 0; n < n_cells; ++n)
  {
    double angle = random_angle(rng); // arbor orientation, fixed per cell
    arbor_cos[n] = std::cos(angle);
    arbor_sin[n] = std::sin(angle);
    soma_x[n] = random_pos(rng); // initial soma position
    soma_y[n] = random_pos(rng);
    // Plating: land cells on scaffold, not bare field. Rejection sampling, so
    // the extra draws only happen when confined.
    if (confined)
    {
      // Rejection sampling: about 1/(scaffold fraction) draws per cell, so
      // 2.5 on the shipped 40 % lattice. The bound is a livelock guard, not a
      // tuning knob, at 200k it covers fractions down to about 1e-4.
      const int kPlacementAttempts = 200000;
      int guard = 0;
      while (!on_scaffold(soma_x[n], soma_y[n]))
      {
        soma_x[n] = random_pos(rng);
        soma_y[n] = random_pos(rng);
        if (++guard > kPlacementAttempts)
          throw std::runtime_error(
              "plating: could not place a cell on the scaffold - is the "
              "geometry nearly empty, or seed_frac too small to reach it?");
      }
    }
  }
  // Only drawn when the pool has more than one, so a single-morphology run uses
  // the same draws it always did.
  if (n_templates > 1)
  {
    std::uniform_int_distribution<int> pick(0, n_templates - 1);
    for (int n = 0; n < n_cells; ++n)
      ctmpl[n] = pick(rng);
  }
  for (int n = 0; n < n_cells; ++n)
    ensure_loaded(ctmpl[n]);

  const bool migrating = p.mig_speed_um_per_h > 0.0;
  if (migrating)
  {
    for (int n = 0; n < n_cells; ++n)
      chead[n] = random_angle(rng); // migration heading, evolves over time
    // each cell independently either seeks or migrates normally
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    for (int n = 0; n < n_cells; ++n)
      seeker[n] = (u01(rng) < p.p_seek) ? 1 : 0;
  }

  Result r;
  const double frac_final = div_to_fraction(p.div_days);
  // each cell reveals a fraction of its own arbor, not a shared length
  auto reveal_of = [&](int n, double frac)
  { return frac * pool[(std::size_t)ctmpl[n]].max_path_distance; };
  const bool cumulative = (p.occupancy_mode == OccupancyMode::Cumulative);

  if (!migrating)
  {
    // one stamp per cell where it was seeded; cum and inst are the same here
    for (int n = 0; n < n_cells; ++n)
      rasterize(ctmpl[n], arbor_cos[n], arbor_sin[n], soma_x[n], soma_y[n],
                reveal_of(n, frac_final), true, true, n);
  }
  else
  {
    const double total_hours = p.div_days * 24.0;
    int n_steps = (int)std::ceil(total_hours / p.dt_hours);
    if (n_steps < 1)
      n_steps = 1;
    // Rotational diffusion. The heading turns by Normal(0, sqrt(2*D*dt)) each
    // step and decorrelates over a path length of mig_persist. It uses the base
    // speed even under speed coupling, so persistence stays one number for the
    // whole population.
    const double heading_diffusion =
        p.mig_speed_um_per_h / std::max(1e-9, p.mig_persist_um);
    const double half_field_um = 0.5 * field_um;

    if (cumulative) // step 0 footprint
      for (int n = 0; n < n_cells; ++n)
        rasterize(ctmpl[n], arbor_cos[n], arbor_sin[n], soma_x[n], soma_y[n],
                  reveal_of(n, div_to_fraction(0.0)), true, false, n);

    // Which stiffness map is in effect: largest div <= now, same carry-forward
    // rule as the display. Recomputed per step since a run crosses boundaries.
    const bool grad_active = p.stiff_weight != 0.0 && !p.stiff_maps.empty();
    const bool speed_active =
        p.stiff_speed_weight > 0.0 && !p.stiff_maps.empty();
    const bool stiff_active = grad_active || speed_active;
    const double stiff_comp_w_um = stiff_active ? field_um / p.stiff_nx : 1.0;
    const double stiff_comp_h_um = stiff_active ? field_um / p.stiff_ny : 1.0;

    // Speed-coupling transfer function: beta(S) = 0.2 + 0.8*sigmoid(k*(S-c)),
    // the adhesion-gain form from Elayan 2026 Eq. 20. Their published centre of
    // 15 kPa sits well above this scaffold's 0.1-1.5 kPa, so the NaN defaults
    // calibrate to the loaded maps instead of using constants that flatline.
    double speed_centre = p.speed_centre_kpa;
    double speed_slope = p.speed_slope_per_kpa;
    if (speed_active)
    {
      if (std::isnan(speed_centre) || std::isnan(speed_slope))
      {
        std::vector<double> all_values;
        for (const auto &map : p.stiff_maps)
          all_values.insert(all_values.end(), map.begin(), map.end());
        std::sort(all_values.begin(), all_values.end());
        if (std::isnan(speed_centre))
          speed_centre = percentile_sorted(all_values, 0.5);
        if (std::isnan(speed_slope))
        {
          const double spread = percentile_sorted(all_values, 0.9) -
                                percentile_sorted(all_values, 0.1);
          // -4/spread puts sigmoid(+-2) at the p10/p90 values. A uniform map
          // has no spread, so slope 0 and the coupling is inert, correctly.
          speed_slope = (spread > 0) ? -4.0 / spread : 0.0;
        }
      }
      r.speed_centre_used = speed_centre;
      r.speed_slope_used = speed_slope;
    }
    auto speed_multiplier = [&](double stiffness) -> double
    {
      // Elayan 2026 Eq. 20 again, calibrated above. The sigmoid spans 0.8 on a
      // floor of 0.2, so a cell in the least adhesive region still moves at a
      // fifth of base speed rather than stopping dead.
      const double gain =
          0.2 +
          0.8 / (1.0 + std::exp(-speed_slope * (stiffness - speed_centre)));
      return (1.0 - p.stiff_speed_weight) + p.stiff_speed_weight * gain;
    };

    // Value of the in-effect stiffness map at a position, compartment
    // granularity like cue_bias. No interpolation, no RNG.
    auto stiffness_at = [&](const std::vector<double> &field, double x,
                            double y) -> double
    {
      int comp_x, comp_y;
      compartment_of(x, y, field_um, p.stiff_nx, p.stiff_ny, comp_x, comp_y);
      return field[(std::size_t)comp_y * p.stiff_nx + comp_x];
    };

    // Unit vector of the signed gradient at (x,y), or (0,0) where the field is
    // flat. direction > 0 steers down-gradient, < 0 up. No RNG.
    auto cue_bias = [&](const std::vector<double> &field, int field_nx,
                        int field_ny, double comp_w_um, double comp_h_um,
                        double x, double y, double direction, double &bias_x,
                        double &bias_y)
    {
      bias_x = bias_y = 0.0;
      int comp_x, comp_y;
      compartment_of(x, y, field_um, field_nx, field_ny, comp_x, comp_y);
      const int x_lo = comp_x > 0 ? comp_x - 1 : comp_x;
      const int x_hi = comp_x < field_nx - 1 ? comp_x + 1 : comp_x;
      const int y_lo = comp_y > 0 ? comp_y - 1 : comp_y;
      const int y_hi = comp_y < field_ny - 1 ? comp_y + 1 : comp_y;
      const double grad_x = (field[(std::size_t)comp_y * field_nx + x_hi] -
                             field[(std::size_t)comp_y * field_nx + x_lo]) /
                            ((x_hi - x_lo) * comp_w_um + 1e-12);
      const double grad_y = (field[(std::size_t)y_hi * field_nx + comp_x] -
                             field[(std::size_t)y_lo * field_nx + comp_x]) /
                            ((y_hi - y_lo) * comp_h_um + 1e-12);
      const double grad_len = std::sqrt(grad_x * grad_x + grad_y * grad_y);
      if (grad_len > 1e-12)
      {
        bias_x = -direction * grad_x / grad_len;
        bias_y = -direction * grad_y / grad_len;
      }
    };

    // Scaffold seek has TWO halves, both weighted by seek_weight, picked per
    // step by whether the soma is on the scaffold (see the walk below):
    // OFF it, steer down the distance-to-scaffold field built here; ON it,
    // follow the strand's local tangent from the depth field below. The second
    // half is what spreads cells ALONG a strand, without it every cell piles
    // onto the one it first touched.
    const bool seeking =
        p.seek_mode == SeekMode::Scaffold && p.seek_weight > 0.0;
    const bool has_geometry =
        p.scaffold_nx > 0 && p.scaffold_ny > 0 && !p.scaffold_cells.empty();
    std::vector<double> scaffold_distance_um;
    // Distance to the nearest NON-scaffold compartment: same transform,
    // inverted mask. Its level sets run ALONG a strand, so the perpendicular
    // to its gradient is the tangent guidance follows. Built only when guiding.
    std::vector<double> scaffold_depth_um;
    double geom_comp_w_um = 1.0, geom_comp_h_um = 1.0;
    // Built whenever a geometry is loaded: the soma-scaffold readout needs it
    // even when nothing steers.
    if (has_geometry)
    {
      geom_comp_w_um = field_um / p.scaffold_nx;
      geom_comp_h_um = field_um / p.scaffold_ny;
      scaffold_distance_um = geometry_edt(p.scaffold_cells, p.scaffold_nx,
                                          p.scaffold_ny, field_um, false);
      if (seeking)
        scaffold_depth_um = geometry_edt(p.scaffold_cells, p.scaffold_nx,
                                         p.scaffold_ny, field_um, true);
    }

    // Cell-cell cohesion. Neurons pull themselves into aggregates along their
    // own neurites (Ballesteros-Esteban 2023), so the cue's reach is twice the
    // revealed arbor's span, derived from the reveal, not fitted.
    //
    // Pairwise, not over a field: an EDT per step costs O(grid^2) whatever N is
    // (11 % of a 150x150 run, 46 % at 300x300), this costs O(N^2) and touches
    // no grid, so the compartment budget stays a rasterising question.
    const bool cohering = p.cohere_weight > 0.0 && n_cells > 1;
    std::vector<std::vector<double>> reach_table;
    if (cohering)
    {
      reach_table.resize(pool.size());
      for (std::size_t t = 0; t < pool.size(); ++t)
        if (!pool[t].segments.empty())
          reach_table[t] =
              arbor_reach_table(pool[t].segments, pool[t].max_path_distance);
    }
    std::vector<double> cohere_bx((std::size_t)(cohering ? n_cells : 0), 0.0);
    std::vector<double> cohere_by((std::size_t)(cohering ? n_cells : 0), 0.0);

    for (int step = 1; step <= n_steps; ++step)
    {
      const double hours_elapsed =
          std::min((double)step * p.dt_hours, total_hours);
      const double reveal_frac_now = div_to_fraction(hours_elapsed / 24.0);
      const bool stamp_now = cumulative || (step == n_steps);
      const double heading_step_sd =
          std::sqrt(2.0 * heading_diffusion * p.dt_hours);

      int stiff_map_index = -1;
      if (stiff_active)
      {
        for (std::size_t i = 0; i < p.stiff_divs.size(); ++i)
          if (p.stiff_divs[i] <= hours_elapsed / 24.0)
            stiff_map_index = (int)i;
        // Counted, not assumed: before the earliest map block nothing is in
        // effect, so a late-starting map can leave a "coupled" run almost
        // entirely uncoupled. The CLI warns off these two numbers.
        ++r.migration_steps;
        if (stiff_map_index >= 0)
          ++r.steps_with_stiff_map;
      }

      // Cohesion biases come off the positions at the START of the step, for
      // every cell, before any of them moves. Do this inside the walk loop
      // below and cell 0's new position would steer cell 1, making the result
      // depend on cell order.
      if (cohering)
      {
        for (int n = 0; n < n_cells; ++n)
        {
          const std::vector<double> &tbl = reach_table[(std::size_t)ctmpl[n]];
          double reach = 0.0;
          if (!tbl.empty())
          {
            int bin = (int)std::lround(reveal_frac_now * (kReachBins - 1));
            bin = bin < 0 ? 0 : (bin > kReachBins - 1 ? kReachBins - 1 : bin);
            reach = 2.0 * tbl[(std::size_t)bin];
          }
          double ax = 0.0, ay = 0.0;
          for (int m = 0; reach > 0.0 && m < n_cells; ++m)
          {
            if (m == n)
              continue;
            const double dx = soma_x[m] - soma_x[n];
            const double dy = soma_y[m] - soma_y[n];
            const double d2 = dx * dx + dy * dy;
            if (d2 > reach * reach || d2 < 1e-12)
              continue;
            const double d = std::sqrt(d2);
            // Linear falloff. A neighbour at the edge of reach is held by one
            // thread, one next door by most of the arbor.
            const double w = 1.0 - d / reach;
            ax += w * dx / d;
            ay += w * dy / d;
          }
          const double len = std::sqrt(ax * ax + ay * ay);
          cohere_bx[(std::size_t)n] = len > 1e-12 ? ax / len : 0.0;
          cohere_by[(std::size_t)n] = len > 1e-12 ? ay / len : 0.0;
        }
      }

      for (int n = 0; n < n_cells; ++n)
      {
        chead[n] += heading_step_sd * random_normal(rng);
        double dir_x = std::cos(chead[n]), dir_y = std::sin(chead[n]);

        // The decision rule, both halves: off the scaffold steer toward it, on
        // it follow the structure rather than resuming an isotropic walk. The
        // second half is contact guidance, Li 2025 attributes the M-Alg
        // advantage to topographical cues and reports networks throughout, not
        // cells piled where they first touched, which is what steering at the
        // nearest material alone produces.
        if (seeking && seeker[n])
        {
          const bool on_it = on_scaffold(soma_x[n], soma_y[n]);
          double bias_x = 0.0, bias_y = 0.0;
          if (!on_it)
          {
            cue_bias(scaffold_distance_um, p.scaffold_nx, p.scaffold_ny,
                     geom_comp_w_um, geom_comp_h_um, soma_x[n], soma_y[n], 1.0,
                     bias_x, bias_y);
          }
          else if (!scaffold_depth_um.empty())
          {
            // Depth's gradient points across the strand, so its perpendicular
            // is the tangent. Take the sign that agrees with the heading:
            // guidance aligns a walk, never reverses it. On a flat depth field
            // the gradient vanishes and the cell just walks, guidance never
            // invents a direction.
            double grad_x, grad_y;
            cue_bias(scaffold_depth_um, p.scaffold_nx, p.scaffold_ny,
                     geom_comp_w_um, geom_comp_h_um, soma_x[n], soma_y[n], 1.0,
                     grad_x, grad_y);
            if (grad_x != 0.0 || grad_y != 0.0)
            {
              double tan_x = -grad_y, tan_y = grad_x;
              if (tan_x * dir_x + tan_y * dir_y < 0.0)
              {
                tan_x = -tan_x;
                tan_y = -tan_y;
              }
              bias_x = tan_x;
              bias_y = tan_y;
            }
          }
          if (bias_x != 0.0 || bias_y != 0.0)
          {
            dir_x = (1.0 - p.seek_weight) * dir_x + p.seek_weight * bias_x;
            dir_y = (1.0 - p.seek_weight) * dir_y + p.seek_weight * bias_y;
            double len = std::sqrt(dir_x * dir_x + dir_y * dir_y) + 1e-12;
            dir_x /= len;
            dir_y /= len;
          }
        }

        // Cohesion: steer toward the neighbours this cell can currently
        // reach. Same blend shape as the cue above and no RNG draw, so at
        // weight 0 the run is bit-identical to an uncoupled one.
        if (cohering)
        {
          const double bias_x = cohere_bx[(std::size_t)n];
          const double bias_y = cohere_by[(std::size_t)n];
          // The gate is the competition: substrate adhesion where the cell has
          // found the scaffold, cell-cell adhesion where it has not. Applied to
          // the weight rather than the bias, so a gated-off cell walks exactly
          // as it would with no cohesion at all.
          const double w = (p.cohere_gate == CohereGate::Scaffold &&
                            on_scaffold(soma_x[n], soma_y[n]))
                               ? 0.0
                               : p.cohere_weight;
          if (w > 0.0 && (bias_x != 0.0 || bias_y != 0.0))
          {
            dir_x = (1.0 - w) * dir_x + w * bias_x;
            dir_y = (1.0 - w) * dir_y + w * bias_y;
            double len = std::sqrt(dir_x * dir_x + dir_y * dir_y) + 1e-12;
            dir_x /= len;
            dir_y /= len;
          }
        }

        // Mechanical cue, gradient half. Blended in by |weight|, with the sign
        // picking the direction. No draws, so the RNG stream is unchanged.
        if (grad_active && stiff_map_index >= 0)
        {
          double bias_x, bias_y;
          cue_bias(p.stiff_maps[(std::size_t)stiff_map_index], p.stiff_nx,
                   p.stiff_ny, stiff_comp_w_um, stiff_comp_h_um, soma_x[n],
                   soma_y[n], p.stiff_weight > 0 ? 1.0 : -1.0, bias_x, bias_y);
          if (bias_x != 0.0 || bias_y != 0.0)
          {
            const double weight = std::fabs(p.stiff_weight);
            dir_x = (1.0 - weight) * dir_x + weight * bias_x;
            dir_y = (1.0 - weight) * dir_y + weight * bias_y;
            double len = std::sqrt(dir_x * dir_x + dir_y * dir_y) + 1e-12;
            dir_x /= len;
            dir_y /= len;
          }
        }

        // Mechanical cue, magnitude half. Local stiffness scales the step
        // length. Slow regions collect cells, since density goes as 1/v, which
        // is how a ring-shaped dose field steers with no gradient to follow.
        double step_um = p.mig_speed_um_per_h * p.dt_hours;
        if (speed_active && stiff_map_index >= 0)
          step_um *= speed_multiplier(
              stiffness_at(p.stiff_maps[(std::size_t)stiff_map_index],
                           soma_x[n], soma_y[n]));

        const double prev_x = soma_x[n], prev_y = soma_y[n];
        soma_x[n] += step_um * dir_x;
        soma_y[n] += step_um * dir_y;

        // bounce off the field edge so cells stay on the substrate
        if (soma_x[n] < -half_field_um)
        {
          soma_x[n] = -half_field_um - (soma_x[n] + half_field_um);
          chead[n] = M_PI - chead[n];
        }
        if (soma_x[n] > half_field_um)
        {
          soma_x[n] = half_field_um - (soma_x[n] - half_field_um);
          chead[n] = M_PI - chead[n];
        }
        if (soma_y[n] < -half_field_um)
        {
          soma_y[n] = -half_field_um - (soma_y[n] + half_field_um);
          chead[n] = -chead[n];
        }
        if (soma_y[n] > half_field_um)
        {
          soma_y[n] = half_field_um - (soma_y[n] - half_field_um);
          chead[n] = -chead[n];
        }
        // clamp too, in case one step overshot the far side entirely
        if (soma_x[n] < -half_field_um)
          soma_x[n] = -half_field_um;
        else if (soma_x[n] > half_field_um)
          soma_x[n] = half_field_um;
        if (soma_y[n] < -half_field_um)
          soma_y[n] = -half_field_um;
        else if (soma_y[n] > half_field_um)
          soma_y[n] = half_field_um;

        // Confinement. A step off the scaffold is cancelled, so the cell stays
        // put while its heading keeps diffusing, waiting at the edge until it
        // points inward again. No RNG.
        if (confined && !on_scaffold(soma_x[n], soma_y[n]))
        {
          soma_x[n] = prev_x;
          soma_y[n] = prev_y;
        }

        if (stamp_now)
          rasterize(ctmpl[n], arbor_cos[n], arbor_sin[n], soma_x[n], soma_y[n],
                    reveal_of(n, reveal_frac_now), true, false, n);
      }
    }
    // final frame for the figure: final positions, full DIV reveal
    for (int n = 0; n < n_cells; ++n)
      rasterize(ctmpl[n], arbor_cos[n], arbor_sin[n], soma_x[n], soma_y[n],
                reveal_of(n, frac_final), false, true, n);
  }

  // dilate once over the padded grid, then measure the middle only
  std::vector<std::uint8_t> padded_mask((std::size_t)padded_px * padded_px, 0);
  dilate_by_edt(centreline_px, padded_mask, padded_px, territory_px);

  r.mask.assign((std::size_t)grid_px * grid_px, 0);
  std::size_t covered = 0;
  for (int y = 0; y < grid_px; ++y)
    for (int x = 0; x < grid_px; ++x)
    {
      std::uint8_t covered_here =
          padded_mask[(std::size_t)(y + margin_px) * padded_px +
                      (x + margin_px)];
      r.mask[(std::size_t)y * grid_px + x] = covered_here;
      covered += covered_here;
    }

  r.occupancy_percent = 100.0 * (double)covered / ((double)grid_px * grid_px);
  r.reveal_fraction = frac_final;
  r.max_path_distance = pool[(std::size_t)ctmpl[0]].max_path_distance;
  r.cell_morphology = ctmpl;
  // 0 for templates never sampled, since those were never parsed
  for (const auto &tmpl : pool)
    r.template_max_path.push_back(tmpl.max_path_distance);
  r.grid_side = grid_px;
  r.covered_pixels = covered;
  r.soma_x.assign(soma_x.begin(), soma_x.end());
  r.soma_y.assign(soma_y.begin(), soma_y.end());

  // Where the somata ended up relative to the scaffold: post-processing, no
  // RNG, no number changed. The only readout separating "reached the scaffold"
  // from "spread along it", which occupancy mixes. Reported whenever a geometry
  // is loaded, so an uncoupled run gives a baseline.
  if (!p.scaffold_cells.empty() && p.scaffold_nx > 0 && p.scaffold_ny > 0 &&
      !soma_x.empty())
  {
    const std::vector<double> dist = geometry_edt(
        p.scaffold_cells, p.scaffold_nx, p.scaffold_ny, field_um, false);
    int on_count = 0;
    double dist_sum = 0.0;
    for (std::size_t i = 0; i < soma_x.size(); ++i)
    {
      if (on_scaffold(soma_x[i], soma_y[i]))
        ++on_count;
      int cx, cy;
      compartment_of(soma_x[i], soma_y[i], field_um, p.scaffold_nx,
                     p.scaffold_ny, cx, cy);
      dist_sum += dist[(std::size_t)cy * p.scaffold_nx + cx];
    }
    r.soma_scaffold_count = on_count;
    r.soma_scaffold_mean_dist_um = dist_sum / (double)soma_x.size();
  }

  if (!p.pgm_path.empty())
  {
    std::ofstream f(p.pgm_path, std::ios::binary);
    f << "P5\n" << grid_px << " " << grid_px << "\n255\n";
    std::vector<std::uint8_t> img((std::size_t)grid_px * grid_px);
    for (std::size_t i = 0; i < img.size(); ++i)
      img[i] = r.mask[i] ? 30 : 245;
    f.write((const char *)img.data(), (std::streamsize)img.size());
  }
  return r;
}

} // namespace occupancy
