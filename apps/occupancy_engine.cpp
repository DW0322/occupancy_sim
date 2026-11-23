// CLI for the engine: argument parsing and reporting only, so the behaviour in
// the occupancy library stays testable without argv.
#include "occupancy/clustering.hpp"
#include "occupancy/compartments.hpp"
#include "occupancy/scaffold_geometry.hpp"
#include "occupancy/schedule.hpp"
#include "occupancy/segment.hpp"
#include "occupancy/simulation.hpp"
#include "occupancy/stiffness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using occupancy::Params;

namespace
{
// App-level options. Compartment output is post-processing over Result.mask;
struct CompartmentOpts
{
  int n = 30;           // compartment grid side
  std::string csv_path; // empty = layer off
  // Externally supplied files. Loading either is inert, validated, reported,
  // displayed, but moves no simulated number.
  std::string stiffness_path;       // carry alongside a run
  std::string stiffness_check_path; // validate a file and exit, no run
  std::string geometry_path;
  std::string geometry_check_path;
  // Whether the caller actually named a grid and a field, as opposed to
  // inheriting the defaults below. A --stiffness-check / --geometry-check with
  // neither given is a question about the FILE, so it must not be answered
  // against 30x30/1500 um that nobody asked for: a well-formed 150x150 map -
  // the size INTEGRATION_NOTES tells the producer to build, would come back
  // INVALID for no fault of its own. Given explicitly, the check is against
  // the run, which is the other thing the flag is for.
  bool grid_explicit = false;
  bool field_explicit = false;
  // The switches that make those files actually do something, kept separate
  // from loading them. All default off; that's what keeps the pinned numbers
  // valid, and each changes results once on.
  bool plate_scaffold = false;
  double stiff_weight = 0.0;          // gradient coupling
  double stiff_speed = 0.0;           // magnitude (speed) coupling
  double speed_centre = std::nan(""); // NaN = auto-calibrate from the maps
  double speed_slope = std::nan("");
  bool cohere_check = false; // report reach vs spacing and exit, no run
};

/**
 * @brief Report whether cohesion can act at this density, and print it.
 *
 * @param p Parameters; uses morphology_paths, n_cells, scaffold_um and
 *        seed_fraction.
 * @return 0 if cohesion can act somewhere in DIV 5-25, 1 if inert throughout.
 */
int report_cohesion(const occupancy::Params &p)
{
  if (p.morphology_paths.empty())
  {
    std::cerr << "error: --cohere-check needs --morph or --morph-dir\n";
    return 2;
  }
  if (p.n_cells < 2)
  {
    std::cerr << "error: --cohere-check needs at least 2 cells\n";
    return 2;
  }
  // Mean over the pool: reach varies by template, and a
  // population statistic is what the spacing is.
  std::vector<double> mean_reach((std::size_t)occupancy::kReachBins, 0.0);
  int loaded = 0;
  for (const std::string &path : p.morphology_paths)
  {
    occupancy::MorphologyTemplate t = occupancy::load_template(path);
    if (t.segments.empty())
      continue;
    std::vector<double> tbl =
        occupancy::arbor_reach_table(t.segments, t.max_path_distance);
    for (std::size_t i = 0; i < mean_reach.size(); ++i)
      mean_reach[i] += tbl[i];
    ++loaded;
  }
  if (loaded == 0)
  {
    std::cerr << "error: no morphology in the pool had any segments\n";
    return 2;
  }
  for (double &v : mean_reach)
    v /= (double)loaded;

  // Spacing over the SEEDED region, not the field: dividing by the field
  // overstated it tenfold at --seed-frac 0.1 and reported INERT for a coupling
  // that measurably moves both readouts. Exact at DIV 0, conservative after, so
  // it can under-call reach but never over-call it.
  const double seeded_um = p.seed_fraction * p.scaffold_um;
  const double spacing = seeded_um / std::sqrt((double)p.n_cells);
  const double per_cm2 =
      (double)p.n_cells / ((seeded_um * 1e-4) * (seeded_um * 1e-4));

  std::cout << "COHESION_CHECK field " << p.scaffold_um << " um  seeded "
            << seeded_um << " um (seed_frac " << p.seed_fraction << ")  N "
            << p.n_cells << " (" << per_cm2 << " cells/cm2)  mean spacing "
            << spacing << " um  pool " << loaded << "\n";
  std::cout
      << "  DIV  reveal   arbor_p95      reach   reach/spacing  verdict\n";
  const double divs[] = {5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0};
  bool any_active = false;
  for (double div : divs)
  {
    const double frac = occupancy::div_to_fraction(div);
    int bin = (int)std::lround(frac * (occupancy::kReachBins - 1));
    bin = bin < 0 ? 0
                  : (bin > occupancy::kReachBins - 1 ? occupancy::kReachBins - 1
                                                     : bin);
    const double arbor = mean_reach[(std::size_t)bin];
    const double reach = 2.0 * arbor;
    const double ratio = spacing > 0 ? reach / spacing : 0.0;
    // PROVISIONAL thresholds. Below 1 a cell reaches nobody at mean spacing;
    // above 2 it reaches past its nearest neighbour to the next one out, which
    // is what lets aggregates grow rather than just pair up.
    const char *verdict = ratio < 1.0   ? "INERT - reaches nobody"
                          : ratio < 2.0 ? "MARGINAL - nearest neighbour only"
                                        : "ACTIVE";
    any_active = any_active || ratio >= 1.0;
    std::printf("  %3.0f   %5.2f   %7.1f um  %7.1f um   %8.2f     %s\n", div,
                frac, arbor, reach, ratio, verdict);
  }
  std::cout << "  thresholds PROVISIONAL: calibrated on the shipped pool only. "
               "Raise N, or shrink\n  the field or --seed-frac, to move the "
               "ratio - the reach is fixed by the morphology\n  and the reveal "
               "schedule.\n";
  return any_active ? 0 : 1;
}

void usage(const char *prog)
{
  // One line per flag. Every flag's reasoning, and worked examples, are in
  // docs/RUNNING.md
  // u can just type occsim --help all commands will show as per the
  // documentation.
  std::cerr
      << "usage:\n  " << prog << " --morph FILE [options]\n\n"
      << "  --morph PATH          morphology; repeat to sample from a pool\n"
      << "  --morph-dir DIR       every .asc/.swc/.h5 in DIR, sorted\n"
      << "  --n INT               cells to seed                 [15]\n"
      << "  --field UM            overall field side            [1500]\n"
      << "                        (--scaffold is a legacy alias)\n"
      << "  --territory UM        dilation radius R             [17.2]\n"
      << "  --res UM              um PER PIXEL                  [0.5]\n"
      << "  --div DAYS            days in vitro                 [15]\n"
      << "  --seed INT            RNG seed                      [1]\n"
      << "  --seed-frac F         seeded fraction of side, 0-1  [1.0]\n"
      << "  --mig-speed UM_PER_H  base soma speed               [0]\n"
      << "  --mig-persist UM      walk persistence length       [50]\n"
      << "  --dt HOURS            timestep                      [1.0]\n"
      << "  --occ-mode cum|inst   occupancy definition          [cum]\n"
      << "  --seek none|scaffold  COUPLING: steer toward the scaffold off it,\n"
      << "                        follow the strand on it        [none]\n"
      << "  --seek-weight F       0-1, how hard seek steers      [0.3]\n"
      << "  --p-seek F            0-1, fraction that seek        [1.0]\n"
      << "  --cohere-weight F     COUPLING: 0-1, steer toward reachable\n"
      << "                        neighbours; radius derived     [0]\n"
      << "  --cohere-gate none|scaffold\n"
      << "                        COUPLING: scaffold = only cells OFF it\n"
      << "                        cohere                         [none]\n"
      << "  --cohere-check        report reach vs cell spacing per DIV, then\n"
      << "                        exit (spacing uses --seed-frac)\n"
      << "  --pgm PATH            write coverage image\n"
      << "  --segments PATH       write per-cell segment CSV\n"
      << "  --compartments INT    compartment grid side          [30]\n"
      << "  --compartments-csv P  write per-compartment occupancy CSV\n"
      << "  --stiffness-map P     carry a stiffness map (docs/STIFFNESS_\n"
      << "                        FORMAT.md); inert without a coupling\n"
      << "  --stiffness-check P   validate a stiffness file and exit\n"
      << "  --scaffold-geometry P scaffold mask (docs/SCAFFOLD_GEOMETRY_\n"
      << "                        FORMAT.md); reports scaffold-restricted\n"
      << "                        occupancy alongside, couples nothing\n"
      << "  --geometry-check P    validate a geometry file and exit\n"
      << "  --plate-scaffold      COUPLING: seed on the scaffold and confine\n"
      << "                        migration to it\n"
      << "  --stiff-weight F      COUPLING: stiffness-GRADIENT bias, -1..1;\n"
      << "                        +ve toward soft. Inert on ring fields\n"
      << "  --stiff-speed F       COUPLING: stiffness-MAGNITUDE speed\n"
      << "                        modulation, 0-1. For ring/dose fields\n"
      << "  --speed-centre KPA    sigmoid centre  [auto: map median]\n"
      << "  --speed-slope PER_KPA signed slope; -ve = slower on stiff\n"
      << "                        [auto: -4/(p90-p10) of the maps]\n"
      << "\nFlag reasoning and examples: docs/RUNNING.md\n";
}

// Guarded numeric conversion, just in case
[[noreturn]] void bad_value(const std::string &what, const char *got,
                            const char *expected)
{
  std::cerr << "error: " << what << " expects " << expected << ", got '" << got
            << "'\n";
  std::exit(2);
}

double need_double(const std::string &what, const char *text)
{
  try
  {
    std::size_t used = 0;
    const double v = std::stod(text, &used);
    // otherwise "15abc" parses as 15 and quietly drops the rest
    if (text[used] != '\0')
      bad_value(what, text, "a number");
    return v;
  }
  catch (const std::exception &)
  {
    bad_value(what, text, "a number");
  }
}

int need_int(const std::string &what, const char *text)
{
  try
  {
    std::size_t used = 0;
    const int v = std::stoi(text, &used);
    if (text[used] != '\0')
      bad_value(what, text, "a whole number");
    return v;
  }
  catch (const std::exception &)
  {
    bad_value(what, text, "a whole number");
  }
}

// The enum flags parse strictly, so a typo throws instead of quietly selecting
// the default. parse_named() runs before main()'s try block, so catch here and
// report it the same way a bad number is reported.
template <typename F>
auto need_enum(const std::string &what, const char *text, const char *expected,
               F parse) -> decltype(parse(std::string()))
{
  try
  {
    return parse(std::string(text));
  }
  catch (const std::exception &)
  {
    bad_value(what, text, expected);
  }
}

unsigned need_unsigned(const std::string &what, const char *text)
{
  // refuse the sign outright. i cant be dealing with c++ unsigned behaviour and
  // loop arounds
  if (text[0] == '-' || text[0] == '+')
    bad_value(what, text, "a non-negative whole number");
  try
  {
    std::size_t used = 0;
    const unsigned long v = std::stoul(text, &used);
    if (text[used] != '\0')
      bad_value(what, text, "a non-negative whole number");
    return (unsigned)v;
  }
  catch (const std::exception &)
  {
    bad_value(what, text, "a non-negative whole number");
  }
}

/**
 * @brief Parse the CLI flags, can never have enough checks, I learnt that the
 * hard way
 * @param argc Argument count from main().
 * @param argv Argument vector from main().
 * @param p Output: simulation parameters, filled in per recognised flag.
 * @param co Output: compartment/stiffness/geometry options.
 * @return false on an unknown flag, in which case print usage() and exit.
 */

// YES, ITS A LOT OF IF STATEMENTS I KNOW
bool parse_named(int argc, char **argv, Params &p, CompartmentOpts &co)
{
  auto need = [&](int i) -> const char *
  {
    if (i + 1 >= argc)
    {
      std::cerr << "error: " << argv[i] << " needs a value\n";
      std::exit(2);
    }
    return argv[i + 1];
  };
  for (int i = 1; i < argc; ++i)
  {
    std::string a = argv[i];
    if (a == "--morph")
      p.morphology_paths.push_back(need(i++));
    else if (a == "--morph-dir")
    {
      // Every morphology in a directory, sorted. Order matters, the sampled
      // index maps into the pool in this order, so sorting is what keeps runs
      // reproducible across machines and filesystems.
      namespace fs = std::filesystem;
      std::string dir = need(i++);
      std::vector<std::string> found;
      std::error_code ec;
      for (const auto &e : fs::directory_iterator(dir, ec))
      {
        std::string ext = e.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".asc" || ext == ".swc" || ext == ".h5")
          found.push_back(e.path().string());
      }
      if (ec || found.empty())
      {
        std::cerr << "error: no .asc/.swc/.h5 morphologies in " << dir << "\n";
        std::exit(2);
      }
      std::sort(found.begin(), found.end());
      p.morphology_paths.insert(p.morphology_paths.end(), found.begin(),
                                found.end());
    }
    else if (a == "--n")
      p.n_cells = need_int(a, need(i++));
    else if (a == "--field" || a == "--scaffold")
    {
      p.scaffold_um = need_double(a, need(i++));
      co.field_explicit = true;
    }
    else if (a == "--territory")
      p.territory_um = need_double(a, need(i++));
    else if (a == "--res")
      p.resolution_um = need_double(a, need(i++));
    else if (a == "--div")
      p.div_days = need_double(a, need(i++));
    else if (a == "--seed")
      p.seed = need_unsigned(a, need(i++));
    else if (a == "--seed-frac")
      p.seed_fraction = need_double(a, need(i++));
    else if (a == "--mig-speed")
      p.mig_speed_um_per_h = need_double(a, need(i++));
    else if (a == "--mig-persist")
      p.mig_persist_um = need_double(a, need(i++));
    else if (a == "--dt")
      p.dt_hours = need_double(a, need(i++));
    else if (a == "--occ-mode")
      p.occupancy_mode = need_enum(a, need(i++), "cum or inst",
                                   occupancy::parse_occupancy_mode);
    else if (a == "--seek")
      p.seek_mode = need_enum(a, need(i++), "none or scaffold",
                              occupancy::parse_seek_mode);
    else if (a == "--seek-weight")
      p.seek_weight = need_double(a, need(i++));
    else if (a == "--p-seek")
      p.p_seek = need_double(a, need(i++));
    else if (a == "--cohere-weight")
      p.cohere_weight = need_double(a, need(i++));
    else if (a == "--cohere-gate")
      p.cohere_gate = need_enum(a, need(i++), "none or scaffold",
                                occupancy::parse_cohere_gate);
    else if (a == "--cohere-check")
      co.cohere_check = true;
    else if (a == "--pgm")
      p.pgm_path = need(i++);
    else if (a == "--segments")
      p.segments_path = need(i++);
    else if (a == "--compartments")
    {
      co.n = need_int(a, need(i++));
      co.grid_explicit = true;
    }
    else if (a == "--compartments-csv")
      co.csv_path = need(i++);
    else if (a == "--stiffness-map")
      co.stiffness_path = need(i++);
    else if (a == "--stiffness-check")
      co.stiffness_check_path = need(i++);
    else if (a == "--scaffold-geometry")
      co.geometry_path = need(i++);
    else if (a == "--geometry-check")
      co.geometry_check_path = need(i++);
    else if (a == "--plate-scaffold")
      co.plate_scaffold = true;
    else if (a == "--stiff-weight")
      co.stiff_weight = need_double(a, need(i++));
    else if (a == "--stiff-speed")
      co.stiff_speed = need_double(a, need(i++));
    else if (a == "--speed-centre")
      co.speed_centre = need_double(a, need(i++));
    else if (a == "--speed-slope")
      co.speed_slope = need_double(a, need(i++));
    else if (a == "-h" || a == "--help")
    {
      usage(argv[0]);
      std::exit(0);
    }
    else
    {
      std::cerr << "error: unknown option " << a << "\n";
      return false;
    }
  }
  return true;
}
} // namespace

namespace
{
// How each report line ends: whether the file just loaded actually reached the
// simulation. Read off the coupling flags, so a run cannot claim a disposition
// it does not have. Like every line except OCCUPANCY, these avoid the
// "= <num> %" pattern the tooling regex matches.
const char *kCarried = " VALID carried-not-coupled\n";
const char *kCheckOnly = " VALID file-check-only, no run\n";

/**
 * @brief Print the one-line stiffness report for an already-loaded map.
 * @param m The loaded stiffness map.
 * @param path Source file path, printed for reference.
 * @param grid_n This run's compartment grid side, for the mismatch check.
 * @param scaffold_um This run's scaffold side, um, for the mismatch check.
 * @param run_div The DIV to report the in-effect map for.
 * @param disposition Trailing text saying whether the map was coupled in. Must
 *        not contain a "= <num> %" pattern.
 * @return The mismatch message, "" when usable.
 */
std::string report_stiffness(const occupancy::StiffnessMap &m,
                             const std::string &path, int grid_n,
                             double scaffold_um, double run_div,
                             const std::string &disposition)
{
  std::string mismatch =
      occupancy::validate_stiffness_against(m, grid_n, scaffold_um);
  std::cout << "STIFFNESS " << path << " " << m.version << " grid " << m.nx
            << "x" << m.ny << " units " << m.units << " divs [";
  for (std::size_t i = 0; i < m.divs.size(); ++i)
    std::cout << (i ? "," : "") << m.divs[i];
  std::cout << "]";
  if (!mismatch.empty())
  {
    std::cout << " INVALID: " << mismatch << "\n";
    return mismatch;
  }
  const int sel = m.select(run_div);
  if (sel < 0)
    std::cout << " no map at or before run DIV " << run_div
              << " (earliest map: DIV " << m.divs.front() << ")";
  else
  {
    double lo = m.maps[(std::size_t)sel][0], hi = lo;
    for (double v : m.maps[(std::size_t)sel])
    {
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
    std::cout << " selected map DIV " << m.divs[(std::size_t)sel]
              << " for run DIV " << run_div << " range [" << lo << ", " << hi
              << "] " << m.units;
  }
  std::cout << disposition;
  return "";
}

/**
 * @brief Print the one-line report for an already-loaded scaffold geometry.
 * @param g The loaded geometry.
 * @param path Source file path, printed for reference.
 * @param grid_n This run's compartment grid side, for the mismatch check.
 * @param scaffold_um This run's scaffold side, um, for the mismatch check.
 * @param disposition Trailing text saying whether the geometry was coupled in.
 *        Must not contain a "= <num> %" pattern.
 * @return The mismatch message, "" when usable.
 */
std::string report_geometry(const occupancy::ScaffoldGeometry &g,
                            const std::string &path, int grid_n,
                            double scaffold_um, const std::string &disposition)
{
  std::string mismatch =
      occupancy::validate_geometry_against(g, grid_n, scaffold_um);
  std::cout << "GEOMETRY " << path << " " << g.version << " grid " << g.nx
            << "x" << g.ny << " scaffold compartments "
            << g.scaffold_compartments() << " of " << g.is_scaffold.size()
            << " (" << (100.0 * g.scaffold_fraction()) << " pct of field)";
  if (!mismatch.empty())
  {
    std::cout << " INVALID: " << mismatch << "\n";
    return mismatch;
  }
  std::cout << disposition;
  return "";
}
} // namespace

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    usage(argv[0]);
    return 2;
  }
  Params p;
  CompartmentOpts co;
  if (!parse_named(argc, argv, p, co))
  {
    usage(argv[0]);
    return 2;
  }

  if (p.segments_path.empty())
    if (const char *env = std::getenv("OCC_SEGMENTS_CSV"))
      p.segments_path = env;

  // --stiffness-check: validate a file and exit. Needs no morphology and runs
  // no simulation, so it goes before Params validation.
  if (!co.stiffness_check_path.empty())
  {
    try
    {
      occupancy::StiffnessMap m =
          occupancy::load_stiffness_map(co.stiffness_check_path);
      // Compare against what the caller named, and against the file's own
      // declaration for whatever they did not, so the bare form is a pure
      // format check that any well-formed map passes at its own size.
      const int grid_n = co.grid_explicit ? co.n : m.nx;
      const double field_um = co.field_explicit ? p.scaffold_um : m.scaffold_um;
      std::string mismatch =
          report_stiffness(m, co.stiffness_check_path, grid_n, field_um,
                           p.div_days, kCheckOnly);
      return mismatch.empty() ? 0 : 1;
    }
    catch (const std::exception &e)
    {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
    }
  }

  if (!co.geometry_check_path.empty())
  {
    try
    {
      occupancy::ScaffoldGeometry g =
          occupancy::load_scaffold_geometry(co.geometry_check_path);
      // Same rule as --stiffness-check above.
      const int grid_n = co.grid_explicit ? co.n : g.nx;
      const double field_um = co.field_explicit ? p.scaffold_um : g.scaffold_um;
      std::string mismatch = report_geometry(g, co.geometry_check_path, grid_n,
                                             field_um, kCheckOnly);
      return mismatch.empty() ? 0 : 1;
    }
    catch (const std::exception &e)
    {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
    }
  }

  // Load and validate before the run, so a bad file fails in milliseconds
  // instead of after the whole simulation. Reported after the OCCUPANCY line.
  // this saved me 34 minutes of run time :,)
  occupancy::StiffnessMap stiff;
  bool have_stiff = false;
  if (!co.stiffness_path.empty())
  {
    try
    {
      stiff = occupancy::load_stiffness_map(co.stiffness_path);
      std::string mismatch =
          occupancy::validate_stiffness_against(stiff, co.n, p.scaffold_um);
      if (!mismatch.empty())
      {
        std::cerr << "error: stiffness map " << co.stiffness_path << ": "
                  << mismatch << "\n";
        return 2;
      }
      have_stiff = true;
    }
    catch (const std::exception &e)
    {
      std::cerr << "error: " << e.what() << "\n";
      return 2;
    }
  }

  // Same deal for the geometry file: carried and reported, only couples to the
  // dynamics under --plate-scaffold or --seek scaffold.
  occupancy::ScaffoldGeometry geom;
  bool have_geom = false;
  if (!co.geometry_path.empty())
  {
    try
    {
      geom = occupancy::load_scaffold_geometry(co.geometry_path);
      std::string mismatch =
          occupancy::validate_geometry_against(geom, co.n, p.scaffold_um);
      if (!mismatch.empty())
      {
        std::cerr << "error: scaffold geometry " << co.geometry_path << ": "
                  << mismatch << "\n";
        return 2;
      }
      have_geom = true;
    }
    catch (const std::exception &e)
    {
      std::cerr << "error: " << e.what() << "\n";
      return 2;
    }
  }

  if (co.plate_scaffold && !have_geom)
  {
    std::cerr << "error: --plate-scaffold needs --scaffold-geometry\n";
    return 2;
  }
  if (p.seek_mode == occupancy::SeekMode::Scaffold && !have_geom)
  {
    std::cerr << "error: --seek scaffold needs --scaffold-geometry\n";
    return 2;
  }
  if (have_geom)
  {
    // The geometry data always rides along once loaded; the two switches that
    // read it are plate_scaffold and seek_mode, and validate() plus run() key
    // on those, never on the data being present.
    p.scaffold_cells = geom.is_scaffold;
    p.scaffold_nx = geom.nx;
    p.scaffold_ny = geom.ny;
    p.plate_scaffold = co.plate_scaffold;
  }

  if ((co.stiff_weight != 0.0 || co.stiff_speed != 0.0) && !have_stiff)
  {
    std::cerr << "error: --stiff-weight/--stiff-speed need --stiffness-map\n";
    return 2;
  }
  if (have_stiff)
  {
    p.stiff_weight = co.stiff_weight;
    p.stiff_speed_weight = co.stiff_speed;
    p.speed_centre_kpa = co.speed_centre;
    p.speed_slope_per_kpa = co.speed_slope;
    if (co.stiff_weight != 0.0 || co.stiff_speed != 0.0)
    {
      p.stiff_divs = stiff.divs;
      p.stiff_maps = stiff.maps;
      p.stiff_nx = stiff.nx;
      p.stiff_ny = stiff.ny;
    }
  }

  std::string bad = occupancy::validate(p);
  if (!bad.empty())
  {
    std::cerr << "error: " << bad << "\n";
    return 2;
  }

  // Screen and exit. No simulation, so it is cheap enough to run before every
  // cohesion sweep rather than after one that turns out to be pointless
  if (co.cohere_check)
    return report_cohesion(p);

  // What each carried file did to this run. Read off the same p.* fields the
  // simulation uses, so a line cannot claim something the run did not do.
  // validate() has already rejected couplings that could never act.
  std::ostringstream stiff_disposition, geom_disposition;
  if (p.stiff_weight != 0.0 || p.stiff_speed_weight != 0.0)
  {
    stiff_disposition << " VALID COUPLED via";
    if (p.stiff_weight != 0.0)
      stiff_disposition << " --stiff-weight " << p.stiff_weight
                        << " (gradient, "
                        << (p.stiff_weight > 0 ? "toward soft" : "toward stiff")
                        << ", direction contested - see INTEGRATION_NOTES)";
    if (p.stiff_speed_weight != 0.0)
      stiff_disposition << " --stiff-speed " << p.stiff_speed_weight
                        << " (magnitude modulates step length)";
    stiff_disposition << "\n";
  }
  else
    stiff_disposition << kCarried;
  if (p.plate_scaffold || p.seek_mode == occupancy::SeekMode::Scaffold)
  {
    geom_disposition << " VALID COUPLED via";
    if (p.plate_scaffold)
      geom_disposition << " --plate-scaffold (cells seeded on and confined to "
                          "the scaffold)";
    if (p.seek_mode == occupancy::SeekMode::Scaffold)
      geom_disposition << " --seek scaffold (off-scaffold cells steer toward "
                          "it, weight "
                       << p.seek_weight << ")";
    geom_disposition << "\n";
  }
  else
    geom_disposition << kCarried;

  if (p.div_days > occupancy::reveal_saturation_div() &&
      p.mig_speed_um_per_h <= 0.0)
    std::cerr << "note: DIV " << p.div_days
              << " with migration off. The reveal "
                 "schedule saturates at DIV "
              << occupancy::reveal_saturation_div()
              << ", so this is identical to DIV "
              << occupancy::reveal_saturation_div() << ".\n";

  try
  {
    occupancy::Result r = occupancy::run(p);
    std::cout << "OCCUPANCY N=" << p.n_cells << " L=" << p.scaffold_um
              << " R=" << p.territory_um << " res=" << p.resolution_um
              << " DIV=" << p.div_days << " frac=" << r.reveal_fraction
              << " seed=" << p.seed << " seed_frac=" << p.seed_fraction
              << " mig_speed=" << p.mig_speed_um_per_h
              << " mig_persist=" << p.mig_persist_um << " dt=" << p.dt_hours
              << " mode=" << occupancy::occupancy_mode_name(p.occupancy_mode)
              << " seek=" << occupancy::seek_mode_name(p.seek_mode)
              << " seek_w=" << p.seek_weight << " p_seek=" << p.p_seek
              << " stiff_weight=" << p.stiff_weight
              << " stiff_speed=" << p.stiff_speed_weight
              << " plate_scaffold=" << (p.plate_scaffold ? "on" : "off")
              << " cohere_w=" << p.cohere_weight
              << " cohere_gate=" << occupancy::cohere_gate_name(p.cohere_gate)
              << " = " << r.occupancy_percent << " %\n";

    if (r.migration_steps > 0 && r.steps_with_stiff_map < r.migration_steps)
    {
      const double active_pct =
          100.0 * (double)r.steps_with_stiff_map / (double)r.migration_steps;
      if (r.steps_with_stiff_map == 0)
        std::cerr << "warning: stiffness coupling never acted - no map in "
                     "effect in "
                  << r.migration_steps << " timesteps (earliest map DIV "
                  << stiff.divs.front() << ", run reaches DIV " << p.div_days
                  << "), so this equals an uncoupled run.\n";
      else if (active_pct < 50.0)
        std::cerr << "warning: stiffness coupling acted on only "
                  << r.steps_with_stiff_map << " of " << r.migration_steps
                  << " timesteps (" << active_pct << " %); no map before DIV "
                  << stiff.divs.front()
                  << ", so a weak effect here is inconclusive, not a null.\n";
    }

    // The transfer function the speed coupling resolved to. The NaN defaults
    // calibrate off the maps, so a logged weight on its own would not reproduce
    // the run. Off the "= <num> %" pattern like every other line.
    if (p.stiff_speed_weight != 0.0)
      std::cout << "STIFFNESS_SPEED centre " << r.speed_centre_used
                << " kPa slope " << r.speed_slope_used
                << " per kPa (beta = 0.2 + 0.8*sigmoid(slope*(S-centre)); "
                << (r.speed_slope_used < 0
                        ? "slower on stiff - cells accumulate in stiff regions"
                    : r.speed_slope_used > 0
                        ? "faster on stiff - cells accumulate in soft "
                          "regions"
                        : "slope 0 - the maps are uniform and the "
                          "coupling is inert")
                << ")\n";

    double reported_pct = r.occupancy_percent;
    const bool need_bins = !co.csv_path.empty() || have_geom;
    if (need_bins && co.n < 1)
    {
      std::cerr << "error: --compartments must be >= 1\n";
      return 2;
    }
    occupancy::CompartmentGrid g;
    if (need_bins)
    {
      g = occupancy::bin_mask(r.mask, r.grid_side, co.n);
      g.resolution_um = p.resolution_um;
      g.scaffold_um = p.scaffold_um;
    }

    if (!co.csv_path.empty())
    {
      // Every CSV carries its own invocation: argv verbatim, the resolved
      // geometry, and the aggregation check.
      std::ostringstream inv;
      for (int i = 0; i < argc; ++i)
        inv << (i ? " " : "") << argv[i];
      std::ostringstream grid_desc;
      grid_desc << "grid " << g.nx << "x" << g.ny << " compartments over "
                << r.grid_side << "x" << r.grid_side << " px, "
                << (p.scaffold_um / g.nx) << " um per compartment side; "
                << "cy=0 is the y=-L/2 edge";
      std::ostringstream agg;
      agg << "aggregate check: sum(covered)=" << g.covered_sum()
          << " covered_pixels=" << r.covered_pixels << " global "
          << r.occupancy_percent << " pct, aggregate " << g.aggregate_percent()
          << " pct, "
          << (g.aggregate_percent() == r.occupancy_percent &&
                      g.covered_sum() == r.covered_pixels
                  ? "EXACT"
                  : "MISMATCH");
      occupancy::write_compartments_csv(
          co.csv_path, g,
          {inv.str(), "occupancy-sim compartments v1", grid_desc.str(),
           agg.str(), "DIV " + std::to_string(p.div_days)});

      const bool exact = g.aggregate_percent() == r.occupancy_percent &&
                         g.covered_sum() == r.covered_pixels;
      std::cout << "COMPARTMENTS " << g.nx << "x" << g.ny << " -> "
                << co.csv_path
                << " aggregate_check=" << (exact ? "EXACT" : "MISMATCH")
                << "\n";
      if (!exact)
        return 1; // an inexact aggregate is a bug
    }

    if (have_stiff)
    {
      report_stiffness(stiff, co.stiffness_path, co.n, p.scaffold_um,
                       p.div_days, stiff_disposition.str());
      const int sel = stiff.select(p.div_days);
      if (sel >= 0 && !r.soma_x.empty())
      {
        const auto &map_in_effect = stiff.maps[(std::size_t)sel];
        double sum = 0;
        for (std::size_t i = 0; i < r.soma_x.size(); ++i)
        {
          int cx, cy;
          occupancy::compartment_of(r.soma_x[i], r.soma_y[i], p.scaffold_um,
                                    stiff.nx, stiff.ny, cx, cy);
          sum += map_in_effect[(std::size_t)cy * stiff.nx + cx];
        }
        std::cout << "SOMA_STIFFNESS mean " << sum / (double)r.soma_x.size()
                  << " " << stiff.units << " over " << r.soma_x.size()
                  << " somata (map DIV " << stiff.divs[(std::size_t)sel]
                  << ")\n";
      }
    }

    if (have_geom)
    {
      report_geometry(geom, co.geometry_path, co.n, p.scaffold_um,
                      geom_disposition.str());
      std::size_t cov = 0, tot = 0;
      for (std::size_t i = 0; i < geom.is_scaffold.size(); ++i)
        if (geom.is_scaffold[i])
        {
          cov += g.covered[i];
          tot += g.total[i];
        }
      const double restricted = tot ? 100.0 * (double)cov / (double)tot : 0.0;
      reported_pct = restricted;
      std::cout << "OCCUPANCY_SCAFFOLD_RESTRICTED " << restricted
                << " pct over " << geom.scaffold_compartments()
                << " scaffold compartments (whole-field headline above: "
                << r.occupancy_percent
                << " pct; denominator quantised to compartments)\n";
    }

    // Where the somata sit relative to the scaffold: the closest thing to what
    // Li 2025 reports, and the readout that shows whether a steering rule
    // worked. Occupancy mixes reaching the scaffold with spreading along it, so
    // piling cells on one strand and spreading them over the lattice can report
    // the same coverage.
    if (r.soma_scaffold_count >= 0)
    {
      const double frac =
          (double)r.soma_scaffold_count / (double)r.soma_x.size();
      std::cout << "SOMA_SCAFFOLD on_scaffold " << 100.0 * frac << " pct ("
                << r.soma_scaffold_count << " of " << r.soma_x.size()
                << " somata), mean distance to scaffold "
                << r.soma_scaffold_mean_dist_um << " um\n";
    }

    {
      // Measured against the region the cells were actually allowed into, NOT
      // the whole field. With --plate-scaffold that is the scaffold; otherwise
      // it is the convex hull of the somata, because a run that seeds a
      // sub-field cluster (--seed-frac < 1) leaves them nowhere near filling
      // the field. Either way, using the field instead would score the
      // confinement, or the seeding, as clustering: a 300 um cluster in a
      // 3000 um field overstates the region 100-fold, which on its own pins the
      // statistic below 1 whatever the cells did.
      double area = p.scaffold_um * p.scaffold_um;
      double perimeter = 4.0 * p.scaffold_um;
      const char *over = "field";
      const occupancy::Region hull = occupancy::soma_hull(r.soma_x, r.soma_y);
      if (hull.area_um2 > 0.0 && !(have_geom && co.plate_scaffold))
      {
        area = hull.area_um2;
        perimeter = hull.perimeter_um;
        over = "soma hull";
      }
      if (have_geom && co.plate_scaffold)
      {
        const double cw = p.scaffold_um / (double)geom.nx;
        const double ch = p.scaffold_um / (double)geom.ny;
        area = (double)geom.scaffold_compartments() * cw * ch;
        // Perimeter of the quantised mask: every edge where a scaffold
        // compartment meets a non-scaffold one or the field boundary.
        std::size_t edges = 0;
        for (int gy = 0; gy < geom.ny; ++gy)
          for (int gx = 0; gx < geom.nx; ++gx)
          {
            if (!geom.is_scaffold[(std::size_t)gy * geom.nx + gx])
              continue;
            const int dx[4] = {1, -1, 0, 0}, dy[4] = {0, 0, 1, -1};
            for (int k = 0; k < 4; ++k)
            {
              const int nx2 = gx + dx[k], ny2 = gy + dy[k];
              if (nx2 < 0 || nx2 >= geom.nx || ny2 < 0 || ny2 >= geom.ny ||
                  !geom.is_scaffold[(std::size_t)ny2 * geom.nx + nx2])
                ++edges;
            }
          }
        perimeter = (double)edges * 0.5 * (cw + ch);
        over = "scaffold";
      }
      const occupancy::ClusteringStats cs =
          occupancy::clustering_stats(r.soma_x, r.soma_y, area, perimeter);
      if (cs.n >= 2)
      {
        std::cout << "CLUSTERING clark_evans " << cs.clark_evans << " (mean_nn "
                  << cs.mean_nn_um << " um vs " << cs.expected_nn_um
                  << " um expected at random over the " << over << ", n "
                  << cs.n << "); <1 clumped, >1 spread\n";
        // A hull is drawn BY the points, so it undercovers the region they
        // came from and the null sits above 1, badly so at small n. Printing
        // it beside the statistic is the difference between a readable number
        // and one that gets quoted against 1.0. Measured over 2000 uniform
        // trials; docs/RUNNING.md section 2 carries the full row.
        if (std::string(over) == "soma hull")
        {
          const int kN[] = {3, 6, 10, 15, 25, 40, 60, 120};
          const double kNull[] = {3.39, 1.95, 1.54, 1.35,
                                  1.23, 1.15, 1.10, 1.06};
          int best = 0;
          for (int i = 1; i < 8; ++i)
            if (std::abs(kN[i] - cs.n) < std::abs(kN[best] - cs.n))
              best = i;
          std::cout << "  note: over a hull the random null is above 1 - about "
                    << kNull[best] << " at n " << kN[best]
                    << ". Read this against that, not against 1.0.\n";
        }
      }
    }

    // Occupancy stops discriminating well before 100 %: uniformly seeded on
    // the old 1500 um field, scaffold coverage runs 57.9 % at N=3, 98.1 % at
    // N=15, 100.0 % at N=60. Past ~95 % a difference between conditions means
    // almost nothing, so warn rather than let it be quoted as if it did. The
    // committed set avoids the ceiling by seeding N=25 as a central cluster in
    // a 3000 um field, the initial condition was the cause, not the count.
    if (reported_pct >= 95.0)
      std::cerr << "warning: occupancy " << reported_pct
                << " % is saturated - past ~95 % it barely responds to N, DIV "
                   "or any coupling. Lower --n, or read CLUSTERING for an "
                   "arrangement question.\n";
  }
  catch (const std::exception &e)
  {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
