#ifndef OCCUPANCY_SIMULATION_HPP
#define OCCUPANCY_SIMULATION_HPP
// Params, Result, and run(): the entry point for one simulation.
#include "occupancy/migration.hpp"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace occupancy
{

/// @brief How occupancy gets counted once somata move.
enum class OccupancyMode
{
  Cumulative, ///< Union of everywhere reached by DIV t. Only goes up. Default.
  Instant     ///< Footprint at DIV t only, so vacated area stops counting.
};

/**
 * @brief Parse an occupancy mode from the CLI.
 *
 * Strict: an unrecognised string is rejected, not defaulted.
 *
 * @param s "cum" or "inst".
 * @return The mode.
 * @throws std::invalid_argument If s is neither.
 */
OccupancyMode parse_occupancy_mode(const std::string &s);

/**
 * @brief Inverse of parse_occupancy_mode().
 * @param m The mode.
 * @return "inst" or "Cumulative".
 */
const char *occupancy_mode_name(OccupancyMode m);

/**
 * @brief Every knob the engine takes.
 *
 * Changing a default breaks the pins in tests/test_regression.cpp. Anything
 * with no literature source behind it is tagged [FREE PARAMETER] below.
 */
struct Params
{
  std::vector<std::string> morphology_paths;
  int n_cells = 15;            ///< Number of cells to seed.
  double scaffold_um = 1500.0; ///< Field side, um. Its area is the denominator.
  double territory_um = 17.2;
  double resolution_um = 0.5; ///< Raster resolution, um per pixel.
  double div_days = 15.0;     ///< Days in vitro to simulate to.
  unsigned seed = 1;          ///< RNG seed.
  double seed_fraction =
      1.0; ///< Fraction of the side cells seed across; 1.0 = whole field.

  /**
   * Soma speed, um/h, same for every cell; only stiff_speed_weight changes it
   * locally. 0 turns migration off without touching the RNG, so such a run
   * still matches the pre-migration engine exactly. Results use 10, inside the
   * literature's 6-60 band, though unmeasured for this system.
   */
  double mig_speed_um_per_h = 0.0;
  double mig_persist_um = 50.0;
  double dt_hours = 1.0; ///< Timestep, hours.
  OccupancyMode occupancy_mode = OccupancyMode::Cumulative;
  SeekMode seek_mode = SeekMode::None; ///< See SeekMode.
  double seek_weight = 0.3;
  double p_seek = 1.0;
  double cohere_weight = 0.0;

  /**
   * What decides whether a cell coheres. See CohereGate. Scaffold makes
   * cohesion the complement of seek, reaggregate where adhesion is weak,
   * spread where it is strong (Ballesteros-Esteban 2023), so one adhesion axis
   * gives both behaviours. Needs scaffold_cells. Default None, since carrying a
   * geometry must not change a run by itself.
   */
  CohereGate cohere_gate = CohereGate::None;

  std::string pgm_path; ///< Optional output path for a coverage image.
  std::string
      segments_path; ///< Optional output path for a per-cell segment CSV.

  /// Which compartments are scaffold, row-major scaffold_nx * scaffold_ny.
  /// Empty means no geometry. Holding it changes nothing: the switches that
  /// read it are plate_scaffold and seek_mode, both off by default.
  std::vector<std::uint8_t> scaffold_cells; ///< Row-major ny*nx, 0/1.
  int scaffold_nx = 0, scaffold_ny = 0; ///< Grid dimensions of scaffold_cells.

  /// Seed onto the scaffold and keep migration inside it. Needs
  /// scaffold_cells. A step that leaves is cancelled while the heading keeps
  /// diffusing, so the cell waits at the edge until it points back in. Li 2025:
  /// off-scaffold adhesion stays transient.
  bool plate_scaffold = false;

  /**
   * Stiffness-GRADIENT migration bias, -1..1. 0 is off. Magnitude blends
   * against the walk like seek_weight; sign is the direction, positive being
   * down-gradient toward soft (Koser 2016, though durotaxis direction is
   * contested, flipping the sign is the whole cost of testing the other way).
   *
   * Needs a large-scale gradient. On a ring or interference field, a Bessel
   * dose map, a cell drops into the nearest local minimum and this is inert;
   * stiff_speed_weight is the mechanism there.
   */
  double stiff_weight = 0.0;

  /**
   * Stiffness-MAGNITUDE speed coupling, 0..1. 0 is off. Scales the step by
   * (1 - w) + w * beta(S), beta the Elayan 2026 adhesion-gain sigmoid (Eq. 20):
   * beta(S) = 0.2 + 0.8 * sigmoid(slope * (S - centre)).
   *
   * A walker with space-dependent speed piles up where it is slow, since
   * steady-state density goes as 1/v (Schnitzer 1993), so this reads local
   * magnitude and works on rings the gradient coupling cannot see. No gradient,
   * no RNG, so w = 0 is bit-identical to uncoupled.
   */
  double stiff_speed_weight = 0.0;

  /// Sigmoid centre, map units (expect kPa). NaN calibrates to the median over
  /// every loaded map. The published constant sits at 15 kPa and would flatline
  /// on a 0.1-1.5 kPa scaffold.
  double speed_centre_kpa = std::nan("");

  /// Signed sigmoid slope, 1/units. NaN calibrates to -4 / (p90 - p10) over the
  /// loaded maps, so beta covers most of its range across the map's spread.
  /// Negative = slower on stiff, so cells build up there. Positive flips it.
  double speed_slope_per_kpa = std::nan("");

  std::vector<double> stiff_divs; ///< Ascending DIVs the maps apply from.
  std::vector<std::vector<double>> stiff_maps; ///< Row-major ny*nx per div.
  int stiff_nx = 0, stiff_ny = 0;              ///< Grid dimensions of each map.
};

/// @brief What one run produces.
struct Result
{
  double occupancy_percent =
      0.0; ///< The headline number: covered area / field area, %.
  double reveal_fraction = 0.0; ///< div_to_fraction(div_days).
  double max_path_distance =
      0.0;           ///< max_path_distance of the first sampled template.
  int grid_side = 0; ///< Side of mask, pixels.
  std::size_t covered_pixels = 0;   ///< Nonzero pixels in mask.
  std::vector<double> soma_x;       ///< Final soma x positions, one per cell.
  std::vector<double> soma_y;       ///< Final soma y positions, one per cell.
  std::vector<int> cell_morphology; ///< Per-cell index into morphology_paths.
  std::vector<double>
      template_max_path; ///< Per template; 0 for ones never sampled.
  std::vector<std::uint8_t>
      mask; ///< grid_side^2 mask, row-major, row 0 at y = -L/2.

  /// The speed-coupling calibration actually used once the NaN defaults
  /// resolve; 0/0 when off. Logged so a run records its transfer function.
  double speed_centre_used = 0.0;
  double speed_slope_used = 0.0;

  /// How much of the run a stiffness map was in effect for. Carry-forward means
  /// no map applies before the earliest one, so a DIV-5-only map on a run to
  /// DIV 5 couples one step in 120 and returns a null that reads like
  /// "stiffness does nothing". Counted so the CLI can warn. Both 0 when not
  /// migrating or not coupled.
  std::size_t migration_steps = 0; ///< Timesteps taken.
  std::size_t steps_with_stiff_map =
      0; ///< Of those, ones with a map in effect.

  /// Where the somata ended up relative to the scaffold: the closest thing to
  /// what Li 2025 measures, and the readout that shows whether a steering rule
  /// worked. Occupancy cannot; it mixes "reached the scaffold" with "spread
  /// along it". Filled in whenever a geometry is loaded, coupled or not, since
  /// reading it changes no simulated number. -1 with no geometry.
  int soma_scaffold_count = -1; ///< Somata inside scaffold compartments.
  double soma_scaffold_mean_dist_um =
      0.0; ///< Mean distance to nearest scaffold, um; 0 for cells on it.
};

/**
 * @brief Check a Params before running.
 *
 * Also rejects couplings that could only ever be inert, a weight with no maps,
 * maps starting after the run ends, a steering rule with migration off, which
 * used to report a null that looked like biology, and a --res so coarse the
 * raster rounds to zero pixels, which used to print "-nan %" and exit 0.
 *
 * @param p The parameters to check.
 * @return Empty string when usable, otherwise a message naming the bad field.
 */
std::string validate(const Params &p);

/**
 * @brief Run the simulation.
 * @param p The parameters for this run.
 * @return Occupancy percent, final positions, and the coverage mask.
 * @throws std::invalid_argument If validate(p) is non-empty.
 * @throws std::runtime_error If a morphology cannot be opened or parsed.
 */
Result run(const Params &p);

} // namespace occupancy

#endif // OCCUPANCY_SIMULATION_HPP
