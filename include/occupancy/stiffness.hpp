#ifndef OCCUPANCY_STIFFNESS_HPP
#define OCCUPANCY_STIFFNESS_HPP
// Per-compartment stiffness maps over DIV: loading, validation, selection.
#include <string>
#include <vector>

namespace occupancy
{

/**
 * @brief Per-compartment stiffness over DIV. Files come from elsewhere; format
 * in docs/STIFFNESS_FORMAT.md.
 */
struct StiffnessMap
{
  int nx = 0, ny = 0;     ///< Compartment grid declared in the file.
  double scaffold_um = 0; ///< Field side declared in the file, um.
  std::string units;      ///< Declared units, carried verbatim (expect kPa).
  std::string version;    ///< Format version, e.g. "v1".

  int div_step = 0;
  std::vector<double> divs; ///< Ascending, at least one.
  std::vector<std::vector<double>> maps;

  int select(double run_div) const;
};

/**
 * @brief Parse and validate a stiffness file.
 *
 * Checks the version line, the metadata, a complete grid per DIV, each
 * (div,cx,cy) exactly once and finite non-negative values.
 *
 * @param path Path to the CSV (docs/STIFFNESS_FORMAT.md).
 * @return The parsed map.
 * @throws std::runtime_error If the file is missing, malformed, or breaks a
 *         rule. The message names the line and the rule.
 */
StiffnessMap load_stiffness_map(const std::string &path);

/**
 * @brief Check a loaded map against this run's geometry.
 * @param m The loaded map.
 * @param compartments_per_side The run's compartment grid side.
 * @param scaffold_um The run's scaffold side, um.
 * @return Empty string when usable, otherwise the mismatch.
 */
std::string validate_stiffness_against(const StiffnessMap &m,
                                       int compartments_per_side,
                                       double scaffold_um);

} // namespace occupancy

#endif // OCCUPANCY_STIFFNESS_HPP
