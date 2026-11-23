#ifndef OCCUPANCY_SCAFFOLD_GEOMETRY_HPP
#define OCCUPANCY_SCAFFOLD_GEOMETRY_HPP
// Scaffold geometry masks: which compartments are scaffold material.
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace occupancy
{

/**
 * Which compartments of the field are scaffold. Format in
 * docs/SCAFFOLD_GEOMETRY_FORMAT.md; the engine validates, never invents.
 */
struct ScaffoldGeometry
{
  int nx = 0, ny = 0;
  double scaffold_um = 0;
  std::string version;
  std::vector<std::uint8_t> is_scaffold;

  std::size_t scaffold_compartments() const;
  double scaffold_fraction() const;
};

ScaffoldGeometry load_scaffold_geometry(const std::string &path);

std::string validate_geometry_against(const ScaffoldGeometry &g,
                                      int compartments_per_side,
                                      double scaffold_um);

} // namespace occupancy

#endif // OCCUPANCY_SCAFFOLD_GEOMETRY_HPP
