#ifndef OCCUPANCY_COMPARTMENTS_HPP
#define OCCUPANCY_COMPARTMENTS_HPP
// Binning the finished coverage mask into a grid of square compartments.
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace occupancy
{

/**
 * @brief Which compartment a position falls in, clamped to the grid.
 *
 * Every cue reads its field at this granularity: scaffold mask, stiffness map,
 * gradient stencil, soma readouts. There used to be seven copies of the
 * mapping. The CLI's and the engine's had to agree for SOMA_STIFFNESS to mean
 * anything, and nothing checked that they did.
 *
 * Do not tidy the arithmetic: it is written exactly as those copies had it, so
 * the result stays bit-identical and the pinned numbers hold. Positions are
 * field-centred; out-of-range clamps to the edge rather than wrapping, which is
 * what a soma on the boundary needs.
 *
 * @param x Position x, um, field-centred.
 * @param y Position y, um, field-centred.
 * @param field_um Field side, um.
 * @param nx Compartments across x; must be >= 1.
 * @param ny Compartments across y; must be >= 1.
 * @param cx Output: column, 0 <= cx < nx.
 * @param cy Output: row, 0 <= cy < ny.
 */
inline void compartment_of(double x, double y, double field_um, int nx, int ny,
                           int &cx, int &cy)
{
  cx = (int)((x + field_um / 2) / field_um * nx);
  cy = (int)((y + field_um / 2) / field_um * ny);
  if (cx < 0)
    cx = 0;
  if (cx >= nx)
    cx = nx - 1;
  if (cy < 0)
    cy = 0;
  if (cy >= ny)
    cy = ny - 1;
}

/**
 * @brief The field split into square compartments, each with its own occupancy.
 *
 * Built after the run from Result.mask, so switching this layer on leaves a run
 * bit-identical. Pixels map by integer division c = p * n / G, so every pixel
 * lands in exactly one compartment even when n does not divide G: sum(covered)
 * == covered_pixels and sum(total) == G*G exactly, at any grid size. A mismatch
 * is a bug, not rounding.
 */
struct CompartmentGrid
{
  int nx = 0, ny = 0; ///< Compartments per side (row-major, cy*nx + cx).
  int grid_side = 0;  ///< Side of the source raster, in pixels.
  double resolution_um = 0.0;       ///< Source raster resolution, um per pixel.
  double scaffold_um = 0.0;         ///< Scaffold side, um.
  std::vector<std::size_t> covered; ///< Covered pixels per compartment.
  std::vector<std::size_t> total;   ///< Total pixels per compartment.

  /// @brief Covered pixels summed over every compartment.
  std::size_t covered_sum() const;

  /// @brief Total pixels summed over every compartment; always grid_side^2.
  std::size_t total_sum() const;

  /**
   * @brief Occupancy of one compartment.
   * @param cx Column, 0 <= cx < nx.
   * @param cy Row, 0 <= cy < ny.
   * @return 100 * covered/total for that compartment.
   */
  double occupancy_percent(int cx, int cy) const;

  /// @brief 100 * covered_sum()/total_sum(); must equal the run's occupancy.
  double aggregate_percent() const;
};

/**
 * @brief Bin a coverage mask into a compartment grid.
 * @param mask A grid_side^2 mask in Result.mask layout: row-major, row 0 at
 *        y = -L/2.
 * @param grid_side Side of mask, pixels.
 * @param n Compartments per side.
 * @return The binned grid.
 * @throws std::invalid_argument If grid_side or n is not positive, n exceeds
 *         grid_side, or mask's size does not match grid_side.
 */
CompartmentGrid bin_mask(const std::vector<std::uint8_t> &mask, int grid_side,
                         int n);

/**
 * @brief Write one run's per-compartment occupancy to CSV.
 *
 * Columns cx,cy,covered_px,total_px,occupancy_percent. cy = 0 is mask row 0,
 * the y = -L/2 edge.
 *
 * @param path Output path.
 * @param g The grid to write.
 * @param header_lines Written first, each prefixed "# "; every CSV records the
 *        invocation that made it.
 * @throws std::runtime_error If path cannot be opened.
 */
void write_compartments_csv(const std::string &path, const CompartmentGrid &g,
                            const std::vector<std::string> &header_lines);

} // namespace occupancy

#endif // OCCUPANCY_COMPARTMENTS_HPP
