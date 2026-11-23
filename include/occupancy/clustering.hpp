#ifndef OCCUPANCY_CLUSTERING_HPP
#define OCCUPANCY_CLUSTERING_HPP
// How aggregated a population of somata is. Post-processing of final soma
// positions; nothing here feeds back into the simulation.
#include <cstddef>
#include <vector>

namespace occupancy
{

/// @brief How clumped a set of somata is, against what chance would give.
struct ClusteringStats
{
  int n = 0;                   ///< Cells the statistic was computed over.
  double mean_nn_um = 0.0;     ///< Observed mean nearest-neighbour distance.
  double expected_nn_um = 0.0; ///< Expected under complete spatial randomness.
  double clark_evans = 0.0; ///< mean_nn / expected_nn. <1 clumped, >1 spread.
};

/**
 * @brief Clark-Evans nearest-neighbour statistic over soma positions. References
 * in readme
 *
 * @param x Soma x positions, um.
 * @param y Soma y positions, um, same length as x.
 * @param area_um2 Area the cells occupy, um^2.
 * @param perimeter_um Perimeter of that area, um. Drives the edge correction.
 * @return The statistic. n < 2 gives all zeros, a lone cell has no neighbour.
 */
ClusteringStats clustering_stats(const std::vector<double> &x,
                                 const std::vector<double> &y, double area_um2,
                                 double perimeter_um);

/// @brief A study region for the Clark-Evans expectation.
struct Region
{
  double area_um2 = 0.0;     ///< Area of the region.
  double perimeter_um = 0.0; ///< Its perimeter, for the edge correction.
};

/**
 * @brief Convex hull of the soma positions, as an area and a perimeter.
 *
 * Clark-Evans needs the region the cells were sampled from, and for a run that
 * seeds a sub-field cluster (--seed-frac < 1) that is neither the whole field
 * nor the seed box: cells are free to leave the box, so the box understates the
 * region and the field overstates it by the ratio of their areas. Seeding a
 * 300 um cluster in a 3000 um field overstates it 100-fold, which alone drives
 * the statistic to "clumped" whatever the cells did.
 *
 * The hull is the standard answer when the study region is "wherever the cells
 * got to". It is defined BY the points, so the outermost somata sit on the
 * boundary, which is exactly the bias Donnelly's perimeter terms correct.
 *
 * @param x Soma x positions, um.
 * @param y Soma y positions, um, same length as x.
 * @return The hull. Fewer than 3 points, or collinear ones, give a zero area;
 *         the caller is expected to fall back on the field.
 */
Region soma_hull(const std::vector<double> &x, const std::vector<double> &y);

} // namespace occupancy

#endif // OCCUPANCY_CLUSTERING_HPP
