#ifndef OCCUPANCY_DISTANCE_TRANSFORM_HPP
#define OCCUPANCY_DISTANCE_TRANSFORM_HPP
// Exact Euclidean distance transform. Used to dilate centrelines into areas.
#include <cstdint>
#include <vector>

namespace occupancy
{

/**
 * @brief 1D squared distance transform: d[q] = min_p(f[p] + (p-q)^2).
 *
 * Felzenszwalb & Huttenlocher 2012, Theory of Computing 8:415. Exact, O(n);
 * single-letter names are the paper's.
 *
 * @param f Input values, length n.
 * @param d Output values, length n.
 * @param n Number of samples.
 * @param v Scratch, length n. Not read on input.
 * @param z Scratch, length n + 1. Not read on input.
 */
void distance_transform_1d(const double *f, double *d, int n, int *v,
                           double *z);

/**
 * @brief Dilate a binary grid by a disk of radius radius_px.
 *
 * dilate(A u B) == dilate(A) u dilate(B), so every timestep's centrelines are
 * collected and dilated once at the end. That is what keeps cumulative mode
 * cheap.
 *
 * @param feature Binary input, row-major, grid_side^2. 1 = feature pixel.
 * @param mask Output, same size, overwritten.
 * @param grid_side Side of the square grid, pixels.
 * @param radius_px Dilation radius, pixels.
 */
void dilate_by_edt(const std::vector<uint8_t> &feature,
                   std::vector<uint8_t> &mask, int grid_side, double radius_px);

} // namespace occupancy

#endif // OCCUPANCY_DISTANCE_TRANSFORM_HPP
