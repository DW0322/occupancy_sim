#include "occupancy/distance_transform.hpp"

namespace occupancy
{

void distance_transform_1d(const double *f, double *d, int n, int *v, double *z)
{
  const double INF = 1e20;
  int k = 0;
  v[0] = 0;
  z[0] = -INF;
  z[1] = INF;
  for (int q = 1; q < n; ++q)
  {
    double s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) /
               (2.0 * q - 2.0 * v[k]);
    while (s <= z[k])
    {
      --k;
      s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) /
          (2.0 * q - 2.0 * v[k]);
    }
    ++k;
    v[k] = q;
    z[k] = s;
    z[k + 1] = INF;
  }
  k = 0;
  for (int q = 0; q < n; ++q)
  {
    while (z[k + 1] < q)
      ++k;
    d[q] = (double)(q - v[k]) * (q - v[k]) + f[v[k]];
  }
}

void dilate_by_edt(const std::vector<uint8_t> &feature,
                   std::vector<uint8_t> &mask, int grid_side, double radius_px)
{
  const double INF = 1e20;
  const double radius_px_sq = radius_px * radius_px;
  std::vector<double> col_dist_sq((std::size_t)grid_side * grid_side);
  std::vector<double> col_in(grid_side), col_out(grid_side);
  std::vector<int> scratch_v(grid_side);
  std::vector<double> scratch_z(grid_side + 1);

  for (int x = 0; x < grid_side; ++x) // pass 1: columns
  {
    for (int y = 0; y < grid_side; ++y)
      col_in[y] = feature[(std::size_t)y * grid_side + x] ? 0.0 : INF;
    distance_transform_1d(col_in.data(), col_out.data(), grid_side,
                          scratch_v.data(), scratch_z.data());
    for (int y = 0; y < grid_side; ++y)
      col_dist_sq[(std::size_t)y * grid_side + x] = col_out[y];
  }
  std::vector<double> row_in(grid_side), row_out(grid_side);
  for (int y = 0; y < grid_side;
       ++y) // pass 2: rows, fed the partial squared distances
  {
    for (int x = 0; x < grid_side; ++x)
      row_in[x] = col_dist_sq[(std::size_t)y * grid_side + x];
    distance_transform_1d(row_in.data(), row_out.data(), grid_side,
                          scratch_v.data(), scratch_z.data());
    for (int x = 0; x < grid_side; ++x)
      mask[(std::size_t)y * grid_side + x] =
          (row_out[x] <= radius_px_sq) ? 1 : 0;
  }
}

} // namespace occupancy
