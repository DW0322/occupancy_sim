#include "occupancy/clustering.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace occupancy
{

ClusteringStats clustering_stats(const std::vector<double> &x,
                                 const std::vector<double> &y, double area_um2,
                                 double perimeter_um)
{
  ClusteringStats s;
  const std::size_t n = x.size() < y.size() ? x.size() : y.size();
  if (n < 2 || area_um2 <= 0.0)
    return s;
  s.n = (int)n;

  double sum_nn = 0.0;
  for (std::size_t i = 0; i < n; ++i)
  {
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < n; ++j)
    {
      if (j == i)
        continue;
      const double dx = x[j] - x[i], dy = y[j] - y[i];
      const double d2 = dx * dx + dy * dy;
      if (d2 < best)
        best = d2;
    }
    sum_nn += std::sqrt(best);
  }
  s.mean_nn_um = sum_nn / (double)n;

  // Donnelly 1978's correction to the Poisson expectation 0.5*sqrt(A/n). The
  // two boundary terms matter most at small n, which is exactly where this
  // model runs.
  const double dn = (double)n;
  s.expected_nn_um = 0.5 * std::sqrt(area_um2 / dn) +
                     0.0514 * perimeter_um / dn +
                     0.041 * perimeter_um / (dn * std::sqrt(dn));
  s.clark_evans =
      s.expected_nn_um > 0.0 ? s.mean_nn_um / s.expected_nn_um : 0.0;
  return s;
}

Region soma_hull(const std::vector<double> &x, const std::vector<double> &y)
{
  Region r;
  const std::size_t n = x.size() < y.size() ? x.size() : y.size();
  if (n < 3)
    return r;

  // Andrew's monotone chain: sort by x then y, sweep the lower hull and the
  // upper hull, concatenate. O(n log n), and n is the cell count.
  std::vector<std::size_t> idx(n);
  for (std::size_t i = 0; i < n; ++i)
    idx[i] = i;
  std::sort(idx.begin(), idx.end(),
            [&](std::size_t a, std::size_t b)
            { return x[a] != x[b] ? x[a] < x[b] : y[a] < y[b]; });

  // >0 when the turn a->b->c is counter-clockwise.
  auto cross = [&](std::size_t a, std::size_t b, std::size_t c)
  {
    return (x[b] - x[a]) * (y[c] - y[a]) - (y[b] - y[a]) * (x[c] - x[a]);
  };

  std::vector<std::size_t> hull(2 * n);
  std::size_t k = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    while (k >= 2 && cross(hull[k - 2], hull[k - 1], idx[i]) <= 0.0)
      --k;
    hull[k++] = idx[i];
  }
  const std::size_t lower = k + 1;
  for (std::size_t i = n - 1; i-- > 0;)
  {
    while (k >= lower && cross(hull[k - 2], hull[k - 1], idx[i]) <= 0.0)
      --k;
    hull[k++] = idx[i];
  }
  // Last point repeats the first, so the vertex count is k-1.
  if (k < 4)
    return r; // every point collinear: no area to speak of
  hull.resize(k);

  // Shoelace for the area, edge lengths for the perimeter, over the closed ring.
  double twice_area = 0.0, perim = 0.0;
  for (std::size_t i = 0; i + 1 < k; ++i)
  {
    const std::size_t a = hull[i], b = hull[i + 1];
    twice_area += x[a] * y[b] - x[b] * y[a];
    perim += std::hypot(x[b] - x[a], y[b] - y[a]);
  }
  r.area_um2 = std::fabs(twice_area) * 0.5;
  r.perimeter_um = perim;
  return r;
}

} // namespace occupancy
