#include "occupancy/compartments.hpp"

#include <fstream>
#include <stdexcept>

namespace occupancy
{

std::size_t CompartmentGrid::covered_sum() const
{
  std::size_t sum = 0;
  for (auto count : covered)
    sum += count;
  return sum;
}

std::size_t CompartmentGrid::total_sum() const
{
  std::size_t sum = 0;
  for (auto count : total)
    sum += count;
  return sum;
}

double CompartmentGrid::occupancy_percent(int cx, int cy) const
{
  const std::size_t idx = (std::size_t)cy * nx + cx;
  if (total[idx] == 0)
    return 0.0;
  return 100.0 * (double)covered[idx] / (double)total[idx];
}

double CompartmentGrid::aggregate_percent() const
{
  // same expression as the engine's global occupancy, over the same integers
  return 100.0 * (double)covered_sum() / (double)total_sum();
}

CompartmentGrid bin_mask(const std::vector<std::uint8_t> &mask, int grid_side,
                         int n)
{
  if (grid_side <= 0)
    throw std::invalid_argument("bin_mask: grid_side must be > 0");
  if (n <= 0)
    throw std::invalid_argument("bin_mask: compartment count must be > 0");
  if (n > grid_side)
    throw std::invalid_argument(
        "bin_mask: more compartments than raster pixels per side");
  if (mask.size() != (std::size_t)grid_side * grid_side)
    throw std::invalid_argument("bin_mask: mask size does not match grid_side");

  CompartmentGrid g;
  g.nx = g.ny = n;
  g.grid_side = grid_side;
  g.covered.assign((std::size_t)n * n, 0);
  g.total.assign((std::size_t)n * n, 0);

  // Integer c = p * n / G is floor(p*n/G): monotone and onto [0, n). Precompute
  // it once per row/column so the inner loop is just an add (~14 ms for 3000^2
  // at n=30).
  std::vector<int> pixel_to_comp((std::size_t)grid_side);
  for (int px = 0; px < grid_side; ++px)
    pixel_to_comp[(std::size_t)px] = (int)((std::int64_t)px * n / grid_side);

  for (int y = 0; y < grid_side; ++y)
  {
    const int cy = pixel_to_comp[(std::size_t)y];
    const std::uint8_t *row = mask.data() + (std::size_t)y * grid_side;
    std::size_t *covered_row = g.covered.data() + (std::size_t)cy * n;
    std::size_t *total_row = g.total.data() + (std::size_t)cy * n;
    for (int x = 0; x < grid_side; ++x)
    {
      const int cx = pixel_to_comp[(std::size_t)x];
      total_row[cx] += 1;
      covered_row[cx] += row[x] ? 1 : 0;
    }
  }
  return g;
}

void write_compartments_csv(const std::string &path, const CompartmentGrid &g,
                            const std::vector<std::string> &header_lines)
{
  std::ofstream f(path);
  if (!f)
    throw std::runtime_error("cannot open " + path + " for writing");
  for (const auto &h : header_lines)
    f << "# " << h << "\n";
  f << "cx,cy,covered_px,total_px,occupancy_percent\n";
  for (int cy = 0; cy < g.ny; ++cy)
    for (int cx = 0; cx < g.nx; ++cx)
    {
      const std::size_t i = (std::size_t)cy * g.nx + cx;
      f << cx << "," << cy << "," << g.covered[i] << "," << g.total[i] << ","
        << g.occupancy_percent(cx, cy) << "\n";
    }
}

} // namespace occupancy
