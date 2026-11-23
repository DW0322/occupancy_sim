// The compartment layer has to be inert and aggregate exactly before anything
// from it is trusted. Four checks:
//   1. On the pinned run, summing compartments reproduces covered_pixels and
//   the
//      global occupancy bitwise on the doubles, not just close.
//   2. Same holds when the compartment count doesn't divide the raster side.
//   3. A tiny synthetic mask with values worked out by hand.
//   4. Binning is deterministic: same mask twice -> identical grids.
//
//  !! This test was made using claude, hand calculated is
//  approximate generations made by claude not proven measurements !!

#include "check.hpp"
#include "occupancy/compartments.hpp"
#include "occupancy/simulation.hpp"

#include <cstdint>
#include <vector>

using occupancy::bin_mask;
using occupancy::CompartmentGrid;

namespace
{
void invariants(const CompartmentGrid &g, const occupancy::Result &r,
                const std::string &tag)
{
  check::ok(g.covered_sum() == r.covered_pixels,
            tag + ": sum(covered) == covered_pixels exactly");
  check::ok(g.total_sum() ==
                (std::size_t)r.grid_side * (std::size_t)r.grid_side,
            tag + ": sum(total) == G*G exactly");
  check::ok(g.aggregate_percent() == r.occupancy_percent,
            tag + ": aggregate == global occupancy bit-exact");
}
} // namespace

int main(int argc, char **argv)
{
  std::cout << "test_compartments\n";
  std::string morph = (argc > 1) ? argv[1] : "data/morph.asc";

  occupancy::Params p;
  p.morphology_paths = {morph};
  p.n_cells = 15;
  p.scaffold_um = 1500.0;
  p.territory_um = 17.2;
  p.resolution_um = 0.5;
  p.div_days = 15.0;
  p.seed = 1;

  occupancy::Result r = occupancy::run(p);
  check::near(r.occupancy_percent, 81.7524, 1e-3,
              "pinned run reproduced before binning");

  for (int n : {30, 7, 31, 1})
  {
    CompartmentGrid g = bin_mask(r.mask, r.grid_side, n);
    invariants(g, r, "n=" + std::to_string(n));
  }

  CompartmentGrid g1 = bin_mask(r.mask, r.grid_side, 1);
  check::ok(g1.occupancy_percent(0, 0) == r.occupancy_percent,
            "n=1 single compartment == global");

  std::vector<std::uint8_t> m = {1, 1, 0, 0, 1, 0, 0, 0,
                                 0, 0, 0, 1, 0, 0, 1, 1};
  CompartmentGrid s = bin_mask(m, 4, 2);
  check::ok(s.covered[0] == 3 && s.total[0] == 4, "synthetic (0,0) = 3/4");
  check::ok(s.covered[1] == 0 && s.total[1] == 4, "synthetic (1,0) = 0/4");
  check::ok(s.covered[2] == 0 && s.total[2] == 4, "synthetic (0,1) = 0/4");
  check::ok(s.covered[3] == 3 && s.total[3] == 4, "synthetic (1,1) = 3/4");
  check::near(s.occupancy_percent(0, 0), 75.0, 0.0, "synthetic percent (0,0)");
  check::near(s.aggregate_percent(), 100.0 * 6.0 / 16.0, 0.0,
              "synthetic aggregate 6/16");

  CompartmentGrid s3 = bin_mask(m, 4, 3);
  check::ok(s3.total_sum() == 16, "n=3 over G=4: every pixel binned once");
  check::ok(s3.covered_sum() == 6, "n=3 over G=4: covered preserved");
  check::ok(s3.total[0] == 4 && s3.total[1] == 2 && s3.total[4] == 1,
            "n=3 over G=4: uneven compartment sizes as expected");

  CompartmentGrid a = bin_mask(r.mask, r.grid_side, 30);
  CompartmentGrid b = bin_mask(r.mask, r.grid_side, 30);
  check::ok(a.covered == b.covered && a.total == b.total,
            "binning is deterministic");

  return check::summary("test_compartments");
}
