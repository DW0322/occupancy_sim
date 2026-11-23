#include "check.hpp"
#include "occupancy/clustering.hpp"
#include <cmath>
#include <random>

int main()
{
  std::cout << "test_clustering\n";

  // --- convex hull -------------------------------------------------------
  {
    // A unit square with an interior point: the point must not reach the hull.
    std::vector<double> x{0, 100, 100, 0, 50};
    std::vector<double> y{0, 0, 100, 100, 50};
    occupancy::Region r = occupancy::soma_hull(x, y);
    check::near(r.area_um2, 10000.0, 1e-6, "hull area of a 100 um square");
    check::near(r.perimeter_um, 400.0, 1e-6, "hull perimeter of that square");
  }
  {
    // Collinear points enclose nothing; the caller falls back on the field.
    std::vector<double> x{0, 10, 20, 30}, y{0, 10, 20, 30};
    occupancy::Region r = occupancy::soma_hull(x, y);
    check::near(r.area_um2, 0.0, 1e-12, "collinear points have no hull area");
  }
  {
    std::vector<double> x{0, 1}, y{0, 1};
    check::near(occupancy::soma_hull(x, y).area_um2, 0.0, 1e-12,
                "two points have no hull area");
  }
  {
    // Right triangle, area 1/2 * b * h, perimeter 3-4-5.
    std::vector<double> x{0, 30, 0}, y{0, 0, 40};
    occupancy::Region r = occupancy::soma_hull(x, y);
    check::near(r.area_um2, 600.0, 1e-6, "hull area of a 3-4-5 triangle");
    check::near(r.perimeter_um, 120.0, 1e-6, "hull perimeter of it");
  }

  // --- the statistic -----------------------------------------------------
  {
    // Four cells on a 100 um lattice, measured over their own square. Spacing
    // is exactly 100 um, so the observed mean NN is 100 by construction.
    std::vector<double> x{0, 100, 0, 100}, y{0, 0, 100, 100};
    occupancy::ClusteringStats s =
        occupancy::clustering_stats(x, y, 10000.0, 400.0);
    check::ok(s.n == 4, "n counts the cells");
    check::near(s.mean_nn_um, 100.0, 1e-9, "mean NN of a 100 um lattice");
    check::ok(s.clark_evans > 1.0, "a lattice reads as spread, not clumped");
  }
  {
    // n < 2 has no neighbour to measure.
    std::vector<double> x{0.0}, y{0.0};
    check::ok(occupancy::clustering_stats(x, y, 1.0, 4.0).n == 0,
              "a lone cell yields no statistic");
  }
  {
    // THE REGRESSION THIS FILE EXISTS FOR: the same cells, scored over their
    // own hull and over a field 100x larger. The field version reports heavy
    // clumping purely because the denominator is wrong.
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-150.0, 150.0);
    std::vector<double> x(25), y(25);
    for (int i = 0; i < 25; ++i)
    {
      x[(std::size_t)i] = u(rng);
      y[(std::size_t)i] = u(rng);
    }
    occupancy::Region h = occupancy::soma_hull(x, y);
    const double ce_hull =
        occupancy::clustering_stats(x, y, h.area_um2, h.perimeter_um)
            .clark_evans;
    const double ce_field =
        occupancy::clustering_stats(x, y, 3000.0 * 3000.0, 4.0 * 3000.0)
            .clark_evans;
    // The hull is defined BY the points, so it undercovers the region they
    // were drawn from and the null sits above 1: measured over 2000 CSR
    // trials it is 1.23 at N=25 (3.39 at N=3, 1.06 at N=120). Read a hull
    // Clark-Evans against that, not against 1, the band here is the N=25 null
    // plus room for one draw's scatter.
    check::ok(ce_hull > 0.9 && ce_hull < 1.6,
              "random scatter over its own hull lands in the N=25 null band");
    check::ok(ce_field < 0.2,
              "the same scatter over the whole field reads as clumped");
    check::ok(ce_hull > 4.0 * ce_field,
              "the denominator, not the arrangement, drives the difference");
  }

  return check::summary("test_clustering");
}
