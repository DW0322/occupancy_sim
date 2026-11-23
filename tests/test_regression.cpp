#include "check.hpp"
#include "occupancy/simulation.hpp"
#include <cstdlib>

int main(int argc, char **argv)
{
  std::cout << "test_regression\n";
  std::string morph = (argc > 1) ? argv[1] : "data/morph.asc";

  occupancy::Params p;
  p.morphology_paths = {morph};
  p.n_cells = 15;
  p.scaffold_um = 1500.0;
  p.territory_um = 17.2;
  p.resolution_um = 0.5;
  p.div_days = 15.0;
  p.seed = 1;
  p.seed_fraction = 1.0;
  p.mig_speed_um_per_h = 0.0;

  occupancy::Result r = occupancy::run(p);
  check::near(r.occupancy_percent, 81.7524, 1e-3,
              "reveal-only N=15 DIV15 seed=1");

  p.mig_speed_um_per_h = 10.0;
  p.mig_persist_um = 50.0;
  p.dt_hours = 1.0;
  p.occupancy_mode = occupancy::OccupancyMode::Cumulative;
  check::near(occupancy::run(p).occupancy_percent, 99.9998, 1e-3,
              "migration cum N=15 DIV15 seed=1");

  // The scaffold-seek decision rule, on an in-code left-half geometry so the
  // pin doesn't depend on a data file. Weight 0 must be the migration pin
  // exactly; that's the inertness contract every coupling honours.
  {
    occupancy::Params s = p;
    const int geometry_side = 10;
    s.scaffold_nx = s.scaffold_ny = geometry_side;
    s.scaffold_cells.assign((std::size_t)geometry_side * geometry_side, 0);
    for (int cy = 0; cy < geometry_side; ++cy)
      for (int cx = 0; cx < geometry_side / 2; ++cx)
        s.scaffold_cells[(std::size_t)cy * geometry_side + cx] = 1;
    s.seek_mode = occupancy::SeekMode::Scaffold;
    s.seek_weight = 0.0;
    s.p_seek = 1.0;
    check::near(occupancy::run(s).occupancy_percent, 99.9998, 1e-3,
                "seek=scaffold weight 0 == migration pin (inert)");
    s.seek_weight = 0.3;
    check::near(occupancy::run(s).occupancy_percent, 94.0087, 1e-3,
                "seek=scaffold weight 0.3, left-half geometry");
  }

  return check::summary("test_regression");
}
