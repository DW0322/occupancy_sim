#include "check.hpp"
#include "occupancy/simulation.hpp"
#include <stdexcept>

int main()
{
  std::cout << "test_params\n";

  occupancy::Params p;
  p.morphology_paths = {"data/morph.asc"};
  check::ok(occupancy::validate(p).empty(), "defaults are valid");

  // the exact mistake above: a seed radius in um passed as a fraction
  occupancy::Params bad = p;
  bad.seed_fraction = 150.0;
  check::ok(!occupancy::validate(bad).empty(),
            "seed_fraction=150 (a radius in um) is rejected");

  occupancy::Params zero = p;
  zero.seed_fraction = 0.0;
  check::ok(!occupancy::validate(zero).empty(), "seed_fraction=0 is rejected");

  occupancy::Params nores = p;
  nores.resolution_um = 0.0;
  check::ok(!occupancy::validate(nores).empty(), "resolution=0 is rejected");

  occupancy::Params nocells = p;
  nocells.n_cells = 0;
  check::ok(!occupancy::validate(nocells).empty(), "n_cells=0 is rejected");

  occupancy::Params dt0 = p;
  dt0.dt_hours = 0.0;
  check::ok(!occupancy::validate(dt0).empty(), "dt=0 is rejected");

  occupancy::Params sw = p;
  sw.seek_weight = 1.5;
  check::ok(!occupancy::validate(sw).empty(), "seek_weight>1 is rejected");

  occupancy::Params ps = p;
  ps.p_seek = -0.1;
  check::ok(!occupancy::validate(ps).empty(), "p_seek<0 is rejected");

  // seek=scaffold is only meaningful with a geometry to steer toward
  occupancy::Params seek_no_geom = p;
  seek_no_geom.seek_mode = occupancy::SeekMode::Scaffold;
  seek_no_geom.mig_speed_um_per_h = 10.0;
  check::ok(!occupancy::validate(seek_no_geom).empty(),
            "seek=scaffold without a geometry is rejected");

  // couplings that could only ever be silently inert are rejected loudly
  occupancy::Params speed_no_maps = p;
  speed_no_maps.stiff_speed_weight = 0.5;
  speed_no_maps.mig_speed_um_per_h = 10.0;
  check::ok(!occupancy::validate(speed_no_maps).empty(),
            "stiff_speed_weight without maps is rejected");

  occupancy::Params late_maps = p;
  late_maps.mig_speed_um_per_h = 10.0;
  late_maps.div_days = 15.0;
  late_maps.stiff_nx = late_maps.stiff_ny = 2;
  late_maps.stiff_divs = {20.0, 25.0}; // every map after the whole run
  late_maps.stiff_maps = {std::vector<double>(4, 1.0),
                          std::vector<double>(4, 2.0)};
  late_maps.stiff_weight = 0.9;
  check::ok(!occupancy::validate(late_maps).empty(),
            "coupling whose maps all start after the run is rejected");

  occupancy::Params no_mig = late_maps;
  no_mig.stiff_divs = {5.0, 10.0}; // maps fine now, but nothing migrates
  no_mig.mig_speed_um_per_h = 0.0;
  check::ok(!occupancy::validate(no_mig).empty(),
            "stiffness coupling with migration off is rejected");

  // A raster that rounds to zero pixels: occupancy is 0/0 and the engine used
  // to print "-nan %" at exit 0. --res is um PER PIXEL, so the wrong unit is
  // the realistic way in.
  {
    occupancy::Params px = p;
    px.scaffold_um = 1500.0;
    px.resolution_um = 5000.0; // e.g. res given in the wrong unit
    const std::string why = occupancy::validate(px);
    check::ok(!why.empty(), "a zero-pixel raster is rejected, not run to NaN");
    check::ok(why.find("PER PIXEL") != std::string::npos,
              "and the message says which way round --res goes");
    px.resolution_um = 1500.0; // exactly one pixel is degenerate but valid
    check::ok(occupancy::validate(px).empty(),
              "a one-pixel raster is still allowed");
  }

  // parse round trips
  check::ok(occupancy::parse_seek_mode("scaffold") ==
                occupancy::SeekMode::Scaffold,
            "parse scaffold");
  check::ok(occupancy::parse_seek_mode("none") == occupancy::SeekMode::None,
            "parse none");
  check::ok(occupancy::parse_occupancy_mode("inst") ==
                occupancy::OccupancyMode::Instant,
            "parse inst");
  check::ok(occupancy::parse_occupancy_mode("cum") ==
                occupancy::OccupancyMode::Cumulative,
            "parse cum");
  check::ok(occupancy::parse_cohere_gate("scaffold") ==
                occupancy::CohereGate::Scaffold,
            "parse cohere gate scaffold");

  auto rejects = [](auto parse, const char *text)
  {
    try
    {
      parse(std::string(text));
      return false;
    }
    catch (const std::invalid_argument &)
    {
      return true;
    }
  };
  check::ok(rejects(occupancy::parse_seek_mode, "scafold"),
            "a misspelled seek mode is rejected, not defaulted to none");
  check::ok(rejects(occupancy::parse_occupancy_mode, "cumm"),
            "a misspelled occupancy mode is rejected, not defaulted to cum");
  check::ok(rejects(occupancy::parse_cohere_gate, "scafold"),
            "a misspelled cohere gate is rejected, not defaulted to none");

  return check::summary("test_params");
}
