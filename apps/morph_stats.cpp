// Morphometry over the engine's own loader. One CSV row per file on stdout:
//   file,segments,max_path_um,total_xy_um,bbox_w_um,bbox_h_um,spacing_um
//
// This is where the territory radius R comes from: half the mean neurite
// spacing, where spacing = bbox area / total projected length. It goes through
// occupancy::load_template, so the projection is exactly what the simulation
// stamps and the estimate can't drift from the engine.

#include "occupancy/segment.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    std::cerr << "usage: " << argv[0] << " <morphology.asc> [more...]\n";
    return 2;
  }
  std::cout << "file,segments,max_path_um,total_xy_um,bbox_w_um,bbox_h_um,"
               "spacing_um\n";
  int failures = 0;
  for (int i = 1; i < argc; ++i)
  {
    try
    {
      occupancy::MorphologyTemplate t = occupancy::load_template(argv[i]);
      double len = 0.0;
      double xmin = std::numeric_limits<double>::infinity(), xmax = -xmin;
      double ymin = xmin, ymax = -xmin;
      for (const auto &s : t.segments)
      {
        len += std::hypot(s.x1 - s.x0, s.y1 - s.y0);
        xmin = std::min({xmin, s.x0, s.x1});
        xmax = std::max({xmax, s.x0, s.x1});
        ymin = std::min({ymin, s.y0, s.y1});
        ymax = std::max({ymax, s.y0, s.y1});
      }
      const double w = xmax - xmin, h = ymax - ymin;
      const double spacing = (len > 0.0) ? (w * h) / len : 0.0;
      std::cout << argv[i] << ',' << t.segments.size() << ','
                << t.max_path_distance << ',' << len << ',' << w << ',' << h
                << ',' << spacing << '\n';
    }
    catch (const std::exception &e)
    {
      std::cerr << "FAILED " << argv[i] << ": " << e.what() << '\n';
      ++failures;
    }
  }
  return failures == 0 ? 0 : 1;
}
