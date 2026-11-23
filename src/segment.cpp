#include "occupancy/segment.hpp"
#include <algorithm>
#include <cmath>
#include <morphio/morphology.h>
#include <morphio/section.h>
#include <morphio/soma.h>
#include <queue>
#include <utility>

namespace occupancy
{

std::vector<Segment> load_segments(const std::string &path,
                                   double &max_path_distance)
{
  morphio::Morphology m(path);
  std::vector<Segment> segs;
  max_path_distance = 0.0;

  // soma centroid, so the arbor comes out relative to the soma
  double sum_x = 0, sum_y = 0;
  std::size_t soma_point_count = 0;
  try
  {
    const auto soma = m.soma();
    const auto &soma_points = soma.points();
    for (const auto &pt : soma_points)
    {
      sum_x += pt[0];
      sum_y += pt[1];
      ++soma_point_count;
    }
  }
  catch (...)
  {
  }
  double soma_x = 0.0, soma_y = 0.0;
  if (soma_point_count > 0)
  {
    soma_x = sum_x / (double)soma_point_count;
    soma_y = sum_y / (double)soma_point_count;
  }

  // BFS over the tree, carrying path distance from the root
  for (const auto &root : m.rootSections())
  {
    std::queue<std::pair<morphio::Section, double>> pending;
    pending.push({root, 0.0});
    while (!pending.empty())
    {
      auto [sec, dist_at_start] = pending.front();
      pending.pop();
      const auto &pts = sec.points();
      double dist_from_soma = dist_at_start;
      for (std::size_t i = 1; i < pts.size(); ++i)
      {
        double x0 = pts[i - 1][0] - soma_x, y0 = pts[i - 1][1] - soma_y;
        double x1 = pts[i][0] - soma_x, y1 = pts[i][1] - soma_y;
        double len = std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
        segs.push_back(
            Segment{x0, y0, x1, y1, dist_from_soma, dist_from_soma + len});
        dist_from_soma += len;
        if (dist_from_soma > max_path_distance)
          max_path_distance = dist_from_soma;
      }
      for (const auto &child : sec.children())
        pending.push({child, dist_from_soma});
    }
  }
  std::sort(segs.begin(), segs.end(),
            [](const Segment &a, const Segment &b) { return a.d0 < b.d0; });
  return segs;
}

std::vector<double> arbor_reach_table(const std::vector<Segment> &segs,
                                      double max_path_distance)
{
  std::vector<double> table((std::size_t)kReachBins, 0.0);
  if (segs.empty() || max_path_distance <= 0.0)
    return table;

  // One radius per segment, taken at its far end, which is how far this piece
  // of neurite reaches into the plane. Paired with d0 so the input's prefix
  // ordering carries over.
  std::vector<double> radius;
  radius.reserve(segs.size());
  for (const Segment &seg : segs)
    radius.push_back(std::sqrt(seg.x1 * seg.x1 + seg.y1 * seg.y1));

  // segs is sorted by d0, so everything revealed at a given fraction is a
  // prefix of it. Further Info on bluebrain documentation
  std::vector<double> revealed;
  revealed.reserve(segs.size());
  std::size_t next = 0;
  for (int bin = 0; bin < kReachBins; ++bin)
  {
    const double frac = (double)bin / (double)(kReachBins - 1);
    const double reveal_um = frac * max_path_distance;
    while (next < segs.size() && segs[next].d0 <= reveal_um)
    {
      revealed.push_back(radius[next]);
      ++next;
    }
    if (revealed.empty())
      continue;
    const std::size_t k = (std::size_t)(0.95 * (double)revealed.size());
    const std::size_t idx = k >= revealed.size() ? revealed.size() - 1 : k;
    std::nth_element(revealed.begin(), revealed.begin() + (std::ptrdiff_t)idx,
                     revealed.end());
    // Running max, not the raw p95. A percentile over a growing set can dip -
    // reveal advances along PATH distance, so a late tortuous dendrite adds
    // points close to the soma and drags the 95th down (0.013 um on the
    // regression morphology). Physically the reach cannot shrink: neurites get
    // revealed, never retracted. Clamping here keeps that exact rather than
    // approximate, so the cue's radius is guaranteed non-decreasing in DIV.
    const double p95 = revealed[idx];
    table[(std::size_t)bin] = bin > 0 && table[(std::size_t)bin - 1] > p95
                                  ? table[(std::size_t)bin - 1]
                                  : p95;
  }
  return table;
}

MorphologyTemplate load_template(const std::string &path)
{
  MorphologyTemplate tmpl;
  tmpl.path = path;
  tmpl.segments = load_segments(path, tmpl.max_path_distance);
  return tmpl;
}

} // namespace occupancy
