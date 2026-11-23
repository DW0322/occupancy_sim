#ifndef OCCUPANCY_SEGMENT_HPP
#define OCCUPANCY_SEGMENT_HPP
// Loading morphologies and projecting them to 2D.
#include <string>
#include <vector>

namespace occupancy
{

/// One straight piece of neurite, projected to the XY plane.
struct Segment
{
  double x0, y0;
  double x1, y1;
  double d0, d1;
};

/**
 * @brief Load a morphology and project it onto the 2D print surface.
 *
 * Drops Z, recentres on the soma, records each segment's path distance from it,
 * and sorts by d0 so a caller scanning by reveal fraction can stop at the
 * front. Morphologies are Blue Brain NMC Portal reconstructions (Markram 2015,
 * Ramaswamy 2015), citing them is a condition of use.
 *
 * @param path Morphology file (.asc/.swc/.h5).
 * @param max_path_distance Output: longest soma-to-tip path, um.
 * @return The segments, sorted ascending by d0.
 */
std::vector<Segment> load_segments(const std::string &path,
                                   double &max_path_distance);

/**
 * @brief One loaded morphology, ready to stamp.
 *
 * Cells sample from a pool so a population is not N rotations of one template.
 * Mean occupancy spans 13.9 points across templates against about 1 point of
 * standard error: which template you load matters more than which seed.
 */
struct MorphologyTemplate
{
  std::string path;               ///< Source file path.
  std::vector<Segment> segments;  ///< Segments, sorted by d0.
  double max_path_distance = 0.0; ///< Longest soma-to-tip path, um.
};

/**
 * @brief Load a morphology into a MorphologyTemplate.
 * @param path Morphology file (.asc/.swc/.h5).
 * @return The loaded template.
 */
MorphologyTemplate load_template(const std::string &path);

/// Bins in an arbor reach table; index = round(reveal_fraction *
/// (kReachBins-1)).
constexpr int kReachBins = 101;

/**
 * @brief How far the revealed arbor actually reaches, per reveal fraction.
 *
 * @param segs Segments, soma-centred, sorted by d0 (as load_segments returns).
 * @param max_path_distance Longest soma-to-tip path, um.
 * @return kReachBins radii in um, ascending in reveal fraction. Empty input
 *         gives all zeros, which reads as "reaches nobody".
 */
std::vector<double> arbor_reach_table(const std::vector<Segment> &segs,
                                      double max_path_distance);

} // namespace occupancy

#endif // OCCUPANCY_SEGMENT_HPP
