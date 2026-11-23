#ifndef OCCUPANCY_SRC_CSV_FIELDS_HPP
#define OCCUPANCY_SRC_CSV_FIELDS_HPP
// Field parsing shared by the stiffness and scaffold-geometry CSV loaders.
// Internal to src/, not part of the public API.
//
// Fields are taken verbatim: no padding, no trailing junk. std::stod and
// std::stoi are laxer than that on their own ("0.5abc" parses as 0.5), so these
// wrappers tighten them up.

#include <cctype>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

namespace occupancy
{
namespace csv
{

/// Split one line on commas, keeping every field exactly as written.
inline std::vector<std::string> split(const std::string &s)
{
  std::vector<std::string> out;
  std::string cur;
  std::istringstream is(s);
  while (std::getline(is, cur, ','))
    out.push_back(cur);
  return out;
}

/// Non-empty and unpadded, i.e. safe to hand to a conversion.
inline bool clean_field(const std::string &s)
{
  return !s.empty() && !std::isspace((unsigned char)s.front()) &&
         !std::isspace((unsigned char)s.back());
}

/// Parse a whole field as a double. False if padded or not fully consumed.
inline bool to_double(const std::string &s, double &out)
{
  if (!clean_field(s))
    return false;
  try
  {
    std::size_t used = 0;
    out = std::stod(s, &used);
    return used == s.size();
  }
  catch (const std::exception &)
  {
    return false;
  }
}

/// Parse a whole field as an int. False if it is padded or not fully consumed.
inline bool to_int(const std::string &s, int &out)
{
  if (!clean_field(s))
    return false;
  try
  {
    std::size_t used = 0;
    out = std::stoi(s, &used);
    return used == s.size();
  }
  catch (const std::exception &)
  {
    return false;
  }
}

} // namespace csv
} // namespace occupancy

#endif // OCCUPANCY_SRC_CSV_FIELDS_HPP
