#include "occupancy/scaffold_geometry.hpp"
#include "csv_fields.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace occupancy
{

using csv::split;
using csv::to_double;
using csv::to_int;

namespace
{
[[noreturn]] void fail(const std::string &path, int line, const std::string &m)
{
  throw std::runtime_error("scaffold geometry " + path + ":" +
                           std::to_string(line) + ": " + m);
}
} // namespace

std::size_t ScaffoldGeometry::scaffold_compartments() const
{
  std::size_t n = 0;
  for (auto v : is_scaffold)
    n += v ? 1 : 0;
  return n;
}

double ScaffoldGeometry::scaffold_fraction() const
{
  if (is_scaffold.empty())
    return 0.0;
  return (double)scaffold_compartments() / (double)is_scaffold.size();
}

ScaffoldGeometry load_scaffold_geometry(const std::string &path)
{
  std::ifstream f(path);
  if (!f)
    throw std::runtime_error("cannot open scaffold geometry " + path);

  ScaffoldGeometry g;
  std::string line;
  int ln = 0;
  bool version_seen = false, header_seen = false;
  bool grid_seen = false, scaffold_seen = false;
  std::map<int, int> cells; // index -> 0/1

  while (std::getline(f, line))
  {
    ++ln;
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (line.empty())
      continue;

    if (!version_seen)
    {
      if (line != "# occupancy-sim scaffold geometry v1")
        fail(path, ln,
             "first non-blank line must be exactly "
             "'# occupancy-sim scaffold geometry v1' (got '" +
                 line + "')");
      g.version = "v1";
      version_seen = true;
      continue;
    }
    if (line[0] == '#')
      continue;

    auto t = split(line);
    if (!header_seen)
    {
      if (t.size() == 3 && t[0] == "grid")
      {
        if (!to_int(t[1], g.nx) || !to_int(t[2], g.ny))
          fail(path, ln, "grid sides must be whole numbers");
        if (g.nx < 1 || g.ny < 1)
          fail(path, ln, "grid sides must be >= 1");
        grid_seen = true;
        continue;
      }
      if (t.size() == 2 && t[0] == "scaffold_um")
      {
        if (!to_double(t[1], g.scaffold_um))
          fail(path, ln, "scaffold_um must be a number");
        if (!(g.scaffold_um > 0))
          fail(path, ln, "scaffold_um must be > 0");
        scaffold_seen = true;
        continue;
      }
      if (t.size() == 3 && t[0] == "cx" && t[1] == "cy" && t[2] == "scaffold")
      {
        if (!grid_seen)
          fail(path, ln, "data header before 'grid,NX,NY' metadata line");
        if (!scaffold_seen)
          fail(path, ln, "data header before 'scaffold_um,L' metadata line");
        header_seen = true;
        continue;
      }
      fail(path, ln,
           "unrecognised metadata line '" + line +
               "' (expected grid,NX,NY / scaffold_um,L / "
               "the cx,cy,scaffold header)");
    }

    if (t.size() != 3)
      fail(path, ln, "expected 3 fields cx,cy,scaffold");
    int cx, cy, v;
    if (!to_int(t[0], cx) || !to_int(t[1], cy) || !to_int(t[2], v))
      fail(path, ln, "unparseable integer field in '" + line + "'");
    if (cx < 0 || cx >= g.nx || cy < 0 || cy >= g.ny)
      fail(path, ln,
           "compartment (" + std::to_string(cx) + "," + std::to_string(cy) +
               ") outside declared grid " + std::to_string(g.nx) + "x" +
               std::to_string(g.ny));
    if (v != 0 && v != 1)
      fail(path, ln,
           "scaffold flag must be exactly 0 or 1 (got '" + t[2] +
               "'). A partial-coverage value is not representable in v1 - "
               "use a finer grid instead.");
    if (!cells.emplace(cy * g.nx + cx, v).second)
      fail(path, ln, "duplicate cell (cx=" + t[0] + ", cy=" + t[1] + ")");
  }

  if (!version_seen)
    throw std::runtime_error("scaffold geometry " + path +
                             ": empty file (no version line)");
  if (!header_seen)
    throw std::runtime_error("scaffold geometry " + path +
                             ": no cx,cy,scaffold data header");

  const int expect = g.nx * g.ny;
  if ((int)cells.size() != expect)
  {
    std::string missing = "?";
    for (int i = 0; i < expect; ++i)
      if (!cells.count(i))
      {
        missing = "(" + std::to_string(i % g.nx) + "," +
                  std::to_string(i / g.nx) + ")";
        break;
      }
    throw std::runtime_error("scaffold geometry " + path + ": " +
                             std::to_string(cells.size()) + " of " +
                             std::to_string(expect) +
                             " cells present (first missing: " + missing +
                             ") - the grid must be complete");
  }

  g.is_scaffold.assign((std::size_t)expect, 0);
  for (const auto &c : cells)
    g.is_scaffold[(std::size_t)c.first] = (std::uint8_t)c.second;

  if (g.scaffold_compartments() == 0)
    throw std::runtime_error("scaffold geometry " + path +
                             ": no compartment is marked as scaffold - the "
                             "occupancy denominator would be zero");
  return g;
}

std::string validate_geometry_against(const ScaffoldGeometry &g,
                                      int compartments_per_side,
                                      double scaffold_um)
{
  std::ostringstream e;
  if (g.nx != compartments_per_side || g.ny != compartments_per_side)
    e << "geometry grid " << g.nx << "x" << g.ny
      << " does not match the run's compartment grid " << compartments_per_side
      << "x" << compartments_per_side
      << " (set --compartments to the file, or regenerate the file); ";
  const double rel = std::fabs(g.scaffold_um - scaffold_um) /
                     std::max(1.0, std::fabs(scaffold_um));
  if (rel > 1e-6)
    e << "geometry scaffold_um " << g.scaffold_um
      << " does not match the run's " << scaffold_um << "; ";
  std::string s = e.str();
  if (s.size() >= 2)
    s.erase(s.size() - 2);
  return s;
}

} // namespace occupancy
