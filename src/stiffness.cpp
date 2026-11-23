#include "occupancy/stiffness.hpp"
#include "csv_fields.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace occupancy
{

using csv::clean_field;
using csv::split;
using csv::to_double;
using csv::to_int;

namespace
{
[[noreturn]] void fail(const std::string &path, int line, const std::string &m)
{
  throw std::runtime_error("stiffness map " + path + ":" +
                           std::to_string(line) + ": " + m);
}

// every stiffness file must open with exactly this line
const char *const kVersionLine = "# occupancy-sim stiffness map v1";
} // namespace

int StiffnessMap::select(double run_div) const
{
  int best = -1;
  for (std::size_t i = 0; i < divs.size(); ++i)
    if (divs[i] <= run_div)
      best = (int)i;
  return best;
}

StiffnessMap load_stiffness_map(const std::string &path)
{
  std::ifstream f(path);
  if (!f)
    throw std::runtime_error("cannot open stiffness map " + path);

  StiffnessMap m;
  std::string line;
  int ln = 0;

  // first non-blank line has to be the v1 marker
  bool version_seen = false;
  bool header_seen = false;
  bool grid_seen = false;
  bool scaffold_seen = false;
  bool units_seen = false;
  bool step_seen = false;

  // div -> (cy*nx+cx) -> value. Completeness is checked at the end so the error
  // can name the missing cell.
  std::map<double, std::map<int, double>> cells;

  while (std::getline(f, line))
  {
    ++ln;
    // strip trailing CR from Windows-written files
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (line.empty())
      continue;

    if (!version_seen)
    {
      if (line != kVersionLine)
        fail(path, ln,
             std::string("first non-blank line must be exactly '") +
                 kVersionLine + "' (got '" + line + "')");
      m.version = "v1";
      version_seen = true;
      continue;
    }
    if (line[0] == '#')
      continue; // skip other comments

    auto t = split(line);
    if (!header_seen)
    {
      // gotta start reading the meta data
      if (t.size() == 3 && t[0] == "grid")
      {
        if (!to_int(t[1], m.nx) || !to_int(t[2], m.ny))
          fail(path, ln, "grid sides must be whole numbers");
        if (m.nx < 1 || m.ny < 1)
          fail(path, ln, "grid sides must be >= 1");
        grid_seen = true;
        continue;
      }
      if (t.size() == 2 && t[0] == "scaffold_um")
      {
        if (!to_double(t[1], m.scaffold_um))
          fail(path, ln, "scaffold_um must be a number");
        if (!(m.scaffold_um > 0))
          fail(path, ln, "scaffold_um must be > 0");
        scaffold_seen = true;
        continue;
      }
      if (t.size() == 2 && t[0] == "units")
      {
        if (!clean_field(t[1]))
          fail(path, ln, "units must not be empty or padded with spaces");
        m.units = t[1];
        units_seen = true;
        continue;
      }
      if (t.size() == 2 && t[0] == "div_step")
      {
        if (!to_int(t[1], m.div_step))
          fail(path, ln, "div_step must be a whole number");
        // any positive cadence: the producer's pipeline sets it, we only
        // demand it be declared so the gap check below has a ruler
        if (m.div_step < 1)
          fail(path, ln,
               "div_step must be >= 1 (days between consecutive maps); got " +
                   t[1]);
        step_seen = true;
        continue;
      }
      if (t.size() == 4 && t[0] == "div" && t[1] == "cx" && t[2] == "cy" &&
          t[3] == "stiffness")
      {
        if (!grid_seen)
          fail(path, ln, "data header before 'grid,NX,NY' metadata line");
        if (!scaffold_seen)
          fail(path, ln, "data header before 'scaffold_um,L' metadata line");
        if (!units_seen)
          fail(path, ln, "data header before 'units,U' metadata line");
        if (!step_seen)
          fail(path, ln, "data header before 'div_step,N' metadata line");
        header_seen = true;
        continue;
      }
      fail(path, ln,
           "unrecognised metadata line '" + line +
               "' (expected grid,NX,NY / scaffold_um,L / units,U / "
               "div_step,N / the div,cx,cy,stiffness header)");
    }

    // data region
    if (t.size() != 4)
      fail(path, ln, "expected 4 fields div,cx,cy,stiffness");
    double dv, val;
    int cx, cy;
    if (!to_double(t[0], dv) || !to_int(t[1], cx) || !to_int(t[2], cy) ||
        !to_double(t[3], val))
      fail(path, ln, "unparseable numeric field in '" + line + "'");
    // isfinite too, not just >= 0: "nan" parses, compares false to everything
    // (so every guard below would silently pass) and is UB as a std::map key
    if (!std::isfinite(dv) || dv < 0)
      fail(path, ln, "div must be finite and >= 0");
    if (cx < 0 || cx >= m.nx || cy < 0 || cy >= m.ny)
      fail(path, ln,
           "compartment (" + std::to_string(cx) + "," + std::to_string(cy) +
               ") outside declared grid " + std::to_string(m.nx) + "x" +
               std::to_string(m.ny));
    if (!std::isfinite(val) || val < 0)
      fail(path, ln, "stiffness must be finite and >= 0");
    auto &grid = cells[dv];
    int idx = cy * m.nx + cx;
    if (!grid.emplace(idx, val).second)
      fail(path, ln,
           "duplicate cell (div=" + t[0] + ", cx=" + t[1] + ", cy=" + t[2] +
               ")");
  }

  // error checks for the map data

  if (!version_seen)
    throw std::runtime_error("stiffness map " + path +
                             ": empty file (no version line)");
  if (!header_seen)
    throw std::runtime_error("stiffness map " + path +
                             ": no div,cx,cy,stiffness data header");
  if (cells.empty())
    throw std::runtime_error("stiffness map " + path + ": no data rows");

  const int expect = m.nx * m.ny;
  for (const auto &kv : cells)
  {
    if ((int)kv.second.size() != expect)
    {
      // name one missing cell so whoever made the file can find the hole
      std::string missing = "?";
      for (int i = 0; i < expect; ++i)
        if (!kv.second.count(i))
        {
          missing = "(" + std::to_string(i % m.nx) + "," +
                    std::to_string(i / m.nx) + ")";
          break;
        }
      throw std::runtime_error(
          "stiffness map " + path + ": div " + std::to_string(kv.first) +
          " has " + std::to_string(kv.second.size()) + " of " +
          std::to_string(expect) + " cells (first missing: " + missing +
          ") - every DIV must carry the complete grid");
    }
    m.divs.push_back(kv.first);
    std::vector<double> flat((std::size_t)expect);
    for (const auto &c : kv.second)
      flat[(std::size_t)c.first] = c.second;
    m.maps.push_back(std::move(flat));
  }
  // std::map iterates ascending, so divs are already sorted
  //
  // Cadence check: every DIV must be a multiple of div_step with no gaps.
  // Without this a dropped block would be invisible, carry-forward would hold
  // the previous map over the hole and the figure would still look fine.
  for (double d : m.divs)
  {
    const double q = d / m.div_step;
    if (std::fabs(q - std::llround(q)) > 1e-9)
      throw std::runtime_error("stiffness map " + path + ": div " +
                               std::to_string(d) +
                               " is not a multiple of the declared div_step " +
                               std::to_string(m.div_step));
  }
  for (std::size_t i = 1; i < m.divs.size(); ++i)
  {
    const double expected_div = m.divs[i - 1] + m.div_step;
    if (std::fabs(m.divs[i] - expected_div) > 1e-9)
      throw std::runtime_error(
          "stiffness map " + path + ": gap in the DIV sequence - div_step " +
          std::to_string(m.div_step) + " after div " +
          std::to_string(m.divs[i - 1]) + " requires " +
          std::to_string(expected_div) + ", found " +
          std::to_string(m.divs[i]));
  }
  return m;
}

std::string validate_stiffness_against(const StiffnessMap &m,
                                       int compartments_per_side,
                                       double scaffold_um)
{
  std::ostringstream e;
  if (m.nx != compartments_per_side || m.ny != compartments_per_side)
    e << "map grid " << m.nx << "x" << m.ny
      << " does not match the run's compartment grid " << compartments_per_side
      << "x" << compartments_per_side
      << " (set --compartments to the map, or regenerate the map); ";
  const double rel = std::fabs(m.scaffold_um - scaffold_um) /
                     std::max(1.0, std::fabs(scaffold_um));
  if (rel > 1e-6)
    e << "map scaffold_um " << m.scaffold_um << " does not match the run's "
      << scaffold_um << "; ";
  std::string s = e.str();
  if (s.size() >= 2)
    s.erase(s.size() - 2);
  return s;
}

} // namespace occupancy
