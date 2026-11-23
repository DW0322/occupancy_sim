// The stiffness loader is the contract with whoever produces the maps, so the
// rules in docs/STIFFNESS_FORMAT.md are pinned here rather than trusted to the
// parser.

// !! Test generated using claude, take structure at face value not a validated
// or proven geometry only for testing!!

#include "check.hpp"
#include "occupancy/stiffness.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

namespace
{
const char *const kGood = "# occupancy-sim stiffness map v1\n"
                          "# provenance: written by test_stiffness, synthetic\n"
                          "grid,2,2\n"
                          "scaffold_um,1500\n"
                          "units,kPa\n"
                          "div_step,5\n"
                          "div,cx,cy,stiffness\n"
                          "5,0,0,1.00\n"
                          "5,1,0,1.25\n"
                          "5,0,1,1.50\n"
                          "5,1,1,1.75\n"
                          "10,0,0,2.00\n"
                          "10,1,0,2.25\n"
                          "10,0,1,2.50\n"
                          "10,1,1,2.75\n";

std::string write_tmp(const std::string &name, const std::string &body)
{
  std::ofstream(name, std::ios::binary) << body;
  return name;
}

void rejects(const std::string &name, const std::string &body,
             const std::string &what)
{
  const std::string path = write_tmp(name, body);
  bool threw = false;
  try
  {
    occupancy::load_stiffness_map(path);
  }
  catch (const std::runtime_error &)
  {
    threw = true;
  }
  check::ok(threw, what);
  std::remove(path.c_str());
}
} // namespace

int main()
{
  std::cout << "test_stiffness\n";

  // --- a valid file loads and carries its metadata ---------------------------
  const std::string good = write_tmp("good.csv", kGood);
  const occupancy::StiffnessMap m = occupancy::load_stiffness_map(good);

  check::ok(m.version == "v1", "version carried");
  check::ok(m.units == "kPa", "units carried verbatim");
  check::ok(m.div_step == 5, "div_step carried");
  check::ok(m.nx == 2 && m.ny == 2, "grid carried");
  check::near(m.scaffold_um, 1500.0, 1e-9, "scaffold_um carried");
  check::ok(m.divs.size() == 2 && m.divs[0] == 5 && m.divs[1] == 10,
            "two DIVs, ascending");

  check::near(m.maps[0][1], 1.25, 1e-9, "value at (1,0) div5");
  check::near(m.maps[0][2], 1.50, 1e-9, "value at (0,1) div5");
  check::near(m.maps[1][3], 2.75, 1e-9, "value at (1,1) div10");

  // --- select(): carry-forward, never interpolation --------------------------
  check::ok(m.select(4.9) == -1, "select: -1 before earliest map");
  check::ok(m.select(5.0) == 0, "select: exact hit on first map");
  check::ok(m.select(7.0) == 0, "select: carries first map forward");
  check::ok(m.select(10.0) == 1, "select: exact hit on second map");
  check::ok(m.select(25.0) == 1, "select: carries last map forward");

  // --- validation against a run ---------------------------------------------
  check::ok(occupancy::validate_stiffness_against(m, 2, 1500.0).empty(),
            "matching geometry validates");
  check::ok(!occupancy::validate_stiffness_against(m, 30, 1500.0).empty(),
            "grid mismatch rejected");
  check::ok(!occupancy::validate_stiffness_against(m, 2, 900.0).empty(),
            "scaffold mismatch rejected");
  std::remove(good.c_str());

  // --- shapes that must still load ------------------------------------------
  {
    std::string crlf;
    for (const char *c = kGood; *c; ++c)
    {
      if (*c == '\n')
        crlf += '\r';
      crlf += *c;
    }
    const std::string p = write_tmp("crlf.csv", crlf);
    bool loaded = false;
    try
    {
      loaded = occupancy::load_stiffness_map(p).divs.size() == 2;
    }
    catch (const std::runtime_error &)
    {
    }
    check::ok(loaded, "CRLF line endings accepted");
    std::remove(p.c_str());
  }

  {
    const std::string p =
        write_tmp("comments.csv", "# occupancy-sim stiffness map v1\n"
                                  "grid,2,2\n"
                                  "scaffold_um,1500\n"
                                  "units,kPa\n"
                                  "div_step,5\n"
                                  "div,cx,cy,stiffness\n"
                                  "5,0,0,1.00\n"
                                  "\n"
                                  "# a mid-data comment\n"
                                  "5,1,0,1.25\n"
                                  "5,0,1,1.50\n"
                                  "5,1,1,1.75\n");
    bool loaded = false;
    try
    {
      loaded = occupancy::load_stiffness_map(p).divs.size() == 1;
    }
    catch (const std::runtime_error &)
    {
    }
    check::ok(loaded, "comments/blank lines inside data tolerated");
    std::remove(p.c_str());
  }

  // --- what must be refused --------------------------------------------------
  rejects("noversion.csv",
          "grid,2,2\nscaffold_um,1500\nunits,kPa\ndiv_step,5\n"
          "div,cx,cy,stiffness\n5,0,0,1.0\n",
          "reject: missing version line");

  rejects("badversion.csv",
          "# occupancy-sim stiffness map v99\ngrid,2,2\nscaffold_um,1500\n"
          "units,kPa\ndiv_step,5\ndiv,cx,cy,stiffness\n5,0,0,1.0\n",
          "reject: unknown version");

  rejects("nounits.csv",
          "# occupancy-sim stiffness map v1\ngrid,2,2\nscaffold_um,1500\n"
          "div_step,5\ndiv,cx,cy,stiffness\n5,0,0,1.0\n",
          "reject: missing units");

  rejects("headerfirst.csv",
          "# occupancy-sim stiffness map v1\ndiv,cx,cy,stiffness\ngrid,2,2\n"
          "scaffold_um,1500\nunits,kPa\ndiv_step,5\n5,0,0,1.0\n",
          "reject: data header before grid metadata");

  rejects("dup.csv", std::string(kGood) + "10,1,1,2.75\n",
          "reject: duplicate cell");

  rejects("incomplete.csv",
          "# occupancy-sim stiffness map v1\ngrid,2,2\nscaffold_um,1500\n"
          "units,kPa\ndiv_step,5\ndiv,cx,cy,stiffness\n"
          "5,0,0,1.0\n5,1,0,1.0\n5,0,1,1.0\n",
          "reject: incomplete grid at a DIV");

  rejects("negative.csv",
          "# occupancy-sim stiffness map v1\ngrid,2,2\nscaffold_um,1500\n"
          "units,kPa\ndiv_step,5\ndiv,cx,cy,stiffness\n"
          "5,0,0,1.0\n5,1,0,-3.0\n5,0,1,1.0\n5,1,1,1.0\n",
          "reject: negative stiffness");

  rejects("nan.csv",
          "# occupancy-sim stiffness map v1\ngrid,2,2\nscaffold_um,1500\n"
          "units,kPa\ndiv_step,5\ndiv,cx,cy,stiffness\n"
          "5,0,0,1.0\n5,1,0,nan\n5,0,1,1.0\n5,1,1,1.0\n",
          "reject: NaN stiffness");

  rejects("outside.csv",
          "# occupancy-sim stiffness map v1\ngrid,2,2\nscaffold_um,1500\n"
          "units,kPa\ndiv_step,5\ndiv,cx,cy,stiffness\n5,2,0,1.0\n",
          "reject: compartment outside declared grid");

  rejects("gap.csv",
          "# occupancy-sim stiffness map v1\ngrid,2,2\nscaffold_um,1500\n"
          "units,kPa\ndiv_step,5\ndiv,cx,cy,stiffness\n"
          "5,0,0,1.0\n5,1,0,1.0\n5,0,1,1.0\n5,1,1,1.0\n"
          "15,0,0,2.0\n15,1,0,2.0\n15,0,1,2.0\n15,1,1,2.0\n",
          "reject: gap in the DIV sequence");

  rejects("empty.csv", "", "reject: empty file");

  rejects("headeronly.csv",
          "# occupancy-sim stiffness map v1\ngrid,2,2\nscaffold_um,1500\n"
          "units,kPa\ndiv_step,5\ndiv,cx,cy,stiffness\n",
          "reject: no data rows");

  return check::summary("test_stiffness");
}
