#ifndef OCCUPANCY_TESTS_CHECK_HPP
#define OCCUPANCY_TESTS_CHECK_HPP
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

// lil helping checker for me self :)
namespace check
{
inline int failures = 0;

inline void ok(bool cond, const std::string &what)
{
  if (cond)
    std::cout << "  pass  " << what << "\n";
  else
  {
    std::cout << "  FAIL  " << what << "\n";
    ++failures;
  }
}

inline void near(double got, double want, double tol, const std::string &what)
{
  bool good = std::fabs(got - want) <= tol;
  std::cout << (good ? "  pass  " : "  FAIL  ") << what << "  got " << got
            << " want " << want << " tol " << tol << "\n";
  if (!good)
    ++failures;
}

inline int summary(const std::string &name)
{
  std::cout << name << ": " << (failures ? "FAILED" : "ok") << "\n";
  return failures ? 1 : 0;
}
} // namespace check

#endif // OCCUPANCY_TESTS_CHECK_HPP
