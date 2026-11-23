#include "occupancy/schedule.hpp"
#include <cstddef>

namespace occupancy
{

namespace
{
const double kDiv[] = {0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 12.0, 15.0};
const double kFrac[] = {0.02, 0.06, 0.10, 0.16, 0.30, 0.58, 0.76, 0.90, 1.00};
const int kN = (int)(sizeof(kDiv) / sizeof(kDiv[0]));
} // namespace

double div_to_fraction(double div_days)
{
  if (div_days <= kDiv[0])
    return kFrac[0];
  if (div_days >= kDiv[kN - 1])
    return kFrac[kN - 1];
  for (int i = 1; i < kN; ++i)
    if (div_days <= kDiv[i])
    {
      double t = (div_days - kDiv[i - 1]) / (kDiv[i] - kDiv[i - 1]);
      return kFrac[i - 1] + t * (kFrac[i] - kFrac[i - 1]);
    }
  return kFrac[kN - 1];
}

double reveal_saturation_div() { return kDiv[kN - 1]; }

} // namespace occupancy
