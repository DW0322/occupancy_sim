#ifndef OCCUPANCY_SCHEDULE_HPP
#define OCCUPANCY_SCHEDULE_HPP
// The DIV reveal schedule: how much of an arbor is visible on a given day.

namespace occupancy
{

/**
 * @brief Fraction of an arbor revealed at a given DIV.
 *
 * @param div_days Days in vitro.
 * @return Fraction 0-1 of each arbor's max path distance that is revealed.
 */
double div_to_fraction(double div_days);

/**
 * @brief The DIV where div_to_fraction() first hits 1.0.
 * @return The saturation DIV (15).
 */
double reveal_saturation_div();

} // namespace occupancy

#endif // OCCUPANCY_SCHEDULE_HPP
