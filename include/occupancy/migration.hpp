#ifndef OCCUPANCY_MIGRATION_HPP
#define OCCUPANCY_MIGRATION_HPP
// Seek modes: optional biases layered on top of the soma random walk.
#include <string>

namespace occupancy
{

/**
 * @brief What a seeking cell steers toward.
 *
 * Scaffold is the decision rule the collaboration asked for (seek the
 * scaffold or carry on migrating), and has two halves off the same geometry:
 *
 * - OFF it, steer at the nearest scaffold compartment. Haptotaxis: somal
 *   translocation in dissociated cortical culture is driven by an adhesion
 *   difference (Biomaterials 2018), and the microstructure is the target.
 * - ON it, follow the strand's local tangent rather than going isotropic
 *   again. Contact guidance: Li 2025 attributes the M-Alg advantage to
 *   topographical cues and reports networks throughout, not cells piled where
 *   they first touched, which is what steering at the nearest material gives.
 *
 * Li's contrast is P-Alg vs M-Alg, two SCAFFOLDS, not the lit and unlit parts
 * of one field, so reading it as a within-field gradient is our
 * extrapolation. It sets the cue's direction, not its spatial form.
 *
 * Params::p_seek picks which cells apply the rule, Params::seek_weight how hard
 * it steers, so "prefer the scaffold" is swept rather than switched.
 */
enum class SeekMode
{
  None,    ///< No bias; plain random walk.
  Scaffold ///< Haptotaxis: steer toward the loaded scaffold geometry.
};

/**
 * @brief Parse a seek mode from the CLI.
 *
 * Strict: a silent fallback turned `--seek scafold` into a clean-looking null
 * run that read exactly like "the rule does nothing".
 *
 * @param s "none" or "scaffold".
 * @return The mode.
 * @throws std::invalid_argument If s is neither.
 */
SeekMode parse_seek_mode(const std::string &s);

/**
 * @brief What decides whether a cell coheres with its neighbours.
 *
 * Aggregation in culture is a competition: cells reaggregate where
 * cell-substrate adhesion is weak and spread where it is strong. Gating on the
 * geometry makes cohesion the complement of seek rather than a second knob
 * beside it. Off by default, so a loaded geometry stays inert until switched
 * on.
 */
enum class CohereGate
{
  None,    ///< Cohesion applies everywhere at full weight.
  Scaffold ///< Only cells off the scaffold cohere.
};

/**
 * @brief Parse a cohesion gate from the CLI. Strict, as parse_seek_mode().
 * @param s "none" or "scaffold".
 * @return The gate.
 * @throws std::invalid_argument If s is neither.
 */
CohereGate parse_cohere_gate(const std::string &s);

/**
 * @brief Inverse of parse_cohere_gate().
 * @param g The gate.
 * @return "scaffold" or "none".
 */
const char *cohere_gate_name(CohereGate g);

/**
 * @brief Inverse of parse_seek_mode().
 * @param m The mode.
 * @return "scaffold" or "none".
 */
const char *seek_mode_name(SeekMode m);

} // namespace occupancy

#endif // OCCUPANCY_MIGRATION_HPP
