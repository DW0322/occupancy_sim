#include "occupancy/migration.hpp"

#include <stdexcept>

namespace occupancy
{

SeekMode parse_seek_mode(const std::string &s)
{
  if (s == "scaffold")
    return SeekMode::Scaffold;
  if (s == "none")
    return SeekMode::None;
  throw std::invalid_argument("unknown seek mode '" + s +
                              "' (expected none or scaffold)");
}

const char *seek_mode_name(SeekMode m)
{
  return (m == SeekMode::Scaffold) ? "scaffold" : "none";
}

CohereGate parse_cohere_gate(const std::string &s)
{
  if (s == "scaffold")
    return CohereGate::Scaffold;
  if (s == "none")
    return CohereGate::None;
  throw std::invalid_argument("unknown cohere gate '" + s +
                              "' (expected none or scaffold)");
}

const char *cohere_gate_name(CohereGate g)
{
  return (g == CohereGate::Scaffold) ? "scaffold" : "none";
}

} // namespace occupancy
