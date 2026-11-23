#include "check.hpp"
#include "occupancy/distance_transform.hpp"
#include <random>
#include <vector>

int main()
{
  std::cout << "test_distance_transform\n";
  const int G = 48;
  std::mt19937 rng(7);
  std::uniform_int_distribution<int> coord(0, G - 1);

  for (double R : {1.0, 3.5, 7.0})
  {
    std::vector<uint8_t> feat((size_t)G * G, 0);
    std::vector<std::pair<int, int>> pts;
    for (int k = 0; k < 12; ++k)
    {
      int x = coord(rng), y = coord(rng);
      feat[(size_t)y * G + x] = 1;
      pts.push_back({x, y});
    }
    std::vector<uint8_t> mask((size_t)G * G, 0);
    occupancy::dilate_by_edt(feat, mask, G, R);

    // brute force: a pixel is covered iff some feature pixel is within R
    bool same = true;
    for (int y = 0; y < G && same; ++y)
      for (int x = 0; x < G && same; ++x)
      {
        bool want = false;
        for (auto &p : pts)
        {
          double dx = x - p.first, dy = y - p.second;
          if (dx * dx + dy * dy <= R * R)
          {
            want = true;
            break;
          }
        }
        if (want != (mask[(size_t)y * G + x] != 0))
          same = false;
      }
    check::ok(same,
              "EDT dilation matches brute force at R=" + std::to_string(R));
  }

  // the union property the accumulate-then-dilate-once design rests on
  std::vector<uint8_t> A((size_t)G * G, 0), B((size_t)G * G, 0),
      U((size_t)G * G, 0);
  for (int k = 0; k < 8; ++k)
  {
    int x = coord(rng), y = coord(rng);
    A[(size_t)y * G + x] = 1;
    U[(size_t)y * G + x] = 1;
    int x2 = coord(rng), y2 = coord(rng);
    B[(size_t)y2 * G + x2] = 1;
    U[(size_t)y2 * G + x2] = 1;
  }
  std::vector<uint8_t> dA((size_t)G * G, 0), dB((size_t)G * G, 0),
      dU((size_t)G * G, 0);
  occupancy::dilate_by_edt(A, dA, G, 4.0);
  occupancy::dilate_by_edt(B, dB, G, 4.0);
  occupancy::dilate_by_edt(U, dU, G, 4.0);
  bool distributes = true;
  for (size_t i = 0; i < dU.size(); ++i)
    if ((dA[i] || dB[i]) != (dU[i] != 0))
      distributes = false;
  check::ok(distributes, "dilate(A u B) == dilate(A) u dilate(B)");

  return check::summary("test_distance_transform");
}
