#ifndef RANDOM_HPP
#define RANDOM_HPP
#include <cstdint>
#include <cmath>
#include <random>

class pcg32_random_t
{
public:
    pcg32_random_t();
    uint32_t rnd_integer();
    float rnd();
private:
    uint64_t state;
    uint64_t inc;
};

class RandomManager
{
public:
    RandomManager();
    float GetRND();
private:
    pcg32_random_t rnd_class;
};

#endif