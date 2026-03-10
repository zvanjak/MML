///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        random_tests.cpp                                                    ///
///  Description: Comprehensive tests for Random.h                                    ///
///               Tests for uniform distributions and random direction generation     ///
///                                                                                   ///
///  Coverage:    Random::UniformReal                                                 ///
///               Random::UniformInt                                                  ///
///               Random::UniformVecDirection2                                        ///
///               Random::UniformVecDirection3                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "../TestPrecision.h"
#include "../../mml/base/Random.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace MML;
using Catch::Approx;

namespace MML::Tests::Base::RandomTests {

///////////////////////////////////////////////////////////////////////////////////////////
///                         UniformReal Tests                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Random::UniformReal - Basic range tests", "[Random][UniformReal]")
{
    SECTION("Values within [0, 1] range")
    {
        constexpr int NUM_SAMPLES = 1000;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real val = Random::UniformReal(0.0, 1.0);
            REQUIRE(val >= 0.0);
            REQUIRE(val < 1.0);  // uniform_real_distribution is [min, max)
        }
    }

    SECTION("Values within negative range [-10, -5]")
    {
        constexpr int NUM_SAMPLES = 1000;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real val = Random::UniformReal(-10.0, -5.0);
            REQUIRE(val >= -10.0);
            REQUIRE(val < -5.0);
        }
    }

    SECTION("Values within mixed range [-5, 5]")
    {
        constexpr int NUM_SAMPLES = 1000;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real val = Random::UniformReal(-5.0, 5.0);
            REQUIRE(val >= -5.0);
            REQUIRE(val < 5.0);
        }
    }

    SECTION("Large range [0, 1e6]")
    {
        constexpr int NUM_SAMPLES = 1000;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real val = Random::UniformReal(0.0, 1e6);
            REQUIRE(val >= 0.0);
            REQUIRE(val < 1e6);
        }
    }

    SECTION("Small range [0.001, 0.002]")
    {
        constexpr int NUM_SAMPLES = 1000;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real val = Random::UniformReal(0.001, 0.002);
            REQUIRE(val >= 0.001);
            REQUIRE(val < 0.002);
        }
    }
}

TEST_CASE("Random::UniformReal - Distribution uniformity", "[Random][UniformReal][Statistics]")
{
    SECTION("Values spread across range [0, 10] - bucket test")
    {
        constexpr int NUM_SAMPLES = 10000;
        constexpr int NUM_BUCKETS = 10;
        std::vector<int> buckets(NUM_BUCKETS, 0);

        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real val = Random::UniformReal(0.0, 10.0);
            int bucket = static_cast<int>(val);  // 0-9
            if (bucket >= 0 && bucket < NUM_BUCKETS)
                buckets[bucket]++;
        }

        // Each bucket should have roughly NUM_SAMPLES/NUM_BUCKETS = 1000
        // Allow 30% deviation for statistical variance
        int expected = NUM_SAMPLES / NUM_BUCKETS;
        int tolerance = expected * 30 / 100;  // 30% tolerance

        for (int i = 0; i < NUM_BUCKETS; ++i)
        {
            INFO("Bucket " << i << " has " << buckets[i] << " samples, expected ~" << expected);
            REQUIRE(buckets[i] >= expected - tolerance);
            REQUIRE(buckets[i] <= expected + tolerance);
        }
    }

    SECTION("Mean and variance approximately correct")
    {
        constexpr int NUM_SAMPLES = 10000;
        Real min = 0.0, max = 10.0;
        Real expected_mean = (min + max) / 2.0;  // 5.0
        Real expected_variance = (max - min) * (max - min) / 12.0;  // 8.333...

        std::vector<Real> samples(NUM_SAMPLES);
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            samples[i] = Random::UniformReal(min, max);
        }

        // Calculate sample mean
        Real sum = std::accumulate(samples.begin(), samples.end(), 0.0);
        Real mean = sum / NUM_SAMPLES;

        // Calculate sample variance
        Real variance_sum = 0.0;
        for (Real s : samples)
        {
            variance_sum += (s - mean) * (s - mean);
        }
        Real variance = variance_sum / (NUM_SAMPLES - 1);

        // Allow 5% deviation
        REQUIRE(mean == Approx(expected_mean).epsilon(0.05));
        REQUIRE(variance == Approx(expected_variance).epsilon(0.10));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         UniformInt Tests                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Random::UniformInt - Basic range tests", "[Random][UniformInt]")
{
    SECTION("Values within [0, 10] range (inclusive)")
    {
        constexpr int NUM_SAMPLES = 1000;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            int val = Random::UniformInt(0, 10);
            REQUIRE(val >= 0);
            REQUIRE(val <= 10);  // uniform_int_distribution is [min, max] inclusive
        }
    }

    SECTION("Values within negative range [-10, -5]")
    {
        constexpr int NUM_SAMPLES = 1000;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            int val = Random::UniformInt(-10, -5);
            REQUIRE(val >= -10);
            REQUIRE(val <= -5);
        }
    }

    SECTION("Single value range [5, 5]")
    {
        constexpr int NUM_SAMPLES = 100;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            int val = Random::UniformInt(5, 5);
            REQUIRE(val == 5);
        }
    }

    SECTION("Two value range [0, 1]")
    {
        constexpr int NUM_SAMPLES = 1000;
        int zeros = 0, ones = 0;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            int val = Random::UniformInt(0, 1);
            REQUIRE((val == 0 || val == 1));
            if (val == 0) zeros++;
            else ones++;
        }
        // Both should appear roughly equally (within 40% of expected 500)
        REQUIRE(zeros >= 300);
        REQUIRE(zeros <= 700);
    }
}

TEST_CASE("Random::UniformInt - Boundary inclusion", "[Random][UniformInt]")
{
    SECTION("Both min and max are reachable")
    {
        constexpr int NUM_SAMPLES = 10000;
        int min_val = 0, max_val = 5;
        bool found_min = false, found_max = false;

        for (int i = 0; i < NUM_SAMPLES && !(found_min && found_max); ++i)
        {
            int val = Random::UniformInt(min_val, max_val);
            if (val == min_val) found_min = true;
            if (val == max_val) found_max = true;
        }

        REQUIRE(found_min);
        REQUIRE(found_max);
    }
}

TEST_CASE("Random::UniformInt - Distribution uniformity", "[Random][UniformInt][Statistics]")
{
    SECTION("Each value appears roughly equally [0, 9]")
    {
        constexpr int NUM_SAMPLES = 10000;
        constexpr int RANGE = 10;  // 0-9
        std::vector<int> counts(RANGE, 0);

        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            int val = Random::UniformInt(0, RANGE - 1);
            counts[val]++;
        }

        // Each value should appear roughly NUM_SAMPLES/RANGE = 1000 times
        int expected = NUM_SAMPLES / RANGE;
        int tolerance = expected * 30 / 100;  // 30% tolerance

        for (int i = 0; i < RANGE; ++i)
        {
            INFO("Value " << i << " appeared " << counts[i] << " times, expected ~" << expected);
            REQUIRE(counts[i] >= expected - tolerance);
            REQUIRE(counts[i] <= expected + tolerance);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         UniformVecDirection2 Tests                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Random::UniformVecDirection2 - Magnitude tests", "[Random][UniformVecDirection2]")
{
    SECTION("Unit vector (abs = 1.0)")
    {
        constexpr int NUM_SAMPLES = 100;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy;
            Real result = Random::UniformVecDirection2(vx, vy, 1.0);
            
            REQUIRE(result == Approx(1.0));
            Real magnitude = std::sqrt(vx * vx + vy * vy);
            REQUIRE(magnitude == Approx(1.0).epsilon(1e-10));
        }
    }

    SECTION("Custom magnitude (abs = 5.0)")
    {
        constexpr int NUM_SAMPLES = 100;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy;
            Real result = Random::UniformVecDirection2(vx, vy, 5.0);
            
            REQUIRE(result == Approx(5.0));
            Real magnitude = std::sqrt(vx * vx + vy * vy);
            REQUIRE(magnitude == Approx(5.0).epsilon(1e-10));
        }
    }

    SECTION("Zero magnitude (abs = 0.0)")
    {
        Real vx, vy;
        Real result = Random::UniformVecDirection2(vx, vy, 0.0);
        
        REQUIRE(result == Approx(0.0));
        REQUIRE(vx == Approx(0.0).margin(1e-15));
        REQUIRE(vy == Approx(0.0).margin(1e-15));
    }

    SECTION("Small magnitude (abs = 1e-10)")
    {
        Real vx, vy;
        Real result = Random::UniformVecDirection2(vx, vy, 1e-10);
        
        REQUIRE(result == Approx(1e-10));
        Real magnitude = std::sqrt(vx * vx + vy * vy);
        REQUIRE(magnitude == Approx(1e-10).epsilon(1e-5));
    }

    SECTION("Large magnitude (abs = 1e6)")
    {
        Real vx, vy;
        Real result = Random::UniformVecDirection2(vx, vy, 1e6);
        
        REQUIRE(result == Approx(1e6));
        Real magnitude = std::sqrt(vx * vx + vy * vy);
        REQUIRE(magnitude == Approx(1e6).epsilon(1e-10));
    }
}

TEST_CASE("Random::UniformVecDirection2 - Angular distribution", "[Random][UniformVecDirection2][Statistics]")
{
    SECTION("Angles uniformly distributed in [0, 2π]")
    {
        constexpr int NUM_SAMPLES = 10000;
        constexpr int NUM_BUCKETS = 8;  // 8 sectors of 45 degrees each
        std::vector<int> buckets(NUM_BUCKETS, 0);

        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy;
            Random::UniformVecDirection2(vx, vy, 1.0);
            
            // Calculate angle in [0, 2π)
            Real angle = std::atan2(vy, vx);
            if (angle < 0) angle += 2 * Constants::PI;
            
            int bucket = static_cast<int>(angle / (2 * Constants::PI / NUM_BUCKETS));
            if (bucket >= NUM_BUCKETS) bucket = NUM_BUCKETS - 1;
            buckets[bucket]++;
        }

        // Each bucket should have roughly NUM_SAMPLES/NUM_BUCKETS = 1250
        int expected = NUM_SAMPLES / NUM_BUCKETS;
        int tolerance = expected * 30 / 100;  // 30% tolerance

        for (int i = 0; i < NUM_BUCKETS; ++i)
        {
            INFO("Angular sector " << i << " has " << buckets[i] << " samples, expected ~" << expected);
            REQUIRE(buckets[i] >= expected - tolerance);
            REQUIRE(buckets[i] <= expected + tolerance);
        }
    }
}

TEST_CASE("Random::UniformVecDirection2 - All quadrants covered", "[Random][UniformVecDirection2]")
{
    SECTION("Vectors appear in all four quadrants")
    {
        constexpr int NUM_SAMPLES = 1000;
        bool q1 = false, q2 = false, q3 = false, q4 = false;

        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy;
            Random::UniformVecDirection2(vx, vy, 1.0);
            
            if (vx > 0 && vy > 0) q1 = true;       // Quadrant I
            else if (vx < 0 && vy > 0) q2 = true;  // Quadrant II
            else if (vx < 0 && vy < 0) q3 = true;  // Quadrant III
            else if (vx > 0 && vy < 0) q4 = true;  // Quadrant IV
        }

        REQUIRE(q1);
        REQUIRE(q2);
        REQUIRE(q3);
        REQUIRE(q4);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         UniformVecDirection3 Tests                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Random::UniformVecDirection3 - Magnitude tests", "[Random][UniformVecDirection3]")
{
    SECTION("Unit vector (abs = 1.0)")
    {
        constexpr int NUM_SAMPLES = 100;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy, vz;
            Real result = Random::UniformVecDirection3(vx, vy, vz, 1.0);
            
            REQUIRE(result == Approx(1.0));
            Real magnitude = std::sqrt(vx * vx + vy * vy + vz * vz);
            REQUIRE(magnitude == Approx(1.0).epsilon(1e-10));
        }
    }

    SECTION("Custom magnitude (abs = 7.5)")
    {
        constexpr int NUM_SAMPLES = 100;
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy, vz;
            Real result = Random::UniformVecDirection3(vx, vy, vz, 7.5);
            
            REQUIRE(result == Approx(7.5));
            Real magnitude = std::sqrt(vx * vx + vy * vy + vz * vz);
            REQUIRE(magnitude == Approx(7.5).epsilon(1e-10));
        }
    }

    SECTION("Zero magnitude (abs = 0.0)")
    {
        Real vx, vy, vz;
        Real result = Random::UniformVecDirection3(vx, vy, vz, 0.0);
        
        REQUIRE(result == Approx(0.0));
        REQUIRE(vx == Approx(0.0).margin(1e-15));
        REQUIRE(vy == Approx(0.0).margin(1e-15));
        REQUIRE(vz == Approx(0.0).margin(1e-15));
    }

    SECTION("Small magnitude (abs = 1e-10)")
    {
        Real vx, vy, vz;
        Real result = Random::UniformVecDirection3(vx, vy, vz, 1e-10);
        
        REQUIRE(result == Approx(1e-10));
        Real magnitude = std::sqrt(vx * vx + vy * vy + vz * vz);
        REQUIRE(magnitude == Approx(1e-10).epsilon(1e-5));
    }

    SECTION("Large magnitude (abs = 1e6)")
    {
        Real vx, vy, vz;
        Real result = Random::UniformVecDirection3(vx, vy, vz, 1e6);
        
        REQUIRE(result == Approx(1e6));
        Real magnitude = std::sqrt(vx * vx + vy * vy + vz * vz);
        REQUIRE(magnitude == Approx(1e6).epsilon(1e-10));
    }
}

TEST_CASE("Random::UniformVecDirection3 - All octants covered", "[Random][UniformVecDirection3]")
{
    SECTION("Vectors appear in all eight octants")
    {
        constexpr int NUM_SAMPLES = 2000;
        std::vector<bool> octants(8, false);

        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy, vz;
            Random::UniformVecDirection3(vx, vy, vz, 1.0);
            
            // Octant index based on signs: (vx>0)*4 + (vy>0)*2 + (vz>0)*1
            int octant = (vx > 0 ? 4 : 0) + (vy > 0 ? 2 : 0) + (vz > 0 ? 1 : 0);
            octants[octant] = true;
        }

        for (int i = 0; i < 8; ++i)
        {
            INFO("Octant " << i << " covered: " << octants[i]);
            REQUIRE(octants[i]);
        }
    }
}

TEST_CASE("Random::UniformVecDirection3 - Spherical distribution uniformity", "[Random][UniformVecDirection3][Statistics]")
{
    SECTION("Z-component uniformly distributed in [-1, 1] (cos(theta))")
    {
        // For a uniform distribution on a sphere, cos(theta) should be uniform in [-1, 1]
        // This is the key test for proper spherical distribution (not biased toward poles)
        constexpr int NUM_SAMPLES = 10000;
        constexpr int NUM_BUCKETS = 10;  // Divide [-1, 1] into 10 buckets
        std::vector<int> z_buckets(NUM_BUCKETS, 0);

        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy, vz;
            Random::UniformVecDirection3(vx, vy, vz, 1.0);
            
            // vz is in [-1, 1] for unit sphere
            int bucket = static_cast<int>((vz + 1.0) / 2.0 * NUM_BUCKETS);
            if (bucket >= NUM_BUCKETS) bucket = NUM_BUCKETS - 1;
            if (bucket < 0) bucket = 0;
            z_buckets[bucket]++;
        }

        // Each bucket should have roughly NUM_SAMPLES/NUM_BUCKETS = 1000
        int expected = NUM_SAMPLES / NUM_BUCKETS;
        int tolerance = expected * 30 / 100;  // 30% tolerance

        for (int i = 0; i < NUM_BUCKETS; ++i)
        {
            INFO("Z-bucket " << i << " (z in [" << (-1.0 + i * 0.2) << ", " << (-1.0 + (i+1) * 0.2) << "]) has " 
                 << z_buckets[i] << " samples, expected ~" << expected);
            REQUIRE(z_buckets[i] >= expected - tolerance);
            REQUIRE(z_buckets[i] <= expected + tolerance);
        }
    }

    SECTION("Azimuthal angle (phi) uniformly distributed")
    {
        constexpr int NUM_SAMPLES = 10000;
        constexpr int NUM_BUCKETS = 8;  // 8 sectors of 45 degrees
        std::vector<int> phi_buckets(NUM_BUCKETS, 0);

        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy, vz;
            Random::UniformVecDirection3(vx, vy, vz, 1.0);
            
            // Calculate azimuthal angle from x-y projection
            Real phi = std::atan2(vy, vx);
            if (phi < 0) phi += 2 * Constants::PI;
            
            int bucket = static_cast<int>(phi / (2 * Constants::PI / NUM_BUCKETS));
            if (bucket >= NUM_BUCKETS) bucket = NUM_BUCKETS - 1;
            phi_buckets[bucket]++;
        }

        // Each bucket should have roughly NUM_SAMPLES/NUM_BUCKETS = 1250
        int expected = NUM_SAMPLES / NUM_BUCKETS;
        int tolerance = expected * 30 / 100;  // 30% tolerance

        for (int i = 0; i < NUM_BUCKETS; ++i)
        {
            INFO("Phi sector " << i << " has " << phi_buckets[i] << " samples, expected ~" << expected);
            REQUIRE(phi_buckets[i] >= expected - tolerance);
            REQUIRE(phi_buckets[i] <= expected + tolerance);
        }
    }
}

TEST_CASE("Random::UniformVecDirection3 - Hemisphere balance", "[Random][UniformVecDirection3][Statistics]")
{
    SECTION("Upper and lower hemispheres equally populated")
    {
        constexpr int NUM_SAMPLES = 10000;
        int upper = 0, lower = 0;

        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            Real vx, vy, vz;
            Random::UniformVecDirection3(vx, vy, vz, 1.0);
            
            if (vz > 0) upper++;
            else if (vz < 0) lower++;
            // vz == 0 is rare, ignore
        }

        // Should be roughly 50/50
        int expected = NUM_SAMPLES / 2;
        int tolerance = expected * 10 / 100;  // 10% tolerance

        REQUIRE(upper >= expected - tolerance);
        REQUIRE(upper <= expected + tolerance);
        REQUIRE(lower >= expected - tolerance);
        REQUIRE(lower <= expected + tolerance);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         Thread Safety Tests                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Random - Thread-local RNG ensures reproducible behavior", "[Random][ThreadSafety]")
{
    SECTION("Multiple calls produce different values")
    {
        // This verifies the RNG is working (not stuck)
        constexpr int NUM_SAMPLES = 100;
        std::vector<Real> values(NUM_SAMPLES);
        
        for (int i = 0; i < NUM_SAMPLES; ++i)
        {
            values[i] = Random::UniformReal(0.0, 1.0);
        }

        // Check that not all values are the same
        Real first = values[0];
        bool all_same = true;
        for (int i = 1; i < NUM_SAMPLES; ++i)
        {
            if (values[i] != first)
            {
                all_same = false;
                break;
            }
        }
        REQUIRE_FALSE(all_same);
    }
}

} // namespace MML::Tests::Base::RandomTests
