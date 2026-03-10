/******************************************************************************
 * MML Example: Formula 1 G-Force Analysis
 * ============================================================================
 * 
 * Demonstrates parametric curve analysis using real F1 track data:
 * 
 *   1. Load telemetry data (position, speed, time)
 *   2. Create parametric curve from track coordinates
 *   3. Calculate curvature along the track
 *   4. Compute G-forces (lateral + longitudinal)
 *   5. Visualize the results
 * 
 * Physics:
 *   - Lateral G = v² · κ / g  (cornering)
 *   - Longitudinal G = (1/g) · dv/dt  (accel/braking)
 *   - Total G = sqrt(G_lat² + G_long²)
 * 
 * Data: Default is Silverstone (~5.8 km, ~89 seconds)
 *       Use --track=monza or -m for Monza (~5.7 km, ~81 seconds)
 * 
 * Build: cmake --build build --target Example03_F1GForces
 * Run:   ./build/src/examples/Release/Example03_F1GForces
 * 
 *****************************************************************************/

#include "MMLBase.h"
#include "mml/base/Vector/Vector.h"
#include "mml/base/Matrix/Matrix.h"
#include "mml/base/InterpolatedFunction.h"
#include "mml/core/Curves.h"
#include "mml/tools/Visualizer.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace MML;

/******************************************************************************
 * DATA STRUCTURES
 *****************************************************************************/

/// @brief Single telemetry sample from the CSV
struct TelemetrySample {
    Real time_s;      ///< Time in seconds
    Real distance_m;  ///< Distance traveled in meters
    Real x_m;         ///< X position in meters
    Real y_m;         ///< Y position in meters
    Real speed_kmh;   ///< Speed in km/h
    
    /// @brief Convert speed to m/s
    Real speed_ms() const { return speed_kmh / 3.6; }
};

/// @brief Complete lap telemetry data
struct LapTelemetry {
    std::vector<TelemetrySample> samples;
    
    // Statistics
    Real minSpeed_kmh = 1e10;
    Real maxSpeed_kmh = -1e10;
    Real totalDistance_m = 0;
    Real totalTime_s = 0;
    
    /// @brief Number of samples
    size_t size() const { return samples.size(); }
    
    /// @brief Access sample by index
    const TelemetrySample& operator[](size_t i) const { return samples[i]; }
    
    /// @brief Calculate statistics after loading
    void computeStatistics() {
        if (samples.empty()) return;
        
        for (const auto& s : samples) {
            minSpeed_kmh = std::min(minSpeed_kmh, s.speed_kmh);
            maxSpeed_kmh = std::max(maxSpeed_kmh, s.speed_kmh);
        }
        totalDistance_m = samples.back().distance_m;
        totalTime_s = samples.back().time_s;
    }
    
    /// @brief Print summary
    void printSummary(const std::string& trackTitle = "Unknown") const {
        std::string header = trackTitle + " LAP TELEMETRY SUMMARY";
        // Center the header in 56-char field
        int padding = (56 - (int)header.size()) / 2;
        if (padding < 0) padding = 0;
        std::cout << "\n";
        std::cout << "╔══════════════════════════════════════════════════════════╗\n";
        std::cout << "║" << std::string(padding, ' ') << std::left << std::setw(56 - padding) << header << "║\n";
        std::cout << "╠══════════════════════════════════════════════════════════╣\n";
        std::cout << "║  Samples:        " << std::setw(8) << samples.size() << "                              ║\n";
        std::cout << "║  Lap Time:       " << std::setw(8) << std::fixed << std::setprecision(2) 
                  << totalTime_s << " s                            ║\n";
        std::cout << "║  Track Length:   " << std::setw(8) << std::fixed << std::setprecision(1) 
                  << totalDistance_m << " m  (" << std::setprecision(2) << totalDistance_m/1000.0 << " km)          ║\n";
        std::cout << "║  Min Speed:      " << std::setw(8) << std::fixed << std::setprecision(1) 
                  << minSpeed_kmh << " km/h                        ║\n";
        std::cout << "║  Max Speed:      " << std::setw(8) << std::fixed << std::setprecision(1) 
                  << maxSpeed_kmh << " km/h                        ║\n";
        std::cout << "║  Avg Speed:      " << std::setw(8) << std::fixed << std::setprecision(1) 
                  << (totalDistance_m / totalTime_s) * 3.6 << " km/h                        ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    }
};

/******************************************************************************
 * CSV LOADING
 *****************************************************************************/

/// @brief Load telemetry data from CSV file
/// @param filename Path to CSV file
/// @return LapTelemetry with all samples loaded
LapTelemetry loadTelemetryCSV(const std::string& filename) {
    LapTelemetry telemetry;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open file: " << filename << std::endl;
        return telemetry;
    }
    
    std::string line;
    int lineNum = 0;
    
    while (std::getline(file, line)) {
        lineNum++;
        
        // Skip header line
        if (lineNum == 1) {
            std::cout << "Header: " << line << std::endl;
            continue;
        }
        
        // Skip empty lines
        if (line.empty()) continue;
        
        // Parse CSV line
        std::stringstream ss(line);
        std::string token;
        TelemetrySample sample;
        
        try {
            std::getline(ss, token, ','); sample.time_s = std::stod(token);
            std::getline(ss, token, ','); sample.distance_m = std::stod(token);
            std::getline(ss, token, ','); sample.x_m = std::stod(token);
            std::getline(ss, token, ','); sample.y_m = std::stod(token);
            std::getline(ss, token, ','); sample.speed_kmh = std::stod(token);
            
            telemetry.samples.push_back(sample);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to parse line " << lineNum << ": " << e.what() << std::endl;
        }
    }
    
    file.close();
    telemetry.computeStatistics();
    
    std::cout << "Loaded " << telemetry.size() << " telemetry samples from " << filename << std::endl;
    
    return telemetry;
}

/******************************************************************************
 * MAIN - Step 1: Load and verify data
 *****************************************************************************/

int main(int argc, char* argv[]) {
    // Parse command line args: --track=monza or --track=silverstone
    std::string trackName = "silverstone";  // Default
    std::string trackFile = "f1_silverstone_lap.csv";
    std::string trackTitle = "Silverstone";
    std::string driverInfo = "Lewis Hamilton";
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--track=silverstone" || arg == "silverstone" || arg == "-s") {
            trackName = "silverstone";
            trackFile = "f1_silverstone_lap.csv";
            trackTitle = "Silverstone";
            driverInfo = "Lewis Hamilton";
        } else if (arg == "--track=monza" || arg == "monza" || arg == "-m") {
            trackName = "monza";
            trackFile = "f1_monza_lap.csv";
            trackTitle = "Monza";
            driverInfo = "Lando Norris";
        }
    }
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║     MML EXAMPLE: FORMULA 1 G-FORCE ANALYSIS                  ║\n";
    std::cout << "║     Track: " << std::left << std::setw(20) << trackTitle 
              << "Driver: " << std::setw(20) << driverInfo << "║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
    
    // =========================================================================
    // STEP 1: Load telemetry data
    // =========================================================================
    std::cout << "\n=== STEP 1: Loading Telemetry Data ===\n\n";
    
    // Try multiple possible paths (depending on where we run from)
    std::vector<std::string> possiblePaths = {
        "src/examples/03_formula_1_sim/data/" + trackFile,
        "../src/examples/03_formula_1_sim/data/" + trackFile,
        "../../src/examples/03_formula_1_sim/data/" + trackFile,
        "data/" + trackFile
    };
    
    LapTelemetry telemetry;
    for (const auto& path : possiblePaths) {
        telemetry = loadTelemetryCSV(path);
        if (!telemetry.samples.empty()) {
            break;
        }
    }
    
    if (telemetry.samples.empty()) {
        std::cerr << "ERROR: Could not load telemetry data from any path!\n";
        return 1;
    }
    
    // Print summary
    telemetry.printSummary(trackTitle);
    
    // =========================================================================
    // STEP 2: Show some sample data points
    // =========================================================================
    std::cout << "=== Sample Data Points ===\n\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Index |   Time(s) | Distance(m) |      X(m) |      Y(m) | Speed(km/h)\n";
    std::cout << "  ------|-----------|-------------|-----------|-----------|------------\n";
    
    // First 5 points
    for (size_t i = 0; i < 5 && i < telemetry.size(); i++) {
        const auto& s = telemetry[i];
        std::cout << "  " << std::setw(5) << i << " | "
                  << std::setw(9) << s.time_s << " | "
                  << std::setw(11) << s.distance_m << " | "
                  << std::setw(9) << s.x_m << " | "
                  << std::setw(9) << s.y_m << " | "
                  << std::setw(10) << s.speed_kmh << "\n";
    }
    std::cout << "  ...   |    ...    |     ...     |    ...    |    ...    |    ...\n";
    
    // Last 5 points
    for (size_t i = telemetry.size() - 5; i < telemetry.size(); i++) {
        const auto& s = telemetry[i];
        std::cout << "  " << std::setw(5) << i << " | "
                  << std::setw(9) << s.time_s << " | "
                  << std::setw(11) << s.distance_m << " | "
                  << std::setw(9) << s.x_m << " | "
                  << std::setw(9) << s.y_m << " | "
                  << std::setw(10) << s.speed_kmh << "\n";
    }
    
    // =========================================================================
    // STEP 3: Find interesting points (min/max speed)
    // =========================================================================
    std::cout << "\n=== Key Track Points ===\n\n";
    
    // Find min and max speed locations
    size_t minSpeedIdx = 0, maxSpeedIdx = 0;
    for (size_t i = 0; i < telemetry.size(); i++) {
        if (telemetry[i].speed_kmh < telemetry[minSpeedIdx].speed_kmh) minSpeedIdx = i;
        if (telemetry[i].speed_kmh > telemetry[maxSpeedIdx].speed_kmh) maxSpeedIdx = i;
    }
    
    const auto& slowest = telemetry[minSpeedIdx];
    const auto& fastest = telemetry[maxSpeedIdx];
    
    std::cout << "  Slowest point (likely chicane/corner):\n";
    std::cout << "    Index: " << minSpeedIdx << ", Time: " << slowest.time_s << "s\n";
    std::cout << "    Position: (" << slowest.x_m << ", " << slowest.y_m << ") m\n";
    std::cout << "    Speed: " << slowest.speed_kmh << " km/h\n\n";
    
    std::cout << "  Fastest point (likely main straight):\n";
    std::cout << "    Index: " << maxSpeedIdx << ", Time: " << fastest.time_s << "s\n";
    std::cout << "    Position: (" << fastest.x_m << ", " << fastest.y_m << ") m\n";
    std::cout << "    Speed: " << fastest.speed_kmh << " km/h\n\n";
    
    std::cout << "=== Step 1 Complete! Data loaded successfully. ===\n\n";
    
    // =========================================================================
    // STEP 2: Create parametric curve from track coordinates
    // =========================================================================
    std::cout << "=== STEP 2: Creating Parametric Curve ===\n\n";
    
    // Extract (x, y) coordinates into a matrix
    // Matrix format: N rows x 2 columns, where each row is (x, y)
    size_t numPoints = telemetry.size();
    MML::Matrix<Real> trackPoints(static_cast<int>(numPoints), 2);
    
    for (size_t i = 0; i < numPoints; i++) {
        trackPoints(static_cast<int>(i), 0) = telemetry[i].x_m;
        trackPoints(static_cast<int>(i), 1) = telemetry[i].y_m;
    }
    
    std::cout << "  Extracted " << numPoints << " track coordinates\n";
    
    // Calculate X and Y ranges manually
    Real minX = telemetry[0].x_m, maxX = telemetry[0].x_m;
    Real minY = telemetry[0].y_m, maxY = telemetry[0].y_m;
    for (size_t i = 1; i < numPoints; i++) {
        if (telemetry[i].x_m < minX) minX = telemetry[i].x_m;
        if (telemetry[i].x_m > maxX) maxX = telemetry[i].x_m;
        if (telemetry[i].y_m < minY) minY = telemetry[i].y_m;
        if (telemetry[i].y_m > maxY) maxY = telemetry[i].y_m;
    }
    std::cout << "  X range: [" << minX << ", " << maxX << "] m\n";
    std::cout << "  Y range: [" << minY << ", " << maxY << "] m\n\n";
    
    // Create spline-interpolated parametric curve
    // Parameter t goes from 0 to total distance (5793 m)
    Real minDist = telemetry[0].distance_m;
    Real maxDist = telemetry.totalDistance_m;
    
    // Use SplineInterpParametricCurve<2> for smooth track representation
    // The curve is parameterized by distance along the track
    MML::SplineInterpParametricCurve<2> trackCurve(minDist, maxDist, trackPoints, false);
    
    std::cout << "  Created SplineInterpParametricCurve<2>\n";
    std::cout << "  Parameter range: [" << minDist << ", " << maxDist << "] m\n\n";
    
    // Verify curve by sampling a few points
    std::cout << "  Curve verification (comparing original vs interpolated):\n";
    std::cout << "  Distance(m) |  Original X |  Original Y |  Interp X |  Interp Y | Error(m)\n";
    std::cout << "  ------------|-------------|-------------|-----------|-----------|--------\n";
    
    std::vector<size_t> testIndices = {0, 100, 300, 500, 700, numPoints - 1};
    Real maxError = 0.0;
    
    for (size_t idx : testIndices) {
        if (idx >= numPoints) continue;
        
        Real dist = telemetry[idx].distance_m;
        Real origX = telemetry[idx].x_m;
        Real origY = telemetry[idx].y_m;
        
        MML::VectorN<Real, 2> interpPoint = trackCurve(dist);
        Real interpX = interpPoint[0];
        Real interpY = interpPoint[1];
        
        Real error = std::sqrt((origX - interpX) * (origX - interpX) + 
                               (origY - interpY) * (origY - interpY));
        if (error > maxError) maxError = error;
        
        std::cout << "  " << std::setw(11) << dist << " | "
                  << std::setw(11) << origX << " | "
                  << std::setw(11) << origY << " | "
                  << std::setw(9) << interpX << " | "
                  << std::setw(9) << interpY << " | "
                  << std::setw(6) << error << "\n";
    }
    std::cout << "\n  Maximum interpolation error: " << maxError << " m\n\n";
    
    // Visualize the track
    std::cout << "  Visualizing track layout...\n";
    std::string vizTitle = trackTitle + " F1 Track - " + driverInfo + " Fastest Lap";
    std::string vizFile = "f1_" + trackName + "_track_layout";
    
    auto result = MML::Visualizer::VisualizeParamCurve2D(
        trackCurve, 
        vizTitle,
        minDist, maxDist, 
        1000,  // Sample 1000 points for smooth visualization
        vizFile
    );
    
    if (result.success) {
        std::cout << "  Track visualization saved!\n\n";
    } else {
        std::cout << "  Note: Visualization skipped (" << result.errorMessage << ")\n\n";
    }
    
    std::cout << "=== Step 2 Complete! Parametric curve created. ===\n\n";
    
    // =========================================================================
    // STEP 3: Create speed interpolation and calculate derivatives
    // =========================================================================
    std::cout << "=== STEP 3: Speed Interpolation & Derivatives ===\n\n";
    
    // Create vectors for distance and speed
    MML::Vector<Real> distances(static_cast<int>(numPoints));
    MML::Vector<Real> speeds_ms(static_cast<int>(numPoints));  // Speed in m/s
    MML::Vector<Real> times(static_cast<int>(numPoints));
    
    for (size_t i = 0; i < numPoints; i++) {
        distances[static_cast<int>(i)] = telemetry[i].distance_m;
        speeds_ms[static_cast<int>(i)] = telemetry[i].speed_ms();  // Convert km/h to m/s
        times[static_cast<int>(i)] = telemetry[i].time_s;
    }
    
    // Create spline interpolation for speed as function of distance: v(s)
    MML::SplineInterpRealFunc speedInterp(distances, speeds_ms);
    
    // Create km/h version for visualization
    MML::Vector<Real> speeds_kmh(static_cast<int>(numPoints));
    for (size_t i = 0; i < numPoints; i++)
        speeds_kmh[static_cast<int>(i)] = telemetry[i].speed_kmh;
    MML::SplineInterpRealFunc speedInterpKmh(distances, speeds_kmh);
    
    // Create spline interpolation for time as function of distance: t(s)
    MML::SplineInterpRealFunc timeInterp(distances, times);
    
    // Find min/max speed
    Real minSpeed = speeds_ms[0], maxSpeed = speeds_ms[0];
    for (int i = 1; i < static_cast<int>(numPoints); i++) {
        if (speeds_ms[i] < minSpeed) minSpeed = speeds_ms[i];
        if (speeds_ms[i] > maxSpeed) maxSpeed = speeds_ms[i];
    }
    
    std::cout << "  Created speed interpolation v(s)\n";
    std::cout << "  Speed range: " << minSpeed << " - " << maxSpeed << " m/s\n";
    std::cout << "              (" << minSpeed * 3.6 << " - " << maxSpeed * 3.6 << " km/h)\n\n";
    
    // =========================================================================
    // STEP 4: Calculate curvature and G-forces at each point
    // =========================================================================
    std::cout << "=== STEP 4: G-Force Calculation ===\n\n";
    
    const Real g = 9.81;  // Gravitational acceleration (m/s²)
    
    // Sample points along the track for G-force calculation
    int numSamples = 500;
    Real ds = (maxDist - minDist) / (numSamples - 1);
    
    // Storage for results
    std::vector<Real> sample_dist(numSamples);
    std::vector<Real> sample_speed(numSamples);
    std::vector<Real> sample_curvature(numSamples);
    std::vector<Real> sample_g_lateral(numSamples);
    std::vector<Real> sample_g_longitudinal(numSamples);
    std::vector<Real> sample_g_total(numSamples);
    
    Real maxGLat = 0, maxGLong = 0, maxGTotal = 0;
    Real minGLong = 0;  // For braking (negative)
    
    for (int i = 0; i < numSamples; i++) {
        Real s = minDist + i * ds;
        sample_dist[i] = s;
        
        // Get speed at this distance
        Real v = speedInterp(s);
        sample_speed[i] = v;
        
        // Calculate curvature using MML's getCurvature2D() function
        // Uses formula: κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
        Real curvature = Curves::getCurvature2D(trackCurve, s);
        sample_curvature[i] = curvature;
        
        // For longitudinal G calculation, we still need numerical derivatives for boundary handling
        Real h = 1.0;  // Step size for speed derivative (1 meter)
        Real s_minus = std::max(minDist, s - h);
        Real s_plus = std::min(maxDist, s + h);
        Real actual_h = (s_plus - s_minus) / 2.0;
        
        // Lateral G-force: G_lat = v² · κ / g
        Real g_lat = (v * v * curvature) / g;
        sample_g_lateral[i] = g_lat;
        if (g_lat > maxGLat) maxGLat = g_lat;
        
        // Longitudinal G-force: G_long = (1/g) · dv/dt = (1/g) · (dv/ds) · (ds/dt) = (1/g) · v · (dv/ds)
        // Using numerical derivative of speed
        Real v_minus = speedInterp(s_minus);
        Real v_plus = speedInterp(s_plus);
        Real dv_ds = (v_plus - v_minus) / (2.0 * actual_h);
        
        Real g_long = (v * dv_ds) / g;
        sample_g_longitudinal[i] = g_long;
        if (g_long > maxGLong) maxGLong = g_long;
        if (g_long < minGLong) minGLong = g_long;
        
        // Total G-force
        Real g_total = std::sqrt(g_lat * g_lat + g_long * g_long);
        sample_g_total[i] = g_total;
        if (g_total > maxGTotal) maxGTotal = g_total;
    }
    
    std::cout << "  Calculated G-forces at " << numSamples << " points\n\n";
    std::cout << "  G-Force Summary:\n";
    std::cout << "  ┌─────────────────────────────────────────┐\n";
    std::cout << "  │  Max Lateral G:      " << std::fixed << std::setprecision(2) 
              << std::setw(6) << maxGLat << " G          │\n";
    std::cout << "  │  Max Acceleration:   " << std::setw(6) << maxGLong << " G          │\n";
    std::cout << "  │  Max Braking:        " << std::setw(6) << -minGLong << " G          │\n";
    std::cout << "  │  Max Total G:        " << std::setw(6) << maxGTotal << " G          │\n";
    std::cout << "  └─────────────────────────────────────────┘\n\n";
    
    // Find peak G-force locations
    std::cout << "  Peak G-force locations:\n";
    for (int i = 1; i < numSamples - 1; i++) {
        // Local maxima for lateral G (corners)
        if (sample_g_lateral[i] > 3.0 && 
            sample_g_lateral[i] > sample_g_lateral[i-1] && 
            sample_g_lateral[i] > sample_g_lateral[i+1]) {
            std::cout << "    Corner at " << std::setw(6) << static_cast<int>(sample_dist[i]) 
                      << "m: " << std::setprecision(2) << sample_g_lateral[i] << " G lateral, "
                      << static_cast<int>(sample_speed[i] * 3.6) << " km/h\n";
        }
    }
    
    // =========================================================================
    // STEP 5: Visualize G-forces
    // =========================================================================
    std::cout << "\n=== STEP 5: Visualizing G-Forces ===\n\n";
    
    // Create function wrappers for visualization
    // We need to create RealFunction objects that can be passed to the visualizer
    
    // For now, let's output the data to a file for external plotting
    std::string gforceFile = "results/f1_" + trackName + "_gforces.csv";
    std::ofstream outFile(gforceFile);
    if (outFile.is_open()) {
        outFile << "distance_m,speed_kmh,curvature,g_lateral,g_longitudinal,g_total\n";
        for (int i = 0; i < numSamples; i++) {
            outFile << std::fixed << std::setprecision(4)
                    << sample_dist[i] << ","
                    << sample_speed[i] * 3.6 << ","
                    << sample_curvature[i] << ","
                    << sample_g_lateral[i] << ","
                    << sample_g_longitudinal[i] << ","
                    << sample_g_total[i] << "\n";
        }
        outFile.close();
        std::cout << "  G-force data saved to: " << gforceFile << "\n";
    }
    
    // Create spline interpolations for G-force visualization
    MML::Vector<Real> dist_vec(numSamples);
    MML::Vector<Real> g_lat_vec(numSamples);
    MML::Vector<Real> g_long_vec(numSamples);
    MML::Vector<Real> g_total_vec(numSamples);
    
    for (int i = 0; i < numSamples; i++) {
        dist_vec[i] = sample_dist[i];
        g_lat_vec[i] = sample_g_lateral[i];
        g_long_vec[i] = std::abs(sample_g_longitudinal[i]);  // Use absolute value for display
        g_total_vec[i] = sample_g_total[i];
    }
    
    // Create spline interpolations for each G-force component
    MML::SplineInterpRealFunc gLatInterp(dist_vec, g_lat_vec);
    MML::SplineInterpRealFunc gLongInterp(dist_vec, g_long_vec);
    MML::SplineInterpRealFunc gTotalInterp(dist_vec, g_total_vec);
    
    // Visualize all 3 G-forces on the same plot
    std::vector<MML::SplineInterpRealFunc> gForceFuncs = {gLatInterp, gLongInterp, gTotalInterp};
    std::vector<std::string> gForceLegend = {"Lateral G", "Long. G (abs)", "Total G"};
    
    std::string gforceVizFile = "f1_" + trackName + "_gforces_plot";
    auto gforceResult = MML::Visualizer::VisualizeMultiRealFunction(
        gForceFuncs,
        trackTitle + " - G-Forces vs Distance (" + driverInfo + ")",
        gForceLegend,
        minDist, maxDist,
        500,
        gforceVizFile
    );
    
    if (gforceResult.success) {
        std::cout << "  G-force graph (3 curves) saved!\n";
    } else {
        std::cout << "  G-force visualization: " << gforceResult.errorMessage << "\n";
    }
    
    // Create a simple visualization showing speed profile
    // Using SplineInterpRealFunc to visualize speed vs distance
    std::string speedVizFile = "f1_" + trackName + "_speed_profile";
    auto speedResult = MML::Visualizer::VisualizeRealFunction(
        speedInterpKmh,
        trackTitle + " - Speed Profile (km/h)",
        minDist, maxDist,
        500,
        speedVizFile
    );
    
    if (speedResult.success) {
        std::cout << "  Speed profile visualization saved!\n";
    }
    
    std::cout << "\n╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║              F1 G-FORCE ANALYSIS COMPLETE!                   ║\n";
    std::cout << "╠══════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Track: " << std::left << std::setw(20) << trackTitle 
              << "Driver: " << std::setw(22) << driverInfo << "║\n";
    std::cout << "║  Lap Time: " << std::setw(8) << std::fixed << std::setprecision(2) 
              << telemetry.totalTime_s << "s" << std::setw(40) << " " << "║\n";
    std::cout << "║  Max Lateral G: " << std::setw(6) << maxGLat << " G (cornering)" 
              << std::setw(24) << " " << "║\n";
    std::cout << "║  Max Braking G: " << std::setw(6) << -minGLong << " G" 
              << std::setw(32) << " " << "║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    
    return 0;
}
