///////////////////////////////////////////////////////////////////////////////////////////
/// @file RigidBodySerializer.h
/// @brief Serialization utilities for rigid body simulation data
/// @details Exports rigid body trajectory data in MML format for visualization.
///
/// FORMAT: RIGID_BODY_SIMULATION_3D
/// - Container dimensions and metadata
/// - Body definitions (dimensions, mass, color)  
/// - Per-step position and orientation data
///
/// COORDINATE SYSTEM:
/// - MML uses right-handed: X right, Y up, Z toward viewer (matches WPF 3D)
/// - Quaternion format: w, x, y, z (scalar-first, matching WPF convention)
///
/// @author MinimalMathLibrary
/// @date January 2026
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_RIGID_BODY_SERIALIZER_H
#define MML_RIGID_BODY_SERIALIZER_H

#include "RigidBodySimulator.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

namespace MPL
{
    using namespace MML;
    
    /// @brief Rigid body metadata for serialization
    struct RigidBodyInfo
    {
        std::string type;               ///< Body type: "BOX" or "SPHERE"
        Real mass;                      ///< Mass in kg
        Real halfA, halfB, halfC;       ///< Half-extents (box) or radius (sphere uses halfA)
        std::string color;              ///< Color name for visualization
        std::string name;               ///< Optional body name/label
        
        RigidBodyInfo() : type("BOX"), mass(1.0), halfA(0.5), halfB(0.5), halfC(0.5), 
                          color("Gray"), name("Body") {}
        
        /// @brief Constructor for box body
        RigidBodyInfo(Real m, Real ha, Real hb, Real hc, 
                      const std::string& col = "Gray", 
                      const std::string& nm = "Body")
            : type("BOX"), mass(m), halfA(ha), halfB(hb), halfC(hc), color(col), name(nm) {}
        
        /// @brief Create a box body info
        static RigidBodyInfo Box(Real m, Real ha, Real hb, Real hc,
                                 const std::string& col = "Gray",
                                 const std::string& nm = "Box")
        {
            RigidBodyInfo info;
            info.type = "BOX";
            info.mass = m;
            info.halfA = ha;
            info.halfB = hb;
            info.halfC = hc;
            info.color = col;
            info.name = nm;
            return info;
        }
        
        /// @brief Create a sphere body info
        static RigidBodyInfo Sphere(Real m, Real radius,
                                    const std::string& col = "Gray",
                                    const std::string& nm = "Sphere")
        {
            RigidBodyInfo info;
            info.type = "SPHERE";
            info.mass = m;
            info.halfA = radius;  // Store radius in halfA
            info.halfB = radius;  // For compatibility
            info.halfC = radius;
            info.color = col;
            info.name = nm;
            return info;
        }
    };

    /// @brief Serialization result with error handling
    struct SerializeResult
    {
        bool success;
        std::string message;
        
        SerializeResult(bool ok, const std::string& msg = "") 
            : success(ok), message(msg) {}
        
        operator bool() const { return success; }
    };

    /// @brief Serializer for rigid body simulation trajectory data
    /// 
    /// Exports simulation history to MML visualization format.
    /// Supports single-body or multi-body output.
    ///
    /// @example
    /// @code
    /// RigidBodySerializer serializer;
    /// serializer.SetContainerSize(5.0);  // 10m cube
    /// serializer.AddBody(RigidBodyInfo(10.0, 1.0, 0.5, 0.3, "Red", "Box1"));
    /// serializer.AddBody(RigidBodyInfo(8.0, 0.75, 0.5, 0.4, "Blue", "Box2"));
    /// serializer.SaveSimulation("trajectory.mml", simulator.GetHistory(), 0.001);
    /// @endcode
    class RigidBodySerializer
    {
    public:
        RigidBodySerializer() : _containerHalfSize(5.0) {}
        
        /// @brief Set cubic container half-size (walls at ±halfSize)
        void SetContainerSize(Real halfSize) { _containerHalfSize = halfSize; }
        
        /// @brief Add body definition (call for each body in order)
        void AddBody(const RigidBodyInfo& info) { _bodies.push_back(info); }
        
        /// @brief Clear all body definitions
        void ClearBodies() { _bodies.clear(); }
        
        /// @brief Save complete simulation to single file
        /// @param filename Output file path (.mml extension recommended)
        /// @param history Simulation snapshot history from RigidBodySimulator
        /// @param dT Time step between frames
        /// @param saveEveryNSteps Save every Nth step to reduce file size (default: 1)
        /// @return SerializeResult with success/failure info
        SerializeResult SaveSimulation(const std::string& filename,
                                       const std::vector<SimulationSnapshot>& history,
                                       Real dT,
                                       int saveEveryNSteps = 1) const
        {
            if (filename.empty())
                return SerializeResult(false, "Filename cannot be empty");
            if (history.empty())
                return SerializeResult(false, "No simulation data to save");
            if (_bodies.empty())
                return SerializeResult(false, "No bodies defined - call AddBody() first");
            if (dT <= 0)
                return SerializeResult(false, "Time step must be positive");
            if (saveEveryNSteps < 1)
                return SerializeResult(false, "saveEveryNSteps must be >= 1");
            
            std::ofstream file(filename);
            if (!file.is_open())
                return SerializeResult(false, "Cannot open file: " + filename);
            
            try
            {
                WriteHeader(file);
                WriteBodies(file);
                WriteFrames(file, history, dT, saveEveryNSteps);
                file.close();
                
                return SerializeResult(true, "Saved " + std::to_string(CountFrames(history, saveEveryNSteps)) 
                                           + " frames to " + filename);
            }
            catch (const std::exception& e)
            {
                return SerializeResult(false, std::string("Write error: ") + e.what());
            }
        }
        
        /// @brief Save single body trajectory to separate file
        /// @param filename Output file path
        /// @param history Simulation snapshot history
        /// @param bodyIndex Which body to save (0-based)
        /// @param dT Time step
        /// @param saveEveryNSteps Frame skip factor
        SerializeResult SaveSingleBody(const std::string& filename,
                                       const std::vector<SimulationSnapshot>& history,
                                       size_t bodyIndex,
                                       Real dT,
                                       int saveEveryNSteps = 1) const
        {
            if (bodyIndex >= _bodies.size())
                return SerializeResult(false, "Body index out of range");
            
            std::ofstream file(filename);
            if (!file.is_open())
                return SerializeResult(false, "Cannot open file: " + filename);
            
            try
            {
                // Single-body format
                file << "RIGID_BODY_TRAJECTORY_3D\n";
                file << "# Single body trajectory data\n";
                file << "# Coordinate system: X-right, Y-up, Z-toward (WPF compatible)\n";
                file << "# Quaternion: w, x, y, z (scalar-first)\n\n";
                
                const auto& body = _bodies[bodyIndex];
                file << "Type: " << body.type << "\n";
                file << "Name: " << body.name << "\n";
                file << "Mass: " << body.mass << "\n";
                if (body.type == "SPHERE")
                    file << "Radius: " << body.halfA << "\n";
                else
                    file << "HalfExtents: " << body.halfA << " " << body.halfB << " " << body.halfC << "\n";
                file << "Color: " << body.color << "\n\n";
                
                file << "Container: " << _containerHalfSize * 2 << "\n\n";
                
                int numFrames = CountFrames(history, saveEveryNSteps);
                file << "NumFrames: " << numFrames << "\n";
                file << "# time, x, y, z, qw, qx, qy, qz, vx, vy, vz, wx, wy, wz\n";
                
                file << std::fixed << std::setprecision(6);
                int frameIdx = 0;
                for (size_t i = 0; i < history.size(); i += saveEveryNSteps)
                {
                    const auto& snap = history[i];
                    if (bodyIndex < snap.bodyStates.size())
                    {
                        const auto& state = snap.bodyStates[bodyIndex];
                        file << snap.time << ", "
                             << state.position.X() << ", " 
                             << state.position.Y() << ", " 
                             << state.position.Z() << ", "
                             << state.orientation.w() << ", "
                             << state.orientation.x() << ", "
                             << state.orientation.y() << ", "
                             << state.orientation.z() << ", "
                             << state.velocity.X() << ", "
                             << state.velocity.Y() << ", "
                             << state.velocity.Z() << ", "
                             << state.angularVel.X() << ", "
                             << state.angularVel.Y() << ", "
                             << state.angularVel.Z() << "\n";
                    }
                    frameIdx++;
                }
                
                file.close();
                return SerializeResult(true, "Saved " + std::to_string(numFrames) + " frames");
            }
            catch (const std::exception& e)
            {
                return SerializeResult(false, std::string("Write error: ") + e.what());
            }
        }
        
    private:
        Real _containerHalfSize;
        std::vector<RigidBodyInfo> _bodies;
        
        void WriteHeader(std::ofstream& file) const
        {
            file << "RIGID_BODY_SIMULATION_3D\n";
            file << "# Rigid body simulation trajectory data\n";
            file << "# Coordinate system: X-right, Y-up, Z-toward (WPF compatible)\n";
            file << "# Quaternion format: w, x, y, z (scalar-first)\n\n";
            
            Real size = _containerHalfSize * 2;
            file << "Container: " << size << " " << size << " " << size << "\n";
            file << "NumBodies: " << _bodies.size() << "\n\n";
        }
        
        void WriteBodies(std::ofstream& file) const
        {
            for (size_t i = 0; i < _bodies.size(); i++)
            {
                const auto& b = _bodies[i];
                file << "Body_" << (i + 1) << " " << b.type << " " << b.name << " " << b.color << " "
                     << b.mass << " " << b.halfA << " " << b.halfB << " " << b.halfC << "\n";
            }
            file << "\n";
        }
        
        void WriteFrames(std::ofstream& file, 
                        const std::vector<SimulationSnapshot>& history,
                        Real dT,
                        int saveEveryNSteps) const
        {
            int numFrames = CountFrames(history, saveEveryNSteps);
            file << "NumFrames: " << numFrames << "\n\n";
            
            file << std::fixed << std::setprecision(6);
            
            int frameIdx = 0;
            for (size_t i = 0; i < history.size(); i += saveEveryNSteps)
            {
                const auto& snap = history[i];
                file << "Frame " << frameIdx << " " << snap.time << "\n";
                
                for (size_t j = 0; j < snap.bodyStates.size() && j < _bodies.size(); j++)
                {
                    const auto& state = snap.bodyStates[j];
                    // Format: bodyIndex x y z qw qx qy qz vx vy vz wx wy wz
                    file << j << " "
                         << state.position.X() << " " 
                         << state.position.Y() << " " 
                         << state.position.Z() << " "
                         << state.orientation.w() << " "
                         << state.orientation.x() << " "
                         << state.orientation.y() << " "
                         << state.orientation.z() << " "
                         << state.velocity.X() << " "
                         << state.velocity.Y() << " "
                         << state.velocity.Z() << " "
                         << state.angularVel.X() << " "
                         << state.angularVel.Y() << " "
                         << state.angularVel.Z() << "\n";
                }
                frameIdx++;
            }
        }
        
        int CountFrames(const std::vector<SimulationSnapshot>& history, int step) const
        {
            return static_cast<int>((history.size() + step - 1) / step);
        }
    };

} // namespace MPL

#endif // MML_RIGID_BODY_SERIALIZER_H
