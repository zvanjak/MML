/******************************************************************************
 * CHAPTER 9: THE N-BODY GRAVITY PROBLEM
 * ============================================================================
 * 
 * The crown jewel of classical mechanics simulation!
 * 
 * This chapter demonstrates gravitational N-body simulation - one of the oldest
 * and most fascinating problems in physics. From Newton's original work on 
 * planetary motion to modern galaxy simulations, N-body gravity is fundamental.
 * 
 * We present three spectacular demonstrations:
 * 
 *   DEMO 1: SOLAR SYSTEM
 *           Our actual solar system with real astronomical data!
 *           Sun + 8 planets with proper masses and orbital distances.
 * 
 *   DEMO 2: MANY BODIES  
 *           100 bodies orbiting a massive central object.
 *           Shows the dynamics of a miniature star cluster.
 * 
 *   DEMO 3: STAR CLUSTER COLLISION  *** THE SPECTACULAR ONE! ***
 *           TWO 100-body star clusters approaching and COLLIDING!
 *           Watch chaos emerge from order as gravitational interactions
 *           scatter stars in all directions!
 * 
 * Physics: Newton's Law of Universal Gravitation
 *          F = G * m1 * m2 / r^2  (attractive force between masses)
 * 
 * Numerical Methods:
 *   - Euler method (simple, illustrative)
 *   - RK5 Cash-Karp (adaptive, production-quality)
 *   - Energy conservation monitoring
 * 
 * Output: Trajectory files + 3D animated visualizations
 * 
 *****************************************************************************/

#include <iostream>
#include <iomanip>
#include <string>

// Demo function declarations
void Demo_FiveMasses();
void Demo_Solar_system();
void Demo_ManyBodies();
void Demo_StarClusterCollision();

void Chapter09_N_body_problem()
{
    std::cout << "\n";
    std::cout << "======================================================================\n";
    std::cout << "                                                                      \n";
    std::cout << "     ██████╗██╗  ██╗ █████╗ ██████╗ ████████╗███████╗██████╗  █████╗  \n";
    std::cout << "    ██╔════╝██║  ██║██╔══██╗██╔══██╗╚══██╔══╝██╔════╝██╔══██╗██╔══██╗ \n";
    std::cout << "    ██║     ███████║███████║██████╔╝   ██║   █████╗  ██████╔╝ ╚█████╔╝\n";
    std::cout << "    ██║     ██╔══██║██╔══██║██╔═══╝    ██║   ██╔══╝  ██╔══██╗██╔═══██╗\n";
    std::cout << "    ╚██████╗██║  ██║██║  ██║██║        ██║   ███████╗██║  ██║╚██████╔╝\n";
    std::cout << "     ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝        ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ \n";
    std::cout << "                                                                      \n";
    std::cout << "          T H E   N - B O D Y   G R A V I T Y   P R O B L E M         \n";
    std::cout << "                                                                      \n";
    std::cout << "======================================================================\n\n";

    std::cout << "  Newton's Law of Universal Gravitation: F = G * m1 * m2 / r^2\n\n";

    std::cout << "  This chapter explores gravitational interactions between multiple\n";
    std::cout << "  bodies - from our solar system to colliding star clusters!\n\n";

    std::cout << "  Demonstrations:\n";
    std::cout << "  ---------------\n";
    std::cout << "    1. Solar System    - Real astronomical data (Sun + 8 planets)\n";
    std::cout << "    2. Many Bodies     - 100 bodies orbiting central mass\n";
    std::cout << "    3. Cluster Collision - TWO 100-body clusters COLLIDING!\n\n";

    std::cout << std::string(70, '-') << "\n\n";

    // Demo 1: Simple 5-body system showing conservation laws  
    // Demo_FiveMasses();  // Uncomment for basic conservation law demo

    // Demo 2: Solar system with real planetary data
    Demo_Solar_system();

    // Demo 3: 100 bodies orbiting a central mass
    Demo_ManyBodies();

    // Demo 4: THE SPECTACULAR ONE! Two star clusters colliding!
    Demo_StarClusterCollision();

    // Chapter summary
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  CHAPTER 9 COMPLETE!\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "\n  All N-body gravity demonstrations finished successfully.\n";
    std::cout << "  Check the 'results/' folder for trajectory and animation files.\n";
    std::cout << "  Use MML Visualizer to view 3D particle simulations!\n\n";
}
