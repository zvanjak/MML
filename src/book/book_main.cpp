// book_main.cpp - Unified entry point for all book chapter examples
//
// This file calls the coordinator function from each chapter.
// Each chapter can also be built and run independently via the
// chapter-specific executables in src/book/chapters/.
//
// Build target: MML_BookApp

// Chapter function declarations

// Chapter 01: Basic Objects
int Chapter01_basic_objects();

// Chapter 02: Visualization
void Chapter02_Visualization();

// Chapter 03: Basic Algorithms
int Chapter_03_basic_algorithms();

// Chapter 04: Collision Simulator
void Chapter04_Collision_simulator();

// Chapter 05: Throwing Things in the Air
void Chapter5_Throwing_things_in_the_air();

// Chapter 06: Pendulum
void Chapter06_Pendulum();

// Chapter 07: Double and Spherical Pendulum
void Chapter07_DoublePendulum();

// Chapter 08: Gravity
void Chapter08_Gravity();

// Chapter 09: N-Body Problem
void Chapter09_N_body_problem();

// Chapter 11: Coordinate Transformations
void Chapter11_Coordinate_transformations();

// Chapter 12: Carousel
void Chapter12_Carousel();

// Chapter 13: Projectile Launch
void Chapter13_projectile_launch();

// Chapter 14: Rotations and Quaternions
void Chapter14_rotations_quaternions();

// Chapter 15: Rigid Body
void Chapter15_Rigid_body_all();

// Chapter 17: Lorentz Transformations
void Chapter17_Lorentz_transformation();

// Chapter 18: Twin Paradox
void Chapter18_Twin_paradox();

// Chapter 19: Static EM Fields
void Chapter19_electric_charge_distribution();

// Chapter 20: Dynamic EM Fields
void Chapter20_EM_field_investigations();

// Chapter 21: Differential Geometry
void Chapter21_Diff_geometry_curves_surfaces();

// Chapter 22: General Relativity
void Chapter22_General_relativity();


int main()
{
    // Part I: Fundamentals
    Chapter01_basic_objects();
    Chapter02_Visualization();
    Chapter_03_basic_algorithms();
    
    // Part II: Basic Mechanics
    Chapter04_Collision_simulator();
    Chapter5_Throwing_things_in_the_air();
    Chapter06_Pendulum();
    Chapter07_DoublePendulum();
    
    // Part III: Gravitation
    Chapter08_Gravity();
    Chapter09_N_body_problem();
    
    // Part IV: Advanced Mechanics
    Chapter11_Coordinate_transformations();
    Chapter12_Carousel();
    Chapter13_projectile_launch();
    Chapter14_rotations_quaternions();
    Chapter15_Rigid_body_all();
    
    // Part V: Relativity
    Chapter17_Lorentz_transformation();
    Chapter18_Twin_paradox();
    
    // Part VI: Electromagnetism
    Chapter19_electric_charge_distribution();
    Chapter20_EM_field_investigations();
    
    // Part VII: Advanced Topics
    Chapter21_Diff_geometry_curves_surfaces();
    Chapter22_General_relativity();

    return 0;
}