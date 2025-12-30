// book_main.cpp - Unified entry point for all book chapter examples
//
// This file calls the coordinator function from each chapter.
// Each chapter can also be built and run independently via the
// chapter-specific executables in src/book/chapters/.
//
// Build target: MML_BookApp

// Chapter function declarations
// Note: Some chapters return int (Ch01, Ch03), others void

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

// Chapter 10: Coordinate Transformations
void Chapter10_Coordinate_transformations();

// Chapter 11: Carousel
void Chapter11_Carousel();

// Chapter 12: Projectile Launch
void Chapter12_projectile_launch();

// Chapter 13: Rigid Body
void Chapter13_Rigid_body_all();

// Chapter 14: Rotations and Quaternions
void Chapter14_rotations_quaternions();

// Chapter 15: Lorentz Transformations
void Chapter15_Lorentz_transformation();

// Chapter 17: Static Electric Fields
void Chapter17_electric_charge_distribution();

// Chapter 18: Static Magnetic Fields
void Chapter18_Infinite_line_magnetic_field();

// Chapter 19: Dynamic EM Fields
void Chapter19_EM_field_investigations();

// Chapter 20: Differential Geometry
void Chapter20_Diff_geometry_curves_surfaces();

// Chapter 21: General Relativity
void Chapter21_General_relativity();


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
    Chapter10_Coordinate_transformations();
    Chapter11_Carousel();
    Chapter12_projectile_launch();
    Chapter13_Rigid_body_all();
    Chapter14_rotations_quaternions();
    
    // Part V: Relativity
    Chapter15_Lorentz_transformation();
    
    // Part VI: Electromagnetism
    Chapter17_electric_charge_distribution();
    Chapter18_Infinite_line_magnetic_field();
    Chapter19_EM_field_investigations();
    
    // Part VII: Advanced Topics
    Chapter20_Diff_geometry_curves_surfaces();
    Chapter21_General_relativity();

    return 0;
}