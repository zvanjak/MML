#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/CoordSystem.h"
#include "basic_types/CoordTransf.h"
#endif

void Demo_Inertial()
{
    // zadati isti origin, neku brzinu, i vidjeti kako ide transformacija
}

void Demo_2DPolar()
{
    std::cout << "***************    OBLIQUE RECTILINEAR   ****************" << std::endl;

    // new basis vectors (keeping z coordinate out of it ;)
    MML::Vector3Cartesian e1_base{1, 3, 0}, e2_base{4, 0, 0}, e3_base{0, 0, 1};      

    MML::CoordTransfRectilinear transf(e1_base, e2_base, e3_base);

    std::cout << "Is right-handed? - " << transf.IsRightHanded() << std::endl;
    
    std::cout << "New basis vectors:\n";
    std::cout << "e1 = "; e1_base.Print(std::cout,10,5) << std::endl;
    std::cout << "e2 = "; e2_base.Print(std::cout,10,5) << std::endl;
    std::cout << "e3 = "; e3_base.Print(std::cout,10,5) << std::endl;

    std::cout << "Calculated dual basis vectors:\n";
    std::cout << "e1 dual = "; transf.Dual(0).Print(std::cout,10,5) << std::endl;
    std::cout << "e2 dual = "; transf.Dual(1).Print(std::cout,10,5) << std::endl;
    std::cout << "e3 dual = "; transf.Dual(2).Print(std::cout,10,5) << std::endl;

    std::cout << "\nTransf. matrix:\n" << transf._inv << std::endl;
    std::cout << "Inverse transf. matrix:\n" << transf._transf << std::endl;

    MML::Vector3Cartesian vec_A{7, 2, 0};
    std::cout << "Vector A (orig)     = "; vec_A.Print(std::cout,10,5) << std::endl;

    std::cout << transf._inv * vec_A << std::endl;
    std::cout << transf._transf * (transf._inv * vec_A) << std::endl;

    std::cout << "\nTransformed to new basis:\n";
    auto contravar_coef = transf.transf(vec_A);
    std::cout << "Contravar. coeff.   = "; contravar_coef.Print(std::cout,10,5) << std::endl;
    auto back_contravar = transf.transfInverse(contravar_coef);
    std::cout << "Back transf coeff   = "; back_contravar.Print(std::cout,10,5) << std::endl;

    std::cout << "\nCalculating contravariant and covariant components in new basis:\n";
    MML::Vector3Cartesian x_dummy{1.0, 1.0, 1.0};          // point of transf. application is irrelevant (it is linear)
     
    auto contravar_vec = transf.contravariantTransf(vec_A, x_dummy);
    std::cout << "Contravar. comp.    = "; contravar_vec.Print(std::cout,10,5) << std::endl;
    auto contravar_vec2 = transf.contravariantTransfInverse(contravar_vec, x_dummy);
    std::cout << "Back transf. contra = "; contravar_vec2.Print(std::cout,10,5) << std::endl;

    auto covar_vec = transf.covariantTransf(vec_A, x_dummy);
    std::cout << "Covariant  comp.    = "; covar_vec.Print(std::cout,10,5) << std::endl;
    auto covar_vec2 = transf.covariantTransfInverse(covar_vec, x_dummy);
    std::cout << "Back transf covar   = "; covar_vec2.Print(std::cout,10,5) << std::endl;
}

void Demo_Coord_Rectilinear()
{
    std::cout << "\n****    RECTILINEAR COORD SYSTEM   ****" << std::endl;

    MML::Vector3Cartesian e1_base{1, 0, 0}, e2_base{0, 1, 0}, e3_base{0, 0.1, 1};      // ovo definira dot koord system

    MML::CoordTransfRectilinear transf(e1_base, e2_base, e3_base);

    std::cout << "e1 = " << e1_base << std::endl;
    std::cout << "e2 = " << e2_base << std::endl;
    std::cout << "e3 = " << e3_base << std::endl;

    std::cout << "e1 dual = " << transf.Dual(0) << std::endl;
    std::cout << "e2 dual = " << transf.Dual(1) << std::endl;
    std::cout << "e3 dual = " << transf.Dual(2) << std::endl;

    MML::Vector3Cartesian vec_A{0.0, 1.0, 1.0};
    std::cout << "\nVector A: " << vec_A << std::endl;

    auto contravar_coef = transf.transf(vec_A);

    std::cout << "Contravar. coeff.: " << contravar_coef << std::endl;
    MML::VectorN<3> contra_expanded = contravar_coef[0] * e1_base +
                                      contravar_coef[1] * e2_base +
                                      contravar_coef[2] * e3_base;
    std::cout << "Expanded to orig.: " << contra_expanded << std::endl;

    MML::Vector3Cartesian covar_coef{ScalarProd(vec_A, MML::Vector3Cartesian(e1_base)),
                                     ScalarProd(vec_A, MML::Vector3Cartesian(e2_base)),
                                     ScalarProd(vec_A, MML::Vector3Cartesian(e3_base))};        
    std::cout << "Covar. coeff.    : " << covar_coef << std::endl;
    
    MML::VectorN<3> covar_expanded = covar_coef[0] * transf.Dual(0) +
                                     covar_coef[1] * transf.Dual(1) +
                                     covar_coef[2] * transf.Dual(2);
    std::cout << "Expanded to orig.: " << covar_expanded << std::endl;

}

void Demo_Coord_OrthogonalCartesian()
{
    std::cout << "\n****    ORTHOGONAL CARTESIAN COORD SYSTEM   ****" << std::endl;

    MML::Vector3Cartesian b1{1, 0, 0}, b2{0, 1, 0}, b3{0, 0, 1};

    MML::CoordSystemOrthogonalCartesian csys(b1, b2, b3);

    std::cout << "Is orthogonal = " << csys.isOrthogonal() << std::endl;
}

void Demo_Coord_ObliqueCartesian()
{
    std::cout << "\n****    OBLIQUE CARTESIAN COORD SYSTEM   ****" << std::endl;

    MML::Vector3Cartesian b1{1, 0, 0}, b2{0, 1, 0}, b3{0, 0, 1};

    MML::Vector3Cartesian  e1_base{1, 3, 0}, e2_base{4, 0, 0}, e3_base{0, 0, 1};

    MML::CoordSystemObliqueCartesian csys(e1_base, e2_base, e3_base);

    std::cout << "e1 dual = " << csys._dual[0] << std::endl;
    std::cout << "e2 dual = " << csys._dual[1] << std::endl;
    std::cout << "e3 dual = " << csys._dual[2] << std::endl;

    std::cout << "Alpha mat  = " << csys._alpha << std::endl;
    std::cout << "Alpha T    = " << csys._transf << std::endl;
    std::cout << "Inv transf = " << csys._inv << std::endl;
}

void Demo_CoordSystem()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         COORD SYSTEM                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_2DPolar();

    Demo_Coord_Rectilinear();

    Demo_Coord_OrthogonalCartesian();

    Demo_Coord_ObliqueCartesian();



    MML::RotatingMovingFrame earthLocalSpherical;    // opisuje kružnicu oko sunca - KLJUČNO JE PORAVNANJE!!!
    // T = 0 - Greenwhich je točno na x osi
    MML::MovingCoordSystem<3> earthRotational;       // kao koordinate ima sferne koordinate iznad površine zemlje
    MML::MovingCoordSystem<3> earthLocal;            // projekcija na 2D coord system, sa zadanim ishodištem

    // ČIM GIBANJE čestice po putanji odstupa od tangente, znači da je prisutna sila

    // najprije u 2D ?

    // imamo statičnu transf - sve zadano, pa da vidimo

    // imamo simulaciju - ovisnost o t, pa se mijenja

    // usporediti
    //  1. kompletnqa simulacija zemlje oko sunca i topovske granate ispaljene
    //  2. zemlja se zadano giva oko sunca, a imamo lokalni rotacijski sustav
}