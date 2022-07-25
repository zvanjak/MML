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
    std::cout << "****    2D DUALS   ****" << std::endl;

   MML::Vector3Cartesian e1_base{1, 3, 0}, e2_base{4, 0, 0}, e3_base{0, 0, 1};      // ovo definira dot koord system

    MML::CoordTransfRectilinear transf(e1_base, e2_base, e3_base);

    std::cout << "e1 = " << e1_base << std::endl;
    std::cout << "e2 = " << e2_base << std::endl;
    std::cout << "e3 = " << e3_base << std::endl;

    // duali
    std::cout << "e1 dual = " << transf.Dual(0) << std::endl;
    std::cout << "e2 dual = " << transf.Dual(1) << std::endl;
    std::cout << "e3 dual = " << transf.Dual(2) << std::endl;

    MML::Vector3Cartesian vec_A{7, 2, 0};
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


        std::cout << "AT POINT:\n";
    MML::Vector3Cartesian x1_cart{1.0, 0.0, 0.0};
    MML::VectorN<3> x1_transf{ transf.transf(x1_cart) };
    std::cout << "Cartesian : " << x1_cart << std::endl << "Transformed : " << x1_transf << std::endl; 
     
    auto contravar_vec = transf.contravariantTransf(vec_A, x1_transf);
    std::cout << "Contra. comp. at " << x1_transf << " = " << contravar_vec << std::endl;

    auto covar_vec = transf.covariantTransf(vec_A, x1_cart);
    std::cout << "Covar. comp. at  " << x1_transf << " = " << covar_vec << std::endl;

    // auto v_back_transf_to_cart = transfCartToSpher.contravariantTransf(v_transf_to_spher, x1_cart);
    // std::cout << "v (cart) at " << x1_cart << " = " << v_back_transf_to_cart << std::endl;

}

void Demo_Covar_contravar()
{
    std::cout << "****    COVAR & CONTRAVAR   ****" << std::endl;

    MML::Vector3Cartesian e1_base{1, 0, 0}, e2_base{0, 1, 0}, e3_base{0, 0.1, 1};      // ovo definira dot koord system

    std::cout << "e1 = " << e1_base << std::endl;
    std::cout << "e2 = " << e2_base << std::endl;
    std::cout << "e3 = " << e3_base << std::endl;

    // idemo napraviti duale
    MML::Vector3Cartesian cross1 = VectorProd(e2_base, e3_base);
    MML::VectorN<3> e1_dual = (1 / (ScalarProd(e1_base, cross1))) * cross1;

    MML::Vector3Cartesian cross2 = VectorProd(e3_base, e1_base);
    MML::VectorN<3> e2_dual = (1 / (ScalarProd(e2_base, cross2))) * cross2;

    MML::Vector3Cartesian cross3 = VectorProd(e1_base, e2_base);
    MML::VectorN<3> e3_dual = (1 / (ScalarProd(e3_base, cross3))) * cross3;

    std::cout << "e1 dual = " << e1_dual << std::endl;
    std::cout << "e2 dual = " << e2_dual << std::endl;
    std::cout << "e3 dual = " << e3_dual << std::endl;

    MML::Vector3Cartesian vec_A{0.0, 1.0, 1.0};
    std::cout << "\nVector A: " << vec_A << std::endl;

    MML::Vector3Cartesian contravar_coef{ScalarProd(vec_A, MML::Vector3Cartesian(e1_dual)),
                                         ScalarProd(vec_A, MML::Vector3Cartesian(e2_dual)),
                                         ScalarProd(vec_A, MML::Vector3Cartesian(e3_dual))};
    std::cout << "Contravar. coeff.: " << contravar_coef << std::endl;
    
    MML::VectorN<3> contra_expanded = contravar_coef[0] * e1_base +
                                      contravar_coef[1] * e2_base +
                                      contravar_coef[2] * e3_base;
    std::cout << "Expanded to orig.: " << contra_expanded << std::endl;

    MML::Vector3Cartesian covar_coef{ScalarProd(vec_A, MML::Vector3Cartesian(e1_base)),
                                     ScalarProd(vec_A, MML::Vector3Cartesian(e2_base)),
                                     ScalarProd(vec_A, MML::Vector3Cartesian(e3_base))};        
    std::cout << "Covar. coeff.    : " << covar_coef << std::endl;
    
    MML::VectorN<3> covar_expanded = covar_coef[0] * e1_dual +
                                     covar_coef[1] * e2_dual +
                                     covar_coef[2] * e3_dual;
    std::cout << "Expanded to orig.: " << covar_expanded << std::endl;

}

void Demo_Coord_Rectilinear()
{
    std::cout << "****    RECTILINEAR COORD SYSTEM   ****" << std::endl;

    MML::Vector3Cartesian e1_base{1, 0, 0}, e2_base{0, 1, 0}, e3_base{0, 0.1, 1};      // ovo definira dot koord system

    MML::CoordTransfRectilinear transf(e1_base, e2_base, e3_base);

    std::cout << "e1 = " << e1_base << std::endl;
    std::cout << "e2 = " << e2_base << std::endl;
    std::cout << "e3 = " << e3_base << std::endl;

    // idemo napraviti duale
    // std::cout << "e1 dual = " << e1_dual << std::endl;
    // std::cout << "e2 dual = " << e2_dual << std::endl;
    // std::cout << "e3 dual = " << e3_dual << std::endl;

    MML::Vector3Cartesian vec_A{0.0, 1.0, 1.0};
    std::cout << "\nVector A: " << vec_A << std::endl;

    auto contravar_coef = transf.transf(vec_A);

    std::cout << "Contravar. coeff.: " << contravar_coef << std::endl;
    MML::VectorN<3> contra_expanded = contravar_coef[0] * e1_base +
                                      contravar_coef[1] * e2_base +
                                      contravar_coef[2] * e3_base;
    std::cout << "Expanded to orig.: " << contra_expanded << std::endl;

    // MML::Vector3Cartesian covar_coef{ScalarProd(vec_A, MML::Vector3Cartesian(e1_base)),
    //                                  ScalarProd(vec_A, MML::Vector3Cartesian(e2_base)),
    //                                  ScalarProd(vec_A, MML::Vector3Cartesian(e3_base))};        
    // std::cout << "Covar. coeff.    : " << covar_coef << std::endl;
    
    // MML::VectorN<3> covar_expanded = covar_coef[0] * e1_dual +
    //                                  covar_coef[1] * e2_dual +
    //                                  covar_coef[2] * e3_dual;
    // std::cout << "Expanded to orig.: " << covar_expanded << std::endl;

}

void Demo_CoordSystem()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         COORD SYSTEM                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_2DPolar();

    Demo_Covar_contravar();

    Demo_Coord_Rectilinear();

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