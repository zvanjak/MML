#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "basic_types/CoordSystem.h"
#include "core/CoordTransf.h"
#endif

using namespace MML;

void Demo_Inertial()
{
    // zadati isti origin, neku brzinu, i vidjeti kako ide transformacija
}

void Demo_2DPolar()
{
    std::cout << "***************    OBLIQUE RECTILINEAR   ****************" << std::endl;

    // new basis vectors (keeping z coordinate out of it ;)
    Vector3Cartesian e1_base{1, 3, 0}, e2_base{4, 0, 0}, e3_base{0, 0, 1};      

    CoordTransfRectilinear transf(e1_base, e2_base, e3_base);

    std::cout << "Is right-handed? - " << transf.IsRightHanded() << std::endl;
    
    std::cout << "New basis vectors:\n";
    std::cout << "e1 = "; e1_base.Print(std::cout,10,5) << std::endl;
    std::cout << "e2 = "; e2_base.Print(std::cout,10,5) << std::endl;
    std::cout << "e3 = "; e3_base.Print(std::cout,10,5) << std::endl;

    std::cout << "Calculated dual basis vectors:\n";
    std::cout << "e1 dual = "; transf.Dual(0).Print(std::cout,10,5) << std::endl;
    std::cout << "e2 dual = "; transf.Dual(1).Print(std::cout,10,5) << std::endl;
    std::cout << "e3 dual = "; transf.Dual(2).Print(std::cout,10,5) << std::endl;

    std::cout << "\nTransf. matrix:\n" << transf.getAlpha() << std::endl;
    std::cout << "Inverse transf. matrix:\n" << transf.getTransf() << std::endl;

    Vector3Cartesian vec_A{7, 2, 0};
    std::cout << "Vector A (orig)      = "; vec_A.Print(std::cout,10,5) << std::endl;

    // std::cout << "alpha * A            = " << transf._alpha * vec_A << std::endl;
    // std::cout << "transf * (alpha * A) = " << transf._transf * (transf._alpha * vec_A) << std::endl;

    std::cout << "\nTransformed to new basis:\n";
    auto contravar_coef = transf.transf(vec_A);
    auto back_contravar = transf.transfInverse(contravar_coef);
    std::cout << "Contravar. coeff.   = "; contravar_coef.Print(std::cout,10,5) << std::endl;
    std::cout << "Back transf coeff   = "; back_contravar.Print(std::cout,10,5) << std::endl;

    std::cout << "\nCalculating contravariant and covariant components in new basis:\n";
    Vector3Cartesian x_dummy{1.0, 1.0, 1.0};          // point of transf. application is irrelevant (it is linear)
     
    auto contravar_vec  = transf.transfVecContravariant(vec_A, x_dummy);
    auto contravar_vec2 = transf.transfInverseVecContravariant(contravar_vec, x_dummy);
    std::cout << "Contravar. comp. of A = "; contravar_vec.Print(std::cout,10,5) << std::endl;
    std::cout << "Back transf.contravar = "; contravar_vec2.Print(std::cout,10,5) << std::endl;

    auto covar_vec  = transf.transfVecCovariant(vec_A, x_dummy);
    auto covar_vec2 = transf.transfInverseVecCovariant(covar_vec, x_dummy);
    std::cout << "Covariant  comp. of A = "; covar_vec.Print(std::cout,10,5) << std::endl;
    std::cout << "Back transf.covariant = "; covar_vec2.Print(std::cout,10,5) << std::endl;
}

void Demo_Coord_Rectilinear()
{
    std::cout << "\n****    RECTILINEAR COORD SYSTEM   ****" << std::endl;

//    Vector3Cartesian e1_base{1, 0, 0}, e2_base{0, 1, 0}, e3_base{0, 0.1, 1};      // ovo definira dot koord system
    Vector3Cartesian e1_base{1, 3, 0}, e2_base{4, 0, 0}, e3_base{0, 0, 1};      

    CoordTransfRectilinear transf(e1_base, e2_base, e3_base);

    std::cout << "e1 = " << e1_base << std::endl;
    std::cout << "e2 = " << e2_base << std::endl;
    std::cout << "e3 = " << e3_base << std::endl;

    std::cout << "e1 dual = " << transf.Dual(0) << std::endl;
    std::cout << "e2 dual = " << transf.Dual(1) << std::endl;
    std::cout << "e3 dual = " << transf.Dual(2) << std::endl;

    Vector3Cartesian vec_A{7.0, 2.0, 0.0};
    std::cout << "\nVector A: " << vec_A << std::endl;

    auto contravar_coef = transf.transf(vec_A);

    std::cout << "Contravar. coeff.: " << contravar_coef << std::endl;
    VectorN<Real, 3> contra_expanded = contravar_coef[0] * e1_base +
                                      contravar_coef[1] * e2_base +
                                      contravar_coef[2] * e3_base;
    std::cout << "Expanded to orig.: " << contra_expanded << std::endl;

    Vector3Cartesian covar_coef{ScalarProd(vec_A, Vector3Cartesian(e1_base)),
                                     ScalarProd(vec_A, Vector3Cartesian(e2_base)),
                                     ScalarProd(vec_A, Vector3Cartesian(e3_base))};        
    std::cout << "Covar. coeff.    : " << covar_coef << std::endl;
    
    VectorN<Real, 3> covar_expanded = covar_coef[0] * transf.Dual(0) +
                                     covar_coef[1] * transf.Dual(1) +
                                     covar_coef[2] * transf.Dual(2);
    std::cout << "Expanded to orig.: " << covar_expanded << std::endl;

}

void Demo_Coord_OrthogonalCartesian()
{
    std::cout << "\n****    ORTHOGONAL CARTESIAN COORD SYSTEM   ****" << std::endl;

    Vector3Cartesian b1{1, 0, 0}, b2{0, 1, 0}, b3{0, 0, 1};

    CoordSystemOrthogonalCartesian csys(b1, b2, b3);

    std::cout << "Is orthogonal = " << csys.isOrthogonal() << std::endl;
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

    InertialMovingFrame train(Vector3Cartesian{0,0,0}, Vector3Cartesian{1,0,0});
    InertialMovingFrame table_moving_on_train(Vector3Cartesian{1,0,0}, Vector3Cartesian{0.5,0,0});
    RotatingFrame disc_spinning_on_table;
    InertialMovingFrame ball_thrown_up;

    // point on disc - polar coordinates
    Vector3Cartesian point_on_disc{0.5, 0.5, 0};
    auto transf0 = disc_spinning_on_table.transf(point_on_disc, 1.0);
    auto transf1 = table_moving_on_train.transf(transf0, 1.0);
    auto transf2 = train.transf(transf1, 1.0);
/*
    // TODO - simulation of Solar system
    
    - sunce nam je u ishodistu i definira koord sustav
    - moving coord system earth_rotation_around_sun, jupier around sun, mars around sun
        - elipsa s centrom mase u zaristu
        - definira tocku centra planete u odnosu na ishodiste
    - rotating coord system earth_rotation_around_earth_center
        - sferni coord system
        - daje transformaciju iz lokalne (lat, long, height) coord u globalnu (x,y,z) coord
    - local observer coord system - "nebo"
        - za danu lat, long tocku na zemlji
        - definira oko te tocke sferni sustav, gdje je azimut poravnan sa smjerom sjevera kako treba
    - SUSTINA - zelimo transformaciju jupitera i marsa na taj sustav, da vidimo gdje se nalaze na nebo
*/

    RotatingFrame earthSolarSystemCM;    // opisuje kružnicu oko sunca - KLJUČNO JE PORAVNANJE!!!
    // T = 0 - Greenwhich je točno na x osi
    // MovingCoordSystem<3> earthRotational;       // kao koordinate ima sferne koordinate iznad površine zemlje
    // MovingCoordSystem<3> earthLocal;            // projekcija na 2D coord system, sa zadanim ishodištem

    // ČIM GIBANJE čestice po putanji odstupa od tangente, znači da je prisutna sila

    // najprije u 2D ?

    // imamo statičnu transf - sve zadano, pa da vidimo

    // imamo simulaciju - ovisnost o t, pa se mijenja

    // usporediti
    //  1. kompletnqa simulacija zemlje oko sunca i topovske granate ispaljene
    //  2. zemlja se zadano giva oko sunca, a imamo lokalni rotacijski sustav
}