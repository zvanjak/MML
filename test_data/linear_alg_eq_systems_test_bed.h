#if !defined __MML_LINEAR_ALG_EQ_SOLVERS_TEST_BED_H
#define __MML_LINEAR_ALG_EQ_SOLVERS_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Matrix.h"
#endif

#include "linear_alg_eq_systems_defs.h"

namespace MML::TestBeds
{
    // TODO - define 10 linear systems, up to 20 x 20 dimension, with complete solutions
    // TODO - implement generation of random linear system
    // TODO - define 5 symmetric linear systems
    class TestLinearSystem
    {
    public:
        int _n;
        MatrixDbl _mat;
        VectorDbl _rhs;
        VectorDbl _sol;
        VectorComplex _eigen_values;
        std::vector<VectorComplex> _eigen_vectors;

        TestLinearSystem(int n, const MML::MatrixDbl &mat, const MML::VectorDbl &rhs, const MML::VectorDbl &sol, const MML::VectorComplex &eigen_values) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values)
        {}
        TestLinearSystem(int n, const MML::MatrixDbl &mat, const MML::VectorDbl &rhs, const MML::VectorDbl &sol, const MML::VectorComplex &eigen_values, const std::vector<MML::VectorComplex> &eigen_vectors) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values), _eigen_vectors(eigen_vectors)
        {}
    };

    class TestLinearSystemMultiRHS
    {
    public:
        int _n;
        MatrixDbl _mat;
        MatrixDbl _rhs;
        MatrixDbl _sol;

        TestLinearSystemMultiRHS(int n, const MML::MatrixDbl &mat, const MML::MatrixDbl &rhs, const MML::MatrixDbl &sol) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol)
        {}
    };    

    class TestLinearSystemSymmetric
    {
    public:
        int _n;
        MatrixDbl _mat;        // MatrixSymm is not implemented yet
        VectorDbl _rhs;
        VectorDbl _sol;
        VectorDbl _eigen_values;
        std::vector<VectorDbl> _eigen_vectors;

        TestLinearSystemSymmetric(int n, const MML::MatrixDbl &mat, const MML::VectorDbl &rhs, const MML::VectorDbl &sol, const MML::VectorDbl &eigen_values) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values)
        {}
        TestLinearSystemSymmetric(int n, const MML::MatrixDbl &mat, const MML::VectorDbl &rhs, const MML::VectorDbl &sol, const MML::VectorDbl &eigen_values, const std::vector<MML::VectorDbl> &eigen_vectors) : 
            _n(n), _mat(mat), _rhs(rhs), _sol(sol), _eigen_values(eigen_values), _eigen_vectors(eigen_vectors)
        {}        
    };

    // known LinAlgSys
    // - matrica, rhs, solution, i eigen
    class LinearAlgEqTestBed
    {
    public:
        static int numLinAlgEqSystems() { return 3; }
        static const TestLinearSystem& getLinAlgEqSystem(int index) { return _listLinearSystems[index].second; }
        static const TestLinearSystem& getLinAlgEqSystem(std::string sysName)
        {
            for (auto &sys : _listLinearSystems)
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("LinearAlgEqTestBed::getLinAlgEqSystem: system with name " + sysName + " not found");
        }
        static int numLinAlgEqSystemsSymmetric() { return 1; }
        static const TestLinearSystemSymmetric& getLinAlgEqSystemSymmetric(int index) { return _listLinearSystemsSym[index].second; }
        static const TestLinearSystemSymmetric& getLinAlgEqSystemSymmetric(std::string sysName)
        {
            for (auto &sys : _listLinearSystemsSym)
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("LinearAlgEqTestBed::getLinAlgEqSystemSymmetric: system with name " + sysName + " not found");
        }
        static int numLinAlgEqSystemsMultiRHS() { return 1; }
        static const TestLinearSystemMultiRHS& getLinAlgEqSystemMultiRHS(int index) { return _listLinearSystemsMultiRHS[index].second; }
        static const TestLinearSystemMultiRHS& getLinAlgEqSystemMultiRHS(std::string sysName)
        {
            for (auto &sys : _listLinearSystemsMultiRHS)
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("LinearAlgEqTestBed::getLinAlgEqSystemMultiRHS: system with name " + sysName + " not found");
        }

    private:
        const static inline std::pair<std::string, TestLinearSystem> _listLinearSystems[] = { 
            {  
                "mat_3x3", { 3, mat_3x3, mat_3x3_rhs0, mat_3x3_rhs0_sol, mat_3x3_eigen_val, mat_3x3_eigen_vecs }
            },
            {  
                "mat_5x5", { 5, mat_5x5, mat_5x5_rhs0, mat_5x5_rhs0_sol, mat_5x5_eigen_val, mat_5x5_eigen_vecs }
            },
            { 
                "mat_3x3_2",
                {
                    3, 
                    MatrixDbl{3, 3, { 1.0, 2.0, -1.0, 
                                     -1.0, 5.0, 6.0, 
                                      3.0, 1.0, 1.0 }},
                    VectorDbl{1.0, 2.0, 1.0},
                    VectorDbl{0.18867924528301885, 0.41509433962264153, 0.018867924528301921},
                    VectorComplex{Complex(1.0,0.0), Complex(2.0,0.0), Complex(1.0,0.0)},
                    {
                        VectorComplex{Complex(1.0,0.0), Complex(2.0,0.0), Complex(1.0,0.0)},
                        VectorComplex{Complex(1.0,0.0), Complex(2.0,0.0), Complex(1.0,0.0)},
                        VectorComplex{Complex(1.0,0.0), Complex(2.0,0.0), Complex(1.0,0.0)}
                    }
                }
            }
        };

        const static inline std::pair<std::string, TestLinearSystemSymmetric> _listLinearSystemsSym[] = { 
            { 
                "mat_sym_4x4",
                {
                    5, 
                    MatrixDbl{5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                    2.1, 1.5, 1.1, 0.7, 5.0,
                                    2.1, 1.1, 9.6, 5.4, 8.8,
                                    7.4, 0.7, 5.4, 0.4, 8.0,
                                    9.6, 5.0, 8.8, 8.0, 7.7}},
                    VectorDbl{1.0, 2.0, 1.0, 0.0, 0.0},
                    VectorDbl{1.0, 0.0, 0.0, 0.0, 0.0},
                    VectorDbl{1.0, 0.0, 2.0, 0.0, 0.0},
                    {
                        VectorDbl{1.0, 0.0, 1.0, 0.0, 0.0},
                        VectorDbl{1.0, 0.0, 1.0, 0.0, 0.0},
                        VectorDbl{1.0, 0.0, 1.0, 0.0, 0.0},
                        VectorDbl{1.0, 0.0, 1.0, 0.0, 0.0},
                        VectorDbl{1.0, 0.0, 1.0, 0.0, 0.0}
                    }
                }
            }
        };        

        const static inline std::pair<std::string, TestLinearSystemMultiRHS> _listLinearSystemsMultiRHS[] = { 
            {
                "mat_5x5_multi_rhs1",
                {
                    5, 
                    MatrixDbl{5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                    1.6, 1.5, 1.1, 0.7, 5.0,
                                    3.8, 8.0, 9.6, 5.4, 8.8,
                                    4.6, 8.2, 8.4, 0.4, 8.0,
                                    2.6, 2.9, 0.1, 9.6, 7.7}},
                    MatrixDbl{5, 2, {1.1, 1.6, 
                                    4.7, 9.1, 
                                    0.1, 4.0, 
                                    9.3, 8.4, 
                                    0.4, 4.1}},
                    MatrixDbl{5, 2, {-3.9032710424808688,  15.643114174796667, 
                                    5.2353433160849479, -11.587503332831671, 
                                    -3.2920957702478550,   4.4111268480786325, 
                                    -1.7183300108528281,   0.21432757972725644, 
                                    1.5832710097423177,  -0.70999930382454612}}
                }
            }
        };        
    };
}

#endif