///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLExceptions.h                                                     ///
///  Description: Custom exception classes for error handling across MML              ///
///               Vector, Matrix, Algorithm, and Numerical errors                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_EXCEPTIONS_H
#define MML_EXCEPTIONS_H

#include <stdexcept>
#include <string>

namespace MML
{
	//////////             Vector error exceptions            ///////////
	class VectorInitializationError : public std::invalid_argument
	{
		int _size1;
	public:
		VectorInitializationError(std::string inMessage, int size1) 
			: std::invalid_argument(inMessage), _size1(size1) { }
		
		int size() const noexcept { return _size1; }
	};
	
	class VectorDimensionError : public std::invalid_argument
	{
		int _size1, _size2;
	public:
		VectorDimensionError(std::string inMessage, int size1, int size2) 
			: std::invalid_argument(inMessage), _size1(size1), _size2(size2) { }
		
		int expected() const noexcept { return _size1; }
		int actual() const noexcept { return _size2; }
	};
	
	class VectorAccessBoundsError : public std::out_of_range
	{
		int _i, _n;
	public:
		VectorAccessBoundsError(std::string inMessage, int i, int n) 
			: std::out_of_range(inMessage), _i(i), _n(n) { }
		
		int index() const noexcept { return _i; }
		int size() const noexcept { return _n; }
	};

	//////////             Matrix error exceptions            ///////////
	class MatrixAllocationError : public std::out_of_range
	{
		int _rows, _cols;
	public:
		MatrixAllocationError(std::string inMessage, int rows, int cols) 
			: std::out_of_range(inMessage), _rows(rows), _cols(cols) { }
		
		int rows() const noexcept { return _rows; }
		int cols() const noexcept { return _cols; }
	};
	class MatrixAccessBoundsError : public std::out_of_range
	{
		int _i, _j, _rows, _cols;
	public:
		MatrixAccessBoundsError(std::string inMessage, int i, int j, int rows, int cols) 
			: std::out_of_range(inMessage), _i(i), _j(j), _rows(rows), _cols(cols) { }
		
		int row() const noexcept { return _i; }
		int col() const noexcept { return _j; }
		int rows() const noexcept { return _rows; }
		int cols() const noexcept { return _cols; }
	};
	class MatrixDimensionError : public std::invalid_argument
	{
		int _rows1, _cols1, _rows2, _cols2;
	public:
		// Constructor with just a message (for simple error cases)
		explicit MatrixDimensionError(std::string inMessage) 
			: std::invalid_argument(inMessage), _rows1(-1), _cols1(-1), _rows2(-1), _cols2(-1) { }
		
		// Constructor with full dimension info (for dimension mismatch errors)
		MatrixDimensionError(std::string inMessage, int r1, int c1, int r2, int c2) 
			: std::invalid_argument(inMessage), _rows1(r1), _cols1(c1), _rows2(r2), _cols2(c2) { }
		
		int expected_rows() const noexcept { return _rows1; }
		int expected_cols() const noexcept { return _cols1; }
		int actual_rows() const noexcept { return _rows2; }
		int actual_cols() const noexcept { return _cols2; }
	};
	class SingularMatrixError : public std::domain_error
	{
		double _determinant;
		int _pivot_row;
	public:
		SingularMatrixError(std::string inMessage, double det = 0.0, int pivot_row = -1) 
			: std::domain_error(inMessage), _determinant(det), _pivot_row(pivot_row) { }
		
		double determinant() const noexcept { return _determinant; }
		int pivot_row() const noexcept { return _pivot_row; }
	};

	//////////             Integration exceptions            ///////////
	class IntegrationTooManySteps : public std::domain_error
	{
		int _steps;
		double _achieved_precision;
		double _required_precision;
	public:
		IntegrationTooManySteps(std::string inMessage, int steps = 0, 
		                       double achieved = 0.0, double required = 0.0) 
			: std::domain_error(inMessage), _steps(steps), 
			  _achieved_precision(achieved), _required_precision(required) { }
		
		int steps() const noexcept { return _steps; }
		double achieved_precision() const noexcept { return _achieved_precision; }
		double required_precision() const noexcept { return _required_precision; }
	};

	//////////           Interpolation exceptions            ///////////
	class RealFuncInterpInitError: public std::domain_error
	{
	public:
		RealFuncInterpInitError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	class RealFuncInterpRuntimeError: public std::runtime_error
	{
	public:
		RealFuncInterpRuntimeError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	////////////             Tensor exceptions             /////////////
	class TensorCovarContravarNumError : public std::invalid_argument
	{
		int _numContra, _numCo;
	public:
		TensorCovarContravarNumError(std::string inMessage, int size1, int size2) 
			: std::invalid_argument(inMessage), _numContra(size1), _numCo(size2) { }
		
		int num_contravariant() const noexcept { return _numContra; }
		int num_covariant() const noexcept { return _numCo; }
	};
	class TensorCovarContravarArithmeticError : public std::invalid_argument
	{
		int _numContra, _numCo;
		int _bContra, _bCo;
	public:
		TensorCovarContravarArithmeticError(std::string inMessage, int contra, int co, int b_contra, int b_co) 
			: std::invalid_argument(inMessage), _numContra(contra), _numCo(co), _bContra(b_contra), _bCo(b_co) { }
		
		int num_contravariant() const noexcept { return _numContra; }
		int num_covariant() const noexcept { return _numCo; }
		int other_contravariant() const noexcept { return _bContra; }
		int other_covariant() const noexcept { return _bCo; }
	};
	class TensorIndexError : public std::invalid_argument
	{
	public:
		TensorIndexError(std::string inMessage) : std::invalid_argument(inMessage)
		{ }
	};    

	//////////           Root finding exceptions            ///////////
	class RootFindingError: public std::runtime_error
	{
	public:
		RootFindingError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////            ODE solver exceptions             ///////////
	class ODESolverError: public std::runtime_error
	{
	public:
		ODESolverError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////            Statistics exceptions             ///////////
	class StatisticsError : public std::runtime_error
	{
	public:
		StatisticsError(std::string inMessage) : std::runtime_error(inMessage)
		{
		}
	};

	//////////      Numerical method exceptions             ///////////
	class NumericalMethodError : public std::runtime_error
	{
	public:
		NumericalMethodError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};
	
	class MatrixNumericalError : public std::runtime_error
	{
	public:
		MatrixNumericalError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////          Convergence exceptions              ///////////
	class ConvergenceError : public std::runtime_error
	{
		int _iterations;
		double _residual;
	public:
		ConvergenceError(std::string inMessage, int iterations = 0, double residual = 0.0) 
			: std::runtime_error(inMessage), _iterations(iterations), _residual(residual) { }
		
		int iterations() const noexcept { return _iterations; }
		double residual() const noexcept { return _residual; }
	};

	//////////            Geometry exceptions               ///////////
	class GeometryError : public std::domain_error
	{
	public:
		GeometryError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	//////////           Quaternion exceptions              ///////////
	class QuaternionError : public std::domain_error
	{
	public:
		QuaternionError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	//////////             File I/O exceptions              ///////////
	class FileIOError : public std::runtime_error
	{
		std::string _filename;
	public:
		FileIOError(std::string inMessage, std::string filename = "") 
			: std::runtime_error(inMessage), _filename(filename) { }
		
		const std::string& filename() const noexcept { return _filename; }
	};

	//////////            Fourier exceptions                ///////////
	class FourierError : public std::invalid_argument
	{
	public:
		FourierError(std::string inMessage) : std::invalid_argument(inMessage)
		{ }
	};

	//////////          Visualization exceptions            ///////////
	class VisualizerError : public std::runtime_error
	{
	public:
		VisualizerError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////        Curve Fitting exceptions              ///////////
	class CurveFittingError : public std::invalid_argument
	{
	public:
		CurveFittingError(std::string inMessage) : std::invalid_argument(inMessage)
		{ }
	};

	//////////             Data I/O exceptions              ///////////
	class DataError : public std::runtime_error
	{
	public:
		DataError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

}

#endif // MML_EXCEPTIONS_H