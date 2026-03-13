///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLExceptions.h                                                     ///
///  Description: Custom exception classes for error handling across MML              ///
///               Vector, Matrix, Algorithm, and Numerical errors                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_EXCEPTIONS_H
#define MML_EXCEPTIONS_H

#include <stdexcept>
#include <string>

namespace MML
{
	/// @brief Base marker class for all MML library exceptions.
	/// All MML exceptions inherit from both this class and a standard exception class,
	/// preserving the std:: exception hierarchy for backward compatibility.
	/// Use catch(const MMLException&) to handle any MML library error.
	/// Call message() to get the error string.
	class MMLException {
	public:
		virtual ~MMLException() = default;

		/// @brief Get error message (delegates to std::exception::what()).
		const char* message() const noexcept {
			if (auto* e = dynamic_cast<const std::exception*>(this))
				return e->what();
			return "unknown MML error";
		}
	};

	//////////             Base error exceptions              ///////////
	/// @brief Generic argument validation error for invalid parameters.
	/// Use when function arguments violate preconditions (negative where positive required,
	/// empty containers where non-empty required, etc.)
	class ArgumentError : public std::invalid_argument, public MMLException
	{
	public:
		explicit ArgumentError(const std::string& message) 
			: std::invalid_argument(message) { }
	};

	/// @brief Error for operations outside valid domain.
	/// Use for mathematical domain violations (e.g., sqrt of negative, log of zero,
	/// parameter outside valid range).
	class DomainError : public std::domain_error, public MMLException
	{
	public:
		explicit DomainError(const std::string& message) 
			: std::domain_error(message) { }
	};

	/// @brief Error for division by zero.
	/// Use when division or modulo by zero is attempted.
	class DivisionByZeroError : public std::domain_error, public MMLException
	{
	public:
		explicit DivisionByZeroError(const std::string& message) 
			: std::domain_error(message) { }
	};

	/// @brief Error for unimplemented methods.
	/// Use for abstract base class methods or features not yet implemented.
	class NotImplementedError : public std::logic_error, public MMLException
	{
	public:
		explicit NotImplementedError(const std::string& message) 
			: std::logic_error(message) { }
	};

	/// @brief Generic index out of bounds error.
	/// Use for array/container index access violations.
	class IndexError : public std::out_of_range, public MMLException
	{
		int _index;
		int _size;
	public:
		IndexError(const std::string& message, int index = -1, int size = -1)
			: std::out_of_range(message), _index(index), _size(size) { }
		
		int index() const noexcept { return _index; }
		int size() const noexcept { return _size; }
	};

	//////////             Vector error exceptions            ///////////
	class VectorInitializationError : public std::invalid_argument, public MMLException
	{
		int _size1;
	public:
		VectorInitializationError(std::string inMessage, int size1) 
			: std::invalid_argument(inMessage), _size1(size1) { }
		
		int size() const noexcept { return _size1; }
	};
	
	class VectorDimensionError : public std::invalid_argument, public MMLException
	{
		int _size1, _size2;
	public:
		VectorDimensionError(std::string inMessage, int size1, int size2) 
			: std::invalid_argument(inMessage), _size1(size1), _size2(size2) { }
		
		int expected() const noexcept { return _size1; }
		int actual() const noexcept { return _size2; }
	};
	
	class VectorAccessBoundsError : public std::out_of_range, public MMLException
	{
		int _i, _n;
	public:
		VectorAccessBoundsError(std::string inMessage, int i, int n) 
			: std::out_of_range(inMessage), _i(i), _n(n) { }
		
		int index() const noexcept { return _i; }
		int size() const noexcept { return _n; }
	};

	//////////             Matrix error exceptions            ///////////
	class MatrixAllocationError : public std::out_of_range, public MMLException
	{
		int _rows, _cols;
	public:
		MatrixAllocationError(std::string inMessage, int rows, int cols) 
			: std::out_of_range(inMessage), _rows(rows), _cols(cols) { }
		
		int rows() const noexcept { return _rows; }
		int cols() const noexcept { return _cols; }
	};
	class MatrixAccessBoundsError : public std::out_of_range, public MMLException
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
	class MatrixDimensionError : public std::invalid_argument, public MMLException
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
	class SingularMatrixError : public std::domain_error, public MMLException
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
	class IntegrationTooManySteps : public std::domain_error, public MMLException
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
	class RealFuncInterpInitError: public std::domain_error, public MMLException
	{
	public:
		RealFuncInterpInitError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	class RealFuncInterpRuntimeError: public std::runtime_error, public MMLException
	{
	public:
		RealFuncInterpRuntimeError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	////////////             Tensor exceptions             /////////////
	class TensorCovarContravarNumError : public std::invalid_argument, public MMLException
	{
		int _numContra, _numCo;
	public:
		TensorCovarContravarNumError(std::string inMessage, int size1, int size2) 
			: std::invalid_argument(inMessage), _numContra(size1), _numCo(size2) { }
		
		int num_contravariant() const noexcept { return _numContra; }
		int num_covariant() const noexcept { return _numCo; }
	};
	class TensorCovarContravarArithmeticError : public std::invalid_argument, public MMLException
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
	class TensorIndexError : public std::invalid_argument, public MMLException
	{
	public:
		TensorIndexError(std::string inMessage) : std::invalid_argument(inMessage)
		{ }
	};    

	//////////           Root finding exceptions            ///////////
	class RootFindingError: public std::runtime_error, public MMLException
	{
	public:
		RootFindingError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////            ODE solver exceptions             ///////////
	class ODESolverError: public std::runtime_error, public MMLException
	{
	public:
		ODESolverError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////            Statistics exceptions             ///////////
	class StatisticsError : public std::runtime_error, public MMLException
	{
	public:
		StatisticsError(std::string inMessage) : std::runtime_error(inMessage)
		{
		}
	};

	//////////      Numerical method exceptions             ///////////
	class NumericalMethodError : public std::runtime_error, public MMLException
	{
	public:
		NumericalMethodError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};
	
	class MatrixNumericalError : public std::runtime_error, public MMLException
	{
	public:
		MatrixNumericalError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////          Convergence exceptions              ///////////
	class ConvergenceError : public std::runtime_error, public MMLException
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
	class GeometryError : public std::domain_error, public MMLException
	{
	public:
		GeometryError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	//////////           Quaternion exceptions              ///////////
	class QuaternionError : public std::domain_error, public MMLException
	{
	public:
		QuaternionError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	//////////             File I/O exceptions              ///////////
	class FileIOError : public std::runtime_error, public MMLException
	{
		std::string _filename;
	public:
		FileIOError(std::string inMessage, std::string filename = "") 
			: std::runtime_error(inMessage), _filename(filename) { }
		
		const std::string& filename() const noexcept { return _filename; }
	};

	//////////            Fourier exceptions                ///////////
	class FourierError : public std::invalid_argument, public MMLException
	{
	public:
		FourierError(std::string inMessage) : std::invalid_argument(inMessage)
		{ }
	};

	//////////          Visualization exceptions            ///////////
	class VisualizerError : public std::runtime_error, public MMLException
	{
	public:
		VisualizerError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};

	//////////        Curve Fitting exceptions              ///////////
	class CurveFittingError : public std::invalid_argument, public MMLException
	{
	public:
		CurveFittingError(std::string inMessage) : std::invalid_argument(inMessage)
		{ }
	};

	//////////             Data I/O exceptions              ///////////
	class DataError : public std::runtime_error, public MMLException
	{
	public:
		DataError(std::string inMessage) : std::runtime_error(inMessage)
		{ }
	};
	//////////      Numeric input validation exceptions     ///////////
	/// @brief Error for non-finite numeric values (NaN, Inf) passed to algorithms.
	/// Use when algorithm input or intermediate values are not finite.
	class NumericInputError : public std::domain_error, public MMLException
	{
	public:
		explicit NumericInputError(const std::string& message)
			: std::domain_error("Numeric input error: " + message) { }
	};
}

#endif // MML_EXCEPTIONS_H