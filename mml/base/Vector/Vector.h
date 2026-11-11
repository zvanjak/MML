///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Vector.h                                                            ///
///  Description: Generic vector class template for numerical computations            ///
///               Arithmetic operations, norms, dot/cross products                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_Vector_H
#define MML_Vector_H

#include <mml/base/PrintFormatting.h>
#include <mml/MMLBase.h>

// Standard headers - include what we use
#include <algorithm>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace MML
{
	/// @brief Generic vector class template for numerical computations
	/// @details Provides a dynamically-sized vector with:
	///          - Standard arithmetic operations (+, -, *, /)
	///          - Vector norms (L1, L2, L∞)
	///          - Element access (checked and unchecked)
	///          - STL-compatible iteration
	///
	///          Wraps std::vector<Type> with mathematical operations.
	///
	/// @threadsafety Thread-safe for const operations. Non-const operations require external
	///               synchronization. See docs/THREADING.md for details.
	///
	/// @tparam Type Element type (typically Real, Complex, or arithmetic types)
	template<class Type>
	class	Vector
	{
	private:
		std::vector<Type> _elems;

	public:
		typedef Type value_type;      ///< Element type alias for STL compatibility

		/////////////////////                  Constructors                ////////////////////
		
		/// @brief Default constructor - creates empty vector
		Vector() {}
		
		/// @brief Construct vector of given size
		/// @param n Number of elements (must be non-negative)
		/// @throws VectorInitializationError if n < 0
		/// @details Elements initialized to zero for arithmetic types, default-constructed otherwise
		explicit Vector(int n) {
			if(n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			// Initialize to zero for numeric types, default construct for others
			if constexpr (MMLArithmetic<Type>) {
				_elems.resize(n, Type{ 0 });
			} else {
				_elems.resize(n);
			}
		}
		
		/// @brief Construct vector of given size with uniform value
		/// @param n Number of elements (must be non-negative)
		/// @param val Value to fill all elements with
		/// @throws VectorInitializationError if n < 0
		explicit Vector(int n, const Type &val) {
			if (n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n, val);
		}
		
		/// @brief Construct vector from C-style array
		/// @param n Number of elements to copy (must be non-negative)
		/// @param vals Pointer to array of values
		/// @throws VectorInitializationError if n < 0
		explicit Vector(int n, Type* vals) 
		{
			if (n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n);
			for (int i = 0; i < n; ++i)
				_elems[i] = vals[i];
		}
		
		/// @brief Construct from std::vector
		explicit Vector(const std::vector<Type> &values) : _elems(values) {}
		
		/// @brief Construct from initializer list (e.g., Vector<double> v{1.0, 2.0, 3.0})
		explicit Vector(std::initializer_list<Type> list) : _elems(list) {}

		/// @brief Create a unit vector (all zeros except one element = 1)
		/// @param dimVec Dimension of the vector
		/// @param indUnit Index of the unit element (0-based)
		/// @return Unit vector e_i where e_i[indUnit] = 1
		/// @throws VectorDimensionError if indUnit is out of range
		static Vector UnitVector(int dimVec, int indUnit) requires MMLArithmetic<Type>
		{
			if (indUnit < 0 || indUnit >= dimVec)
				throw VectorDimensionError("Vector::UnitVector - wrong unit index", dimVec, indUnit);

			Vector ret(dimVec);
			ret[indUnit] = Type{ 1.0 };
			return ret;
		}

		// not really needed, but let's be explicit
		Vector(const Vector& b) = default; 
		Vector(Vector&& b) noexcept = default;
		Vector& operator=(const Vector& b) = default; 
		Vector& operator=(Vector&& b) noexcept = default;

		////////////////            std::vector forwarding                 ////////////////////
		
		/// @brief Get number of elements
		int  size()			const noexcept { return (int)_elems.size(); }
		
		/// @brief Check if vector is empty (STL-compatible)
		bool empty()	const noexcept { return _elems.empty(); }

		// front and back
		Type& front()							noexcept { return _elems.front(); }
		const Type& front() const noexcept { return _elems.front(); }
		Type& back()							noexcept { return _elems.back(); }
		const Type& back() const	noexcept { return _elems.back(); }

    auto begin()	noexcept { return _elems.begin(); }
    auto end()		noexcept { return _elems.end(); }
    auto begin()	const noexcept { return _elems.begin(); }
    auto end()		const noexcept { return _elems.end(); }
    auto cbegin() const noexcept { return _elems.cbegin(); }
    auto cend()		const noexcept { return _elems.cend(); }

		void push_back(const Type& val) { _elems.push_back(val); }
		void push_back(Type&& val)			{ _elems.push_back(std::move(val)); }
		
		void insert(int pos, const Type& val) {
			if (pos < 0 || pos > static_cast<int>(_elems.size()))
				throw VectorAccessBoundsError("Vector::insert - position out of range", pos, static_cast<int>(_elems.size()));
			_elems.insert(_elems.begin() + pos, val);
		}
		void insert(int pos, Type&& val) {
			if (pos < 0 || pos > static_cast<int>(_elems.size()))
				throw VectorAccessBoundsError("Vector::insert - position out of range", pos, static_cast<int>(_elems.size()));
			_elems.insert(_elems.begin() + pos, std::move(val));
		}
		
		void erase(int pos) {
			if (pos < 0 || pos >= static_cast<int>(_elems.size()))
				throw VectorAccessBoundsError("Vector::erase - position out of range", pos, static_cast<int>(_elems.size()));
			_elems.erase(_elems.begin() + pos);
		}
		void erase(int start, int end) {
			if (start < 0 || end > static_cast<int>(_elems.size()) || start > end)
				throw VectorAccessBoundsError("Vector::erase - range out of bounds", start, static_cast<int>(_elems.size()));
			_elems.erase(_elems.begin() + start, _elems.begin() + end);
		}
		void erase(const Type& val)			{ _elems.erase(std::remove(_elems.begin(), _elems.end(), val), _elems.end()); }

		/// @brief Remove all elements (STL spelling)
		void clear()	{ _elems.clear(); }
		/// @deprecated Use clear() (CONTRACT_POLICY §11: STL-shaped protocol is spelled as the STL spells it).
		[[deprecated("use clear()")]]
		void Clear()	{ clear(); }
		
		/// @brief Resize the vector (STL spelling)
		/// @param newLen New size
		/// @param preserveElements Ignored (std::vector::resize always preserves existing elements)
		void resize(int newLen, [[maybe_unused]] bool preserveElements = false)	
		{ 
			_elems.resize(newLen); 
		}

		/// @deprecated Use resize() (CONTRACT_POLICY §11: STL-shaped protocol is spelled as the STL spells it).
		[[deprecated("use resize()")]]
		void Resize(int newLen, bool preserveElements = false) { resize(newLen, preserveElements); }

		/////////////////////            Accessing elements             ///////////////////////
		
		/// @brief Unchecked element access
		inline Type&       operator[](int n)       noexcept { return _elems[n]; }
		/// @brief Unchecked element access (const)
		inline const Type& operator[](int n) const noexcept { return _elems[n]; }

		/// @brief Checked element access
		/// @param n Index (0-based)
		/// @return Reference to element at index n
		/// @throws VectorAccessBoundsError if index out of bounds
		Type& at(int n)	{
			if(n < 0 || n >= size())
				throw VectorAccessBoundsError("Vector::at - index out of bounds", n, size());
			else
				return _elems[n];
		}
		/// @brief Checked element access (const)
		/// @throws VectorAccessBoundsError if index out of bounds
		const Type& at(int n) const { 
			if(n < 0 || n >= size())
				throw VectorAccessBoundsError("Vector::at - index out of bounds", n, size());
			else
				return _elems[n];
		}

		///////////////////////          Arithmetic operators         ////////////////////////
		
		/// @brief Unary negation
		[[nodiscard]] Vector  operator-() const
		{
			Vector ret(size());
			for (int i = 0; i < size(); i++)
				ret._elems[i] = Type{ -1 } * (*this)[i];
			return ret;
		}
		
		/// @brief Vector addition
		/// @throws VectorDimensionError if sizes don't match
		[[nodiscard]] Vector  operator+(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator+() - vectors must be equal size", size(), b.size());

			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = (*this)[i] + b._elems[i];
			return ret;
		}
		
		/// @brief In-place vector addition
		/// @throws VectorDimensionError if sizes don't match
		Vector& operator+=(const Vector& b)
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator+=() - vectors must be equal size", size(), b.size());
			
			for (int i = 0; i < b.size(); i++)
				_elems[i] += b._elems[i];
			return *this;
		}
		
		/// @brief Vector subtraction
		/// @throws VectorDimensionError if sizes don't match
		[[nodiscard]] Vector  operator-(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator-() - vectors must be equal size", size(), b.size());

			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = (*this)[i] - b._elems[i];
			return ret;
		}
		
		/// @brief In-place vector subtraction
		/// @throws VectorDimensionError if sizes don't match
		Vector& operator-=(const Vector& b)
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator-=() - vectors must be equal size", size(), b.size());
			
			for (int i = 0; i < b.size(); i++)
				_elems[i] -= b._elems[i];
			return *this;
		}

		/// @brief Scalar multiplication (vector * scalar)
		[[nodiscard]] Vector  operator*(Type b) const
		{
			Vector ret(size());;
			for (int i = 0; i < size(); i++)
				ret._elems[i] = b * _elems[i];
			return ret;
		}
		
		/// @brief In-place scalar multiplication
		Vector& operator*=(Type b)
		{
			for (int i = 0; i < size(); i++)
				_elems[i] *= b;
			return *this;
		}
		
		/// @brief Scalar division (unchecked — caller must ensure b != 0)
		[[nodiscard]] Vector  operator/(Type b) const
		{
			Vector ret(size());
			for (int i = 0; i < size(); i++)
				ret._elems[i] = _elems[i] / b;
			return ret;
		}
		
		/// @brief In-place scalar division (unchecked — caller must ensure b != 0)
		Vector& operator/=(Type b)
		{
			for (int i = 0; i < size(); i++)
				_elems[i] /= b;
			return *this;
		}

		/// @brief Safe scalar division with zero-guard
		/// @param b Scalar divisor
		/// @return New vector with elements divided by b
		/// @throws DivisionByZeroError if |b| < DivisionSafetyThreshold
		[[nodiscard]] Vector safeDivide(Type b) const
		{
			if (std::abs(b) < Precision::DivisionSafetyThreshold)
				throw DivisionByZeroError("Vector::safeDivide - division by near-zero scalar");
			return *this / b;
		}
		
		/// @brief Scalar multiplication (scalar * vector)
		friend Vector operator*(Type a, const Vector& b)
		{
			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = a * b._elems[i];
			return ret;
		}

		///////////////////////             Testing equality             ////////////////////////
		
		/// @brief Exact equality comparison; vectors of different sizes are not equal
		bool operator==(const Vector& b) const
		{
			if (size() != b.size())
				return false;

			for (int i = 0; i < size(); i++)
				if ((*this)[i] != b[i])
					return false;

			return true;
		}
		
		/// @brief Exact inequality comparison
		bool operator!=(const Vector& b) const
		{
			return !(*this == b);
		}
		
		/// @brief Approximate equality with tolerance
		/// @param b Vector to compare with
		/// @param eps Maximum allowed difference per element
		/// @return true if |this[i] - b[i]| <= eps for all i
		/// @throws VectorDimensionError if sizes don't match
		bool IsEqualTo(const Vector& b, Real eps = Precision::VectorIsEqualTolerance) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::IsEqual - vectors must be equal size", size(), b.size());

			for (int i = 0; i < size(); i++)
			{
				if (Abs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}	
		
		/// @brief Check if all elements are zero
		bool IsZero() const noexcept
		{
			for (int i = 0; i < size(); i++)
				if (Abs((*this)[i]) != 0.0 )
					return false;
			return true;
		}

		/// @brief Check if all elements are near zero within tolerance
		/// @param eps Maximum absolute value for an element to be considered near-zero
		/// @return true if |this[i]| <= eps for all i
		bool IsNearZero(Real eps = Precision::VectorIsEqualTolerance) const noexcept
		{
			for (int i = 0; i < size(); i++)
				if (Abs((*this)[i]) > eps)
					return false;
			return true;
		}

		/// @deprecated Use IsZero() (CONTRACT_POLICY §11: mathematical predicates are PascalCase).
		[[deprecated("use IsZero()")]]
		bool isZero() const noexcept { return IsZero(); }
		/// @deprecated Use IsNearZero().
		[[deprecated("use IsNearZero()")]]
		bool isNearZero(Real eps = Precision::VectorIsEqualTolerance) const noexcept { return IsNearZero(eps); }

		//////////////////////                 Operations                 ///////////////////////
		
		/// @brief L1 norm (sum of absolute values): ||v||_1 = Σ|v_i|
		Real NormL1() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < size(); i++)
				norm += Abs((*this)[i]);
			return norm;
		}
		
		/// @brief L2 norm (Euclidean length): ||v||_2 = √(Σ|v_i|²)
		/// @details Scaled accumulation (StableNormL2Accumulator): no overflow for components up to
		///          Real max and no underflow to 0 for tiny components; |z| is used for complex entries.
		Real NormL2() const
		{
			StableNormL2Accumulator acc;
			for (int i = 0; i < size(); i++)
				acc.add(Abs((*this)[i]));
			return acc.value();
		}
		
		/// @brief L∞ norm (maximum absolute value): ||v||_∞ = max|v_i|
		Real NormLInf() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < size(); i++)
				norm = std::max(norm, Abs((*this)[i]));
			return norm;
		}

		///////////////////////////               I/O                 ///////////////////////////
		
		/// @brief Print vector to stream with formatting
		/// @param stream Output stream
		/// @param width Field width for each element
		/// @param precision Decimal precision
		/// @param zeroThreshold Values with |x| <= threshold printed as 0
		std::ostream& Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			return PrintFormatting::Vector(stream, size(), [&](int i) { return (*this)[i]; }, width, precision, zeroThreshold, true);
		}
		
		/// @brief Print vector to stream with formatting
		std::ostream& Print(std::ostream& stream, int width, int precision) const
		{
			return Print(stream, width, precision, 0.0);
		}
		
		/// @brief Print vector with message prefix and newline
    std::ostream& PrintLine(std::ostream& stream, const std::string &msg, int width, int precision) const
		{
			stream << msg;
			Print(stream, width, precision);
			stream << std::endl;

			return stream;
		}
		
		/// @brief Print vector as column (one element per line)
		std::ostream& PrintCol(std::ostream& stream, int width, int precision) const
		{
			return PrintFormatting::Column(stream, size(), [&](int i) { return (*this)[i]; }, width, precision);
		}		

		/// @brief Convert to string representation
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		
		/// @brief Convert to string with default formatting
		std::string to_string() const
		{
			return to_string(Defaults::VectorPrintWidth, Defaults::VectorPrintPrecision);
		}

		/// @brief Stream output operator
		friend std::ostream& operator<<(std::ostream &stream, const Vector &a)
		{
			a.Print(stream, Defaults::VectorPrintWidth, Defaults::VectorPrintPrecision);

			return stream;
		}
		
	};

	// Common type aliases
	typedef Vector<int>     VectorInt;      ///< Integer vector
	typedef Vector<float>   VectorFlt;      ///< Single-precision vector
	typedef Vector<double>  VectorDbl;      ///< Double-precision vector
	typedef Vector<Complex> VectorComplex;  ///< Complex vector

	// Short aliases
	typedef Vector<int>     VecI;   ///< Short alias for VectorInt
	typedef Vector<float>   VecF;   ///< Short alias for VectorFlt
	typedef Vector<double>  VecD;   ///< Short alias for VectorDbl
	typedef Vector<Complex> VecC;   ///< Short alias for VectorComplex

	// Verify noexcept move operations enable STL optimizations
	static_assert(std::is_nothrow_move_constructible_v<Vector<double>>,
	              "Vector<double> should be nothrow move constructible");
	static_assert(std::is_nothrow_move_assignable_v<Vector<double>>,
	              "Vector<double> should be nothrow move assignable");
}

#endif // MML_Vector_H