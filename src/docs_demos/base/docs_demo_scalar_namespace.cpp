#include <mml/MMLTypeDefs.h>
#include <mml/base/Vector/Vector.h>
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>

#include <iostream>
#include <stdexcept>
#include <type_traits>

void Docs_Demo_ScalarNamespace()
{
	MML::Real tolerance = MML::Precision::LooseTolerance;
	MML::Complex z(MML::Real(3), MML::Real(4));
	MML::Vector<MML::Real> values{MML::Real(3), MML::Real(4)};
	if (!MML::IsNearlyEqual(values.NormL2(), MML::Real(5), tolerance, tolerance))
		throw std::runtime_error("Namespaced scalar numerical behavior changed");
	{
		using MML::Real;
		using MML::Complex;
		const Real value = static_cast<Real>(1.0L / 3.0L);
		const Complex local(value, Real(0));
		std::cout << "Locally imported scalar: " << local << '\n';
	}
	static_assert(std::is_same_v<MML::Complex, std::complex<MML::Real>>);
	const MML::Real precise = static_cast<MML::Real>(3.141592653589793238462643383279502884L);
	std::cout << "Qualified scalars: |" << z << "|=" << MML::Abs(z)
	          << ", vector norm=" << values.NormL2() << ", typed pi=" << precise << '\n';
}
