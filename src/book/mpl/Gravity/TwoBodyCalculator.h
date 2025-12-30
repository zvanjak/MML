#if !defined MPL_TWO_BODY_CALCULATOR_H
#define MPL_TWO_BODY_CALCULATOR_H

#include "MMLBase.h"

#include "base/VectorTypes.h"
#include "base/Geometry3D.h"

#include "GravityBase.h"

#include "Base/PhysicalConstants.h"

using namespace MML;

namespace MPL
{
	class TwoBodyGravityCalculator
	{
		public:
		TwoBodyGravityCalculator() = default;

		static Vec3Cart calculateForce(Real G, double mass1, const Vec3Cart& position1, double mass2, const Vec3Cart& position2 )
		{
			Vec3Cart r = position2 - position1;
			double distance = r.NormL2();
			
			if (distance == 0.0) return Vec3Cart(0.0, 0.0, 0.0);					// Avoid division by zero
			
			double forceMagnitude = (G * mass1 * mass2) / (distance * distance);
			
			return r.GetAsUnitVector() * forceMagnitude;
		}
		static Vec3Cart calculateForce(double mass1, const Vec3Cart& position1, double mass2, const Vec3Cart& position2)
		{
			return calculateForce(PhyConst::GravityConstant, mass1, position1, mass2, position2);
		}

		static Vec3Cart calculateForce(Real G, const GravityBodyState& body1, const GravityBodyState& body2)
		{
			return calculateForce(G, body1.Mass(), body1.R(), body2.Mass(), body2.R());
		}
		static Vec3Cart calculateForce(const GravityBodyState &body1, const GravityBodyState &body2)
		{
			return calculateForce(PhyConst::GravityConstant, body1.Mass(), body1.R(), body2.Mass(), body2.R());
		}

		// get trajectory type - ellipse, parabolic, hyperbolic
		static std::string GetTrajectoryType(Real G, const GravityBodyState& body1, const GravityBodyState& body2)
		{
			Real reducedMass = body1.Mass() * body2.Mass() / (body1.Mass() + body2.Mass());

			Real distance = (body2.R() - body1.R()).NormL2();
			Real velocity1 = body1.V().NormL2();
			Real velocity2 = body2.V().NormL2();

			Real totalEnergy = 0.5 * reducedMass * (POW2(velocity1) + POW2(velocity2)) - G * body1.Mass() * body2.Mass() / distance;

			if (totalEnergy < 0)
				return "Ellipse";
			else if (totalEnergy == 0)
				return "Parabolic";
			else
				return "Hyperbolic";
		}

		// get plane of motion
		static Plane3D GetPlaneOfMotion(Real G, const GravityBodyState& body1, const GravityBodyState& body2)
		{
			Vec3Cart r = body2.R() - body1.R();
			Vec3Cart v = body2.V() - body1.V();
			
			Vec3Cart normal = VectorProduct(r, v).GetAsUnitVector();
			
			return Plane3D(Pnt3Cart(r.X(),r.Y(), r.Z()), normal);
		}

		// CHECK!!!! calculate eccentricity vector
		//Vec3Cart CalculateEccentricity(const TwoBodyGravitySimConfig &sysConfig)
		//{
		//	Real mu = sysConfig.Mass1() * sysConfig.Mass2() / (sysConfig.Mass1() + sysConfig.Mass2());
		//	Vec3Cart r1 = sysConfig._initState.Body1().R();
		//	Vec3Cart r2 = sysConfig._initState.Body2().R();
		//	Vec3Cart v1 = sysConfig._initState.Body1().V();
		//	Vec3Cart v2 = sysConfig._initState.Body2().V();
		//	Vec3Cart r12 = r2 - r1;
		//	Vec3Cart v12 = v2 - v1;
		//	
		//	Real energy = 0.5 * mu * (POW2(v1.NormL2()) + POW2(v2.NormL2())) - sysConfig._G * sysConfig.Mass1() * sysConfig.Mass2() / r12.NormL2();
		//	
		//	Real h = VectorProduct(r12, v12).NormL2(); // specific angular momentum
		//	
		//	return VectorProduct(v12, h) / (mu * energy);
		//}
	};

	class TwoBodyGravityCalculator2D
	{
		public:
		TwoBodyGravityCalculator2D() = default;

		static Vec2Cart calculateForce(Real G, const Vec2Cart& position1, const Vec2Cart& position2, double mass1, double mass2)
		{
			Vec2Cart r = position2 - position1;
			double distance = r.NormL2();
			
			if (distance == 0.0) return Vec2Cart(0.0, 0.0);					// Avoid division by zero
			
			double forceMagnitude = (G * mass1 * mass2) / (distance * distance);
			
			return r.GetAsUnitVector() * forceMagnitude;
		}
		static Vec2Cart calculateForce(const Vec2Cart& position1, const Vec2Cart& position2, double mass1, double mass2)
		{
			return calculateForce(PhyConst::GravityConstant, position1, position2, mass1, mass2);
		}
		static Vec2Cart calculateForce(Real G, const GravityBodyState2D& body1, const GravityBodyState2D& body2)
		{
			return calculateForce(G, body1.R(), body2.R(), body1.Mass(), body2.Mass());
		}
		static Vec2Cart calculateForce(const GravityBodyState2D &body1, const GravityBodyState2D &body2)
		{
			return calculateForce(PhyConst::GravityConstant, body1.R(), body2.R(), body1.Mass(), body2.Mass());
		}

		static std::string GetTrajectoryType(Real G, const GravityBodyState2D& body1, const GravityBodyState2D& body2)
		{
			Real reducedMass = body1.Mass() * body2.Mass() / (body1.Mass() + body2.Mass());
			
			Real distance = (body2.R() - body1.R()).NormL2();
			Real velocity1 = body1.V().NormL2();
			Real velocity2 = body2.V().NormL2();
			
			Real totalEnergy = 0.5 * (body1.Mass() * POW2(velocity1) + body2.Mass() * POW2(velocity2)) - G * body1.Mass() * body2.Mass() / distance;
			
			if (totalEnergy < 0)
				return "Ellipse";
			else if (totalEnergy == 0)
				return "Parabolic";
			else
				return "Hyperbolic";
		}
	};
}

#endif // MPL_TWO_BODY_CALCULATOR_H