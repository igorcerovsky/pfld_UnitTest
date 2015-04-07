#include "stdafx.h"
#include "CppUnitTest.h"

#include "../pfld/facet.hpp"
#include "../pfld/facet.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest_pfd
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(Test_Point3D_Unit)
		{
			using point = pfld::Point3D < double > ;
			point pt1(0.0, 2.0, 0.0);
			pt1.Unit();
			Assert::IsTrue(pt1 == point(0,1,0));
			point pt2(3.0, 0.0, 0.0);
			pt2.Unit();
			Assert::IsTrue(pt2==point(1,0,0));
			point pt3(0.0, 0.0, 4.0);
			pt3.Unit();
			Assert::IsTrue(pt3==point(0,0,1));
		}

		TEST_METHOD(Test_Facet)
		{
			using point = pfld::Point3D < double >;
			using ptvec = pfld::ptvec;
			pfld::Facet fct;
			ptvec v{ point(0, 0, -1000), point(1000, 0, 0), point(0, 1000, 0) };
			fct.Init( v );
			point r(0, 0, 1), g, g2, g3;
			fct.FldVlado(r, g);
			const double resval = 4.8207079871718046e-008;
			point result(-resval, -resval, resval);
			Assert::IsTrue(g.IsEqualEps(result));

			pfld::Facet fctCopy(fct);
			fctCopy.FldVlado(r, g2);
			Assert::IsTrue(g2.IsEqualEps(result));

			pfld::Facet fctAssign = fct;
			fctAssign.FldVlado(r, g3);
			Assert::IsTrue(g3.IsEqualEps(result));
		}
	};
}