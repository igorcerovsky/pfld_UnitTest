#include "stdafx.h"
#include "CppUnitTest.h"
#include <algorithm>

#include "../pfld/facet_t.hpp"
#include "../pfld/facet_t.cpp"


using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest_pfd
{

	template<typename T>
	void AssertPoints(pfld::Point3D<T> pt1, pfld::Point3D<T> pt2, T eps)
	{
		Assert::AreEqual(pt1.x, pt2.x, eps);
		Assert::AreEqual(pt1.y, pt2.y, eps);
		Assert::AreEqual(pt1.z, pt2.z, eps);
	};
	template<typename T>
	void TestBody(const T eps)
	{
		using flt = T;
		using facet = pfld::Facet < flt > ;
		using point = facet::point;
		using ptvec = facet::ptvec;
		facet fct;
		flt a{ 1000. };
		ptvec v{ point(0, 0, 0), point(0, 0, -a), point(a, 0, -a), point(0, a, -a) };
		std::vector<ptvec> vf{
			{ v[0], v[1], v[2] },
			{ v[0], v[2], v[3] },
			{ v[0], v[3], v[1] },
			{ v[1], v[3], v[2] } };
		std::vector<facet> facets(4);
		auto itVf = vf.begin();
		for (auto it = facets.begin(); it != facets.end(); ++it, ++itVf)
			it->Init(*itVf);

		point r(500, 500, 1), g;
		point M{ 1, 10, 100 }, gGS, m;
		flt gzGS{ 0. }, gzV{ 0. };
		for (auto it = facets.begin(); it != facets.end(); ++it)
		{
			it->Fld_G(r, g);
			it->Fld_Gz(r, gzV);
			it->FldGS(r, M, m, gGS);
			it->FldGS_Gz(r, gzGS);
		}

		const point result(T(5.1341030021201644e-009),
			T(5.1341030021201644e-009),
			T(1.2401175118216113e-008));
		AssertPoints(g, result, eps);
		AssertPoints(gGS, result, eps);
		AssertPoints(g, result, eps);
		Assert::AreEqual(g.z, gGS.z, eps);
		Assert::AreEqual(gGS.z, gzGS, eps);
	}


	TEST_CLASS(UnitTest2)
	{
	public:

		TEST_METHOD(Test_Body_T_float)
		{
			TestBody<float>(1.0e-13f);
		}

		TEST_METHOD(Test_Body_T_double)
		{
			TestBody<double>(1.0e-20);
		}
	};
}