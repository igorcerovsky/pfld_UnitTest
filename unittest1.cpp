#include "stdafx.h"
#include "CppUnitTest.h"

#include "../pfld/facet.hpp"
#include "../pfld/facet.cpp"
#include "../pfld/pfld_compute.hpp"
#include "../pfld/pfld_compute.cpp"
#include "../pfld/pfld_test_io.h"
#include "../pfld/pfld_test_io.cpp"
#include "../pfld/pfld_compute_prll_ft.cpp"


using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest_pfd
{		
	template<class T>
	typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
		almost_equal(T x, T y, int ulp)
	{
		return std::abs(x - y) < std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
			|| std::abs(x - y) < std::numeric_limits<T>::min();
	}

	void Compute(void(*FieldFn)(pfld::facet_vec&, pfld::ptvec&, pfld::valvec&),
		pfld::facet_vec& facets, pfld::ptvec& fldPoints, pfld::valvec& outFld,
		std::string message);

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
			fct.Init(v);
			point r(0, 0, 1), g, g2, g3;
			fct.Fld_G(r, g);
			const double resval = 4.8207079871718046e-008;
			point result(-resval, -resval, resval);
			Assert::IsTrue(g.IsEqualEps(result));

			pfld::Facet fctCopy(fct);
			fctCopy.Fld_G(r, g2);
			Assert::IsTrue(g2.IsEqualEps(result));

			pfld::Facet fctAssign = fct;
			fctAssign.Fld_G(r, g3);
			Assert::IsTrue(g3.IsEqualEps(result));
		}

		TEST_METHOD(Test_Facet_Lin0)
		{
			using point = pfld::Point3D < double >;
			using ptvec = pfld::ptvec;
			pfld::Facet fct;
			ptvec v{ point(0, 0, -1000), point(1000, 0, 0), point(0, 1000, 0) };
			const double ro0{ 1. };
			point ro{ 0., 0., 0. };
			const double resval = 4.8207079871718046e-008;
			point result(-resval, -resval, resval);

			fct.Init(v);
			point r(0, 0, 1), g, g2, g3;
			fct.Fld_G(r, ro, ro0, g);
			Assert::IsTrue(g.IsEqualEps(result));

			pfld::Facet fctCopy(fct);
			fctCopy.Fld_G(r, g2);
			Assert::IsTrue(g2.IsEqualEps(result));

			pfld::Facet fctAssign = fct;
			fctAssign.Fld_G(r, g3);
			Assert::IsTrue(g3.IsEqualEps(result));
		}

		TEST_METHOD(Test_Facet_Lin)
		{
			using point = pfld::Point3D < double >;
			using ptvec = pfld::ptvec;
			pfld::Facet fct;
			const double ro0{ 1000.};
			point ro{ 0., 0., 1. };
			ptvec v{ point(0, 0, -1000), point(1000, 0, 0), point(0, 1000, 0) };
			fct.Init(v);
			point r(0, 0, 1), g, g2, g3;
			fct.Fld_G(r, ro, ro0, g);
			point result(-3.2142476436014269e-005, -3.2142476436014269e-005, 5.6270119911809142e-005);
			Assert::IsTrue(g.IsEqualEps(result));

			pfld::Facet fctCopy(fct);
			fctCopy.Fld_G(r, ro, ro0, g2);
			Assert::IsTrue(g2.IsEqualEps(result));

			pfld::Facet fctAssign = fct;
			fctAssign.Fld_G(r, ro, ro0, g3);
			Assert::IsTrue(g3.IsEqualEps(result));

			pfld::Facet fctGz = fct;
			double gz=0.0;
			fctGz.Fld_Gz(r, ro, ro0, gz);
			Assert::IsTrue(almost_equal(gz, result.z, 4));
		}

		const char* file_facets =  "./../../../pfld_UnitTest/test_data/pfld_test_facets.txt";
		const char* file_points =  "./../../../pfld_UnitTest/test_data/pfld_test_points.txt";
		const char* file_results = "./../../../pfld_UnitTest/test_data/pfld_test_results.txt";
		const int max_facets_to_load = 100;
		const int max_points_to_load = 100;

		TEST_METHOD(Test_Facet_Parallell)
		{
			pfld::facet_vec facets;
			pfld::GetFacets(facets, file_facets, max_facets_to_load, false);

			pfld::ptvec fldPts;
			pfld::GetFieldPoints(fldPts, file_points, max_points_to_load, false);
			void(*FieldFn)(pfld::facet_vec&, pfld::ptvec&, pfld::valvec&);
			pfld::valvec outFld; // uninitialized for this version
			FieldFn = pfld::Field_Gz_;
			Compute(FieldFn, facets, fldPts, outFld, "computing facets parallel future approach...");

			pfld::valvec res;
			pfld::LoadResults(res, file_results, (int)outFld.size());
			auto itCmp = outFld.begin();
			for (auto it = res.begin(); it != res.end(); ++it, ++itCmp)
			{
				bool bEQ = almost_equal(*it, *itCmp, 4);
				Assert::IsTrue(bEQ);
			}

			pfld::valvec outFld1(fldPts.size()); // uninitialized for this version
			FieldFn = pfld::Field_Gz;
			Compute(FieldFn, facets, fldPts, outFld1, "computing facets parallel approach...");
			itCmp = outFld1.begin();
			for (auto it = res.begin(); it != res.end(); ++it, ++itCmp)
			{
				bool bEQ = almost_equal(*it, *itCmp, 4);
				Assert::IsTrue(bEQ);
			}
		}


	};



	void Compute(void(*FieldFn)(pfld::facet_vec&, pfld::ptvec&, pfld::valvec&),
		pfld::facet_vec& facets, pfld::ptvec& fldPoints, pfld::valvec& outFld,
		std::string message)
	{
		using namespace std::chrono;
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		FieldFn(facets, fldPoints, outFld);
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	}
} // namespace