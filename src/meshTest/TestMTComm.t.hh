//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/TestMTComm.t.hh
 * \author Randy M. Roberts
 * \date   Mon Aug 23 15:09:07 1999
 * \brief  Implementation file for the TestMTComm class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestMTComm.hh"
#include <iostream>
#include <functional>
#include <algorithm>

namespace rtt_meshTest
{

template<class T>
inline std::binder2nd<std::not_equal_to<T> > notEqualTo(T val)
{
    return std::bind2nd(std::not_equal_to<T>(), val);
}

template<class MTFactory>
void TestMTComm<MTFactory>::run()
{
    // Run the tests in this test class.

    os() << "Begin Running....... TestMTComm tests." << std::endl;

    setPassed(true);
    
    t3();

    t4();

    os() << "Completed Running... TestMTComm tests." << std::endl;
}

// Test the MT scatters and gathers

template<class MTFactory>
void TestMTComm<MTFactory>::t3()
{
    os() << "t3: beginning.\n";

    typedef MT::ccsf ccsf;
    typedef MT::fcdsf fcdsf;
    typedef MT::ncsf ncsf;
    typedef MT::vcsf vcsf;
    typedef MT::bssf bssf;

    // Create a field constructor with which to test the scatters and gathers.
    
    MTFactoryProduct meshProduct = meshFactory_m.create();
    FieldConstructor &fCtor = meshProduct.fieldConstructor();
    
    ccsf c(fCtor);
    fcdsf f(fCtor);
    fcdsf f2(fCtor);
    vcsf v(fCtor);
    ncsf n(fCtor);
    bssf b(fCtor);

    // Cell centered to face centered scatter

    bool passed;
    std::string test("t3() failed testing MT::scatter(f, c, MT::OpAssign()).");
    
    c = 1.0;
    f = 0.0;
    MT::scatter(f, c, MT::OpAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    test = "t3() failed testing MT::scatter(f, c, MT::OpAddAssign()).";

    c = 1.0;
    f = 0.0;
    b = 0.0;
    MT::scatter(f, c, MT::OpAddAssign());
    MT::gather(b, f, MT::OpAssign());
    MT::gather(f, b, MT::OpAddAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(2.0));
    testassert(passed, test, __FILE__, __LINE__);

    test = "t3() failed testing MT::scatter(f, c, MT::OpSubAssign()).";

    c = 1.0;
    f = 0.0;
    b = 0.0;
    MT::scatter(f, c, MT::OpSubAssign());
    MT::gather(b, f, MT::OpAssign());
    MT::gather(f, b, MT::OpAddAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(-2.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    test = "t3() failed testing MT::scatter(f, c, MT::OpMultAssign()).";

    c = 3.0;
    f = 1.0;
    b = 0.0;
    MT::scatter(f, c, MT::OpMultAssign());
    MT::gather(b, f, MT::OpAssign());
    MT::gather(f, b, MT::OpMultAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(9.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Face centered to cell centered scatter

    test = "t3() failed testing MT::scatter(c, f, MT::OpAddAssign()).";

    f = 1.0;
    c = 0.0;
    MT::scatter(c, f, MT::OpAddAssign());

    passed = c.end() == std::find_if(c.begin(), c.end(), notEqualTo(6.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    test = "t3() failed testing MT::scatter(c, f, MT::OpMultAssign()).";

    f = 3.0;
    c = 1.0;
    MT::scatter(c, f, MT::OpMultAssign());

    passed = c.end() == std::find_if(c.begin(), c.end(), notEqualTo(729.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Vertex centered to face centered scatter

    test = "t3() failed testing MT::scatter(f, v, MT::OpAddAssign()).";

    v = 1.0;
    f = 0.0;
    MT::scatter(f, v, MT::OpAddAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(4.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    test = "t3() failed testing MT::scatter(f, v, MT::OpMultAssign()).";

    v = 3.0;
    f = 1.0;
    MT::scatter(f, v, MT::OpMultAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(81.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Face centered to vertex centered scatter

    test = "t3() failed testing MT::scatter(v, f, MT::OpAddAssign()).";

    f = 1.0;
    v = 0.0;
    MT::scatter(v, f, MT::OpAddAssign());

    passed = v.end() == std::find_if(v.begin(), v.end(), notEqualTo(3.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Vertex centered to node centered scatter

    test = "t3() failed testing MT::scatter(n, v, MT::OpAssign()).";

    v = 1.0;
    n = 0.0;
    MT::scatter(n, v, MT::OpAssign());

    passed = n.end() == std::find_if(n.begin(), n.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    test = "t3() failed testing MT::scatter(n, v, MT::OpMultAssign()).";

    v = 1.0;
    n = 1.0;
    MT::scatter(n, v, MT::OpMultAssign());

    passed = n.end() == std::find_if(n.begin(), n.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Vertex centered to cell centered scatter

    test = "t3() failed testing MT::scatter(c, v, MT::OpAddAssign()).";

    v = 1.0;
    c = 0.0;
    MT::scatter(c, v, MT::OpAddAssign());

    passed = c.end() == std::find_if(c.begin(), c.end(), notEqualTo(8.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Cell centered to face centered gather

    test = "t3() failed testing MT::gather(f, c, MT::OpAssign()).";

    c = 1.0;
    f = 0.0;
    MT::gather(f, c, MT::OpAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    test = "t3() failed testing MT::gather(f, c, MT::OpMinAssign()).";

    c = 5.0;
    f = 3.0;
    MT::gather(f, c, MT::OpMinAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(3.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    c = 5.0;
    f = 6.0;
    MT::gather(f, c, MT::OpMinAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(5.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    test = "t3() failed testing MT::gather(f, c, MT::OpMaxAssign()).";

    c = 5.0;
    f = 3.0;
    MT::gather(f, c, MT::OpMaxAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(5.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    c = 5.0;
    f = 6.0;
    MT::gather(f, c, MT::OpMaxAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(6.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Face centered to boundary specified gather

    test = "t3() failed testing MT::gather(b, f, MT::OpAssign()).";

    f = 1.0;
    b = 0.0;
    MT::gather(b, f, MT::OpAssign());

    passed = b.end() == std::find_if(b.begin(), b.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Boundary specified to face centered gather

    test = "t3() failed testing MT::gather(f, b, MT::OpAssign()).";

    b = 1.0;
    f = 0.0;
    MT::gather(f, b, MT::OpAssign());
    b = 0.0;
    MT::gather(b, f, MT::OpAssign());

    passed = b.end() == std::find_if(b.begin(), b.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Node centered to vertex centered gather

    test = "t3() failed testing MT::gather(v, n, MT::OpAssign()).";

    n = 1.0;
    v = 0.0;
    MT::gather(v, n, MT::OpAssign());

    passed = v.end() == std::find_if(v.begin(), v.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    test = "t3() failed testing MT::gather(v, n, MT::OpAddAssign()).";

    n = 1.0;
    v = 0.0;
    MT::gather(v, n, MT::OpAddAssign());

    passed = v.end() == std::find_if(v.begin(), v.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Cell centered to vertex centered gather

    test = "t3() failed testing MT::gather(v, c, MT::OpAssign()).";

    c = 1.0;
    v = 0.0;
    MT::gather(v, c, MT::OpAssign());

    passed = v.end() == std::find_if(v.begin(), v.end(), notEqualTo(1.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    // Face swaps

    test = "t3() failed testing MT::swap_faces(f, f2).";
    
    f2 = 3.0;
    f = 2.0;
    MT::swap_faces(f, f2);
    b = 1.0;
    MT::gather(b, f, MT::OpAssign());

    passed = b.end() == std::find_if(b.begin(), b.end(), notEqualTo(0.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    b = 3.0;
    MT::gather(f, b, MT::OpAssign());

    passed = f.end() == std::find_if(f.begin(), f.end(), notEqualTo(3.0));
    testassert(passed, test, __FILE__, __LINE__);
    
    os() << "t3: end\n";
}



// Test the MT global reduction functions

template<class MTFactory>
void TestMTComm<MTFactory>::t4()
{
    os() << "t4: beginning.\n";

    typedef MT::ccif ccif;
    typedef MT::fcdif fcdif;
    typedef MT::ncif ncif;
    typedef MT::vcif vcif;
    typedef MT::bsif bsif;

    // Create a mesh and a field constructor with which to test
    // the scatters and gathers.
    
    MTFactoryProduct meshProduct = meshFactory_m.create();
    MT &mesh = meshProduct.mesh();
    FieldConstructor &fCtor = meshProduct.fieldConstructor();
    
    ccif c(fCtor);
    fcdif f(fCtor);
    vcif v(fCtor);
    ncif n(fCtor);
    bsif b(fCtor);

    int value = 3;

    // Cell centered sum

    std::string test("t4() failed testing MT::sum(c).");

    c = value;
    ccif::value_type csum = MT::sum(c);
    testassert(!(csum != value*mesh.get_total_ncells()), test,
	       __FILE__, __LINE__);
    
    // Face centered sum

    test = "t4() failed testing MT::sum(f).";
    
    f = value;
    fcdif::value_type fsum = MT::sum(f);
    testassert(!(fsum != value*6*mesh.get_total_ncells()), test,
	       __FILE__, __LINE__);
    
    // Node centered sum

    test = "t4() failed testing MT::sum(n).";
    
    n = value;
    ncif::value_type nsum = MT::sum(n);
    testassert(!(nsum !=
		 value*((mesh.get_ncx() + 1)*(mesh.get_ncy() + 1)*
			(mesh.get_ncz() + 1))), test, __FILE__, __LINE__);
    
    // Vertex centered sum

    test = "t4() failed testing MT::sum(v).";
    
    v = value;
    vcif::value_type vsum = MT::sum(v);
    testassert(!(vsum != value*8*mesh.get_total_ncells()), test,
	       __FILE__, __LINE__);
    
    // Boundary specified sum

    test = "t4() failed testing MT::sum(b).";
    
    b = value;
    bsif::value_type bsum = MT::sum(b);
    testassert(!(bsum !=
		 value*2*((mesh.get_ncx()*mesh.get_ncy()) +
			  (mesh.get_ncx()*mesh.get_ncz()) +
			  (mesh.get_ncy()*mesh.get_ncz()))),
	       test, __FILE__, __LINE__);
    
    // Cell centered minimum

    test = "t4() failed testing MT::min(c).";
    
    value = 1;
    for (ccif::iterator iter = c.begin(); iter != c.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    ccif::value_type cmin = MT::min(c);
    testassert(!(cmin != 1), test, __FILE__, __LINE__);
    
    // Face centered minimum

    test = "t4() failed testing MT::min(f).";
    
    value = 1;
    for (fcdif::iterator iter = f.begin(); iter != f.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    fcdif::value_type fmin = MT::min(f);
    testassert(!(fmin != 1), test, __FILE__, __LINE__);
    
    // Node centered minimum

    test = "t4() failed testing MT::min(n).";
    
    value = 1;
    for (ncif::iterator iter = n.begin(); iter != n.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    ncif::value_type nmin = MT::min(n);
    testassert(!(nmin != 1), test, __FILE__, __LINE__);
    
    // Vertex centered minimum

    test = "t4() failed testing MT::min(v).";
    
    value = 1;
    for (vcif::iterator iter = v.begin(); iter != v.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    vcif::value_type vmin = MT::min(v);
    testassert(!(vmin != 1), test, __FILE__, __LINE__);
    
    // Boundary specified minimum

    test = "t4() failed testing MT::min(b).";
    
    value = 1;
    for (bsif::iterator iter = b.begin(); iter != b.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    bsif::value_type bmin = MT::min(b);
    testassert(!(bmin != 1), test, __FILE__, __LINE__);

    // Cell centered maximum

    test = "t4() failed testing MT::max(c).";
    
    value = 1;
    for (ccif::iterator iter = c.begin(); iter != c.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    ccif::value_type cmax = MT::max(c);
    testassert(!(cmax != 1), test, __FILE__, __LINE__);
    
    // Face centered maximum

    test = "t4() failed testing MT::max(f).";
    
    value = 1;
    for (fcdif::iterator iter = f.begin(); iter != f.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    fcdif::value_type fmax = MT::max(f);
    testassert(!(fmax != 1), test, __FILE__, __LINE__);
    
    // Node centered maximum

    test = "t4() failed testing MT::max(n).";
    
    value = 1;
    for (ncif::iterator iter = n.begin(); iter != n.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    ncif::value_type nmax = MT::max(n);
    testassert(!(nmax != 1), test, __FILE__, __LINE__);
    
    // Vertex centered maximum

    test = "t4() failed testing MT::max(v).";
    
    value = 1;
    for (vcif::iterator iter = v.begin(); iter != v.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    vcif::value_type vmax = MT::max(v);
    testassert(!(vmax != 1), test, __FILE__, __LINE__);
    
    // Boundary specified maximum

    test = "t4() failed testing MT::max(b).";
    
    value = 1;
    for (bsif::iterator iter = b.begin(); iter != b.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    bsif::value_type bmax = MT::max(b);
    testassert(!(bmax != 1), test, __FILE__, __LINE__);
    
    os() << "t4: end\n";
}


} // end namespace rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of TestMTComm.t.hh
//---------------------------------------------------------------------------//
