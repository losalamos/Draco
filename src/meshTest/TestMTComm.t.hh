//----------------------------------*-C++-*----------------------------------//
// TestMTComm.t.hh
// Randy M. Roberts
// Mon Aug 23 15:09:07 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Test the MT Model for its communication methods as an MT Concept.
//---------------------------------------------------------------------------//

#include "TestMTComm.hh"
#include <iostream>

namespace rtt_meshTest
{
 
template<class MTFactory>
void TestMTComm<MTFactory>::run()
{
    // Run the tests in this test class.

    os_m << "Begin Running....... TestMTComm tests." << std::endl;

    passed_m = true;
    
    t3();

    t4();

    os_m << "Completed Running... TestMTComm tests." << std::endl;
}

// Test the MT scatters and gathers

template<class MTFactory>
void TestMTComm<MTFactory>::t3()
{
    os_m << "t3: beginning.\n";

    typedef MT::ccsf ccsf;
    typedef MT::fcdsf fcdsf;
    typedef MT::ncsf ncsf;
    typedef MT::vcsf vcsf;
    typedef MT::bssf bssf;

    // need a temporary passed variable for each subset of tests.

    bool passed = true;

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

    c = 1.0;
    f = 0.0;
    MT::scatter(f, c, MT::OpAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(f, c, MT::OpAssign())."
	     << std::endl;

    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    c = 1.0;
    f = 0.0;
    b = 0.0;
    MT::scatter(f, c, MT::OpAddAssign());
    MT::gather(b, f, MT::OpAssign());
    MT::gather(f, b, MT::OpAddAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 2.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(f, c, MT::OpAddAssign()) ..."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    c = 1.0;
    f = 0.0;
    b = 0.0;
    MT::scatter(f, c, MT::OpSubAssign());
    MT::gather(b, f, MT::OpAssign());
    MT::gather(f, b, MT::OpAddAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != -2.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(f, c, MT::OpSubAssign()) ..."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    c = 3.0;
    f = 1.0;
    b = 0.0;
    MT::scatter(f, c, MT::OpMultAssign());
    MT::gather(b, f, MT::OpAssign());
    MT::gather(f, b, MT::OpMultAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 9.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(f, c, MT::OpMultAssign()) ..."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Face centered to cell centered scatter

    f = 1.0;
    c = 0.0;
    MT::scatter(c, f, MT::OpAddAssign());
    for (ccsf::const_iterator citer = c.begin(); citer != c.end();
         ++citer)
      if (*citer != 6.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(c, f, MT::OpAddAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    f = 3.0;
    c = 1.0;
    MT::scatter(c, f, MT::OpMultAssign());
    for (ccsf::const_iterator citer = c.begin(); citer != c.end();
         ++citer)
      if (*citer != 729.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(c, f, MT::OpMultAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Vertex centered to face centered scatter

    v = 1.0;
    f = 0.0;
    MT::scatter(f, v, MT::OpAddAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 4.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(f, v, MT::OpAddAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    v = 3.0;
    f = 1.0;
    MT::scatter(f, v, MT::OpMultAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 81.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(f, v, MT::OpMultAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Face centered to vertex centered scatter

    f = 1.0;
    v = 0.0;
    MT::scatter(v, f, MT::OpAddAssign());
    for (vcsf::const_iterator citer = v.begin(); citer != v.end();
         ++citer)
      if (*citer != 3.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(v, f, MT::OpAddAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Vertex centered to node centered scatter

    v = 1.0;
    n = 0.0;
    MT::scatter(n, v, MT::OpAssign());
    for (ncsf::const_iterator citer = n.begin(); citer != n.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(n, v, MT::OpAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    v = 1.0;
    n = 1.0;
    MT::scatter(n, v, MT::OpMultAssign());
    for (ncsf::const_iterator citer = n.begin(); citer != n.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(n, v, MT::OpMultAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Vertex centered to cell centered scatter

    v = 1.0;
    c = 0.0;
    MT::scatter(c, v, MT::OpAddAssign());
    for (ccsf::const_iterator citer = c.begin(); citer != c.end();
         ++citer)
      if (*citer != 8.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::scatter(c, v, MT::OpAddAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Cell centered to face centered gather

    c = 1.0;
    f = 0.0;
    MT::gather(f, c, MT::OpAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::gather(f, c, MT::OpAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    c = 5.0;
    f = 3.0;
    MT::gather(f, c, MT::OpMinAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 3.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "(1) MT::gather(f, c, MT::OpMinAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    c = 5.0;
    f = 6.0;
    MT::gather(f, c, MT::OpMinAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 5.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "(2) MT::gather(f, c, MT::OpMinAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    c = 5.0;
    f = 3.0;
    MT::gather(f, c, MT::OpMaxAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 5.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "(1) MT::gather(f, c, MT::OpMaxAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    c = 5.0;
    f = 6.0;
    MT::gather(f, c, MT::OpMaxAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 6.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "(2) MT::gather(f, c, MT::OpMaxAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Face centered to boundary specified gather

    f = 1.0;
    b = 0.0;
    MT::gather(b, f, MT::OpAssign());
    for (bssf::const_iterator citer = b.begin(); citer != b.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::gather(b, f, MT::OpAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Boundary specified to face centered gather

    b = 1.0;
    f = 0.0;
    MT::gather(f, b, MT::OpAssign());
    b = 0.0;
    MT::gather(b, f, MT::OpAssign());
    for (bssf::const_iterator citer = b.begin(); citer != b.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::gather(f, b, MT::OpAssign()) ..."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Node centered to vertex centered gather

    n = 1.0;
    v = 0.0;
    MT::gather(v, n, MT::OpAssign());
    for (vcsf::const_iterator citer = v.begin(); citer != v.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::gather(v, n, MT::OpAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    n = 1.0;
    v = 0.0;
    MT::gather(v, n, MT::OpAddAssign());
    for (vcsf::const_iterator citer = v.begin(); citer != v.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::gather(v, n, MT::OpAddAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Cell centered to vertex centered gather

    c = 1.0;
    v = 0.0;
    MT::gather(v, c, MT::OpAssign());
    for (vcsf::const_iterator citer = v.begin(); citer != v.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::gather(v, c, MT::OpAssign())."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Face swaps

    f2 = 3.0;
    f = 2.0;
    MT::swap_faces(f, f2);
    b = 1.0;
    MT::gather(b, f, MT::OpAssign());
    for (bssf::const_iterator citer = b.begin(); citer != b.end();
         ++citer)
      if (*citer != 0.0)
        passed = false;
    b = 3.0;
    MT::gather(f, b, MT::OpAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 3.0)
        passed = false;

    if (!passed)
	os_m << "t3() failed testing "
	     << "MT::swap_faces(f, f2)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    os_m << "t3: end\n";
}



// Test the MT global reduction functions

template<class MTFactory>
void TestMTComm<MTFactory>::t4()
{
    os_m << "t4: beginning.\n";

    typedef MT::ccif ccif;
    typedef MT::fcdif fcdif;
    typedef MT::ncif ncif;
    typedef MT::vcif vcif;
    typedef MT::bsif bsif;

    // need a temporary passed variable for each subset of tests.

    bool passed = true;

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

    c = value;
    ccif::value_type csum = MT::sum(c);
    if (csum != value*mesh.get_total_ncells())
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::sum(c)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Face centered sum

    f = value;
    fcdif::value_type fsum = MT::sum(f);
    if (fsum != value*6*mesh.get_total_ncells())
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::sum(f)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Node centered sum

    n = value;
    ncif::value_type nsum = MT::sum(n);
    if (nsum !=
        value*((mesh.get_ncx() + 1)*(mesh.get_ncy() + 1)*(mesh.get_ncz() + 1)))
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::sum(n)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Vertex centered sum

    v = value;
    vcif::value_type vsum = MT::sum(v);
    if (vsum != value*8*mesh.get_total_ncells())
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::sum(v)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Boundary specified sum

    b = value;
    bsif::value_type bsum = MT::sum(b);
    if (bsum !=
        value*2*((mesh.get_ncx()*mesh.get_ncy()) +
             (mesh.get_ncx()*mesh.get_ncz()) +
             (mesh.get_ncy()*mesh.get_ncz())))
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::sum(b)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Cell centered minimum

    value = 1;
    for (ccif::iterator iter = c.begin(); iter != c.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    ccif::value_type cmin = MT::min(c);
    if (cmin != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::min(c)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Face centered minimum

    value = 1;
    for (fcdif::iterator iter = f.begin(); iter != f.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    fcdif::value_type fmin = MT::min(f);
    if (fmin != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::min(f)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Node centered minimum

    value = 1;
    for (ncif::iterator iter = n.begin(); iter != n.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    ncif::value_type nmin = MT::min(n);
    if (nmin != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::min(n)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Vertex centered minimum

    value = 1;
    for (vcif::iterator iter = v.begin(); iter != v.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    vcif::value_type vmin = MT::min(v);
    if (vmin != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::min(v)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Boundary specified minimum

    value = 1;
    for (bsif::iterator iter = b.begin(); iter != b.end(); ++iter)
    {
        *iter = value;
        ++value;
    }
    bsif::value_type bmin = MT::min(b);
    if (bmin != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::min(b)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Cell centered maximum

    value = 1;
    for (ccif::iterator iter = c.begin(); iter != c.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    ccif::value_type cmax = MT::max(c);
    if (cmax != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::max(c)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Face centered maximum

    value = 1;
    for (fcdif::iterator iter = f.begin(); iter != f.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    fcdif::value_type fmax = MT::max(f);
    if (fmax != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::max(f)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Node centered maximum

    value = 1;
    for (ncif::iterator iter = n.begin(); iter != n.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    ncif::value_type nmax = MT::max(n);
    if (nmax != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::max(n)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Vertex centered maximum

    value = 1;
    for (vcif::iterator iter = v.begin(); iter != v.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    vcif::value_type vmax = MT::max(v);
    if (vmax != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::max(v)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    // Boundary specified maximum

    value = 1;
    for (bsif::iterator iter = b.begin(); iter != b.end(); ++iter)
    {
        *iter = value;
        --value;
    }
    bsif::value_type bmax = MT::max(b);
    if (bmax != 1)
        passed = false;

    if (!passed)
	os_m << "t4() failed testing "
	     << "MT::max(b)."
	     << std::endl;
    
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;
    
    os_m << "t4: end\n";
}


} // end namespace rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of TestMTComm.t.hh
//---------------------------------------------------------------------------//
