//----------------------------------*-C++-*----------------------------------//
// txyz.cc
// Geoffrey M. Furnish
// Wed May 13 09:53:44 1998
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the Mesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "../Mesh_XYZ.hh"

#include "nml/Group.hh"

#include "c4/global.hh"

#include <iostream>
using std::cout;
using std::endl;

bool passed = true;

typedef Mesh_XYZ MT;

// The following class exists to test the MT field types with a
// non-trivial type.

class DoubleContainer
{
  public:
    double data;

    DoubleContainer() : data(0.) {}

    DoubleContainer(double _data) : data(_data) {}

    DoubleContainer(const DoubleContainer& dc) : data(dc.data) {}

    DoubleContainer& operator=(const DoubleContainer& dc)
    {
        data = dc.data;
        return *this;
    }
};



// Test the required MT typedefs

void t1()
{
    cout << "t1: beginning.\n";

    typedef MT::size_type size_type;
    typedef MT::FieldConstructor FieldConstructor;

    typedef MT::cctf<DoubleContainer> ccDCf;
    typedef MT::fcdtf<DoubleContainer> fcdDCf;
    typedef MT::nctf<DoubleContainer> ncDCf;
    typedef MT::vctf<DoubleContainer> vcDCf;
    typedef MT::bstf<DoubleContainer> bsDCf;

    typedef MT::ccsf ccsf;
    typedef MT::fcdsf fcdsf;
    typedef MT::ncsf ncsf;
    typedef MT::vcsf vcsf;
    typedef MT::bssf bssf;

    typedef MT::ccif ccif;
    typedef MT::fcdif fcdif;
    typedef MT::ncif ncif;
    typedef MT::vcif vcif;
    typedef MT::bsif bsif;

    typedef MT::ccvsf ccvsf;
    typedef MT::fcdvsf fcdvsf;
    typedef MT::ncvsf ncvsf;
    typedef MT::vcvsf vcvsf;
    typedef MT::bsvsf bsvsf;

    typedef MT::OpAssign OpAssign;
    typedef MT::OpAddAssign OpAddAssign;
    typedef MT::OpSubAssign OpSubAssign;
    typedef MT::OpMultAssign OpMultAssign;
    typedef MT::OpMinAssign OpMinAssign;
    typedef MT::OpMaxAssign OpMaxAssign;

    cout << "t1: end\n";
}



// Test the general MT functions

void t2()
{
    cout << "t2: beginning.\n";

    typedef MT::size_type size_type;
    typedef MT::FieldConstructor FieldConstructor;
    typedef MT::ccsf ccsf;
    typedef MT::fcdsf fcdsf;
    typedef MT::ncsf ncsf;
    typedef MT::vcsf vcsf;
    typedef MT::fcdvsf fcdvsf;

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );
    FieldConstructor FC = spm;

    {
    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

        dsxx::SP<const MT> spm1 = new MT( mdb );
        dsxx::SP<const MT> spm2 = new MT( mdb );
        dsxx::SP<const MT> spm3 = new MT( mdb );

    // Test equivalence relations.

        if (*spm1 == *spm2)
            passed = false;
        if ((*spm1 != *spm2) != !(*spm1 == *spm2))
            passed = false;
        spm2 = spm1;
        if (!(*spm1 == *spm2))
            passed = false;
        if ((*spm1 != *spm2) != !(*spm1 == *spm2))
            passed = false;

    // Invariants

        spm2 = spm1;
        spm3 = spm2;
        if ((spm2 == spm1) && !(*spm2 == *spm1))
            passed = false;
        if (*spm1 != *spm1)
            passed = false;
        if ((*spm1 == *spm2) && !(*spm2 == *spm1))
            passed = false;
        if (((*spm1 == *spm2) && (*spm2 == *spm3)) && !(*spm1 == *spm3))
            passed = false;
    }

    size_type ncells = spm->get_ncells();
    size_type total_ncells = spm->get_total_ncells();
    size_type ncx = spm->get_ncx();
    size_type ncy = spm->get_ncy();
    size_type ncz = spm->get_ncz();
    if (!(ncells <= total_ncells))
        passed = false;
    if (total_ncells != ncx*ncy*ncz)
        passed = false;
    if (!(ncx <= total_ncells))
        passed = false;
    if (!(ncy <= total_ncells))
        passed = false;
    if (!(ncz <= total_ncells))
        passed = false;

    ccsf c(FC);
    fcdsf f(FC);
    fcdvsf fv(FC);
    vcsf v(FC);
    ncsf n(FC);

    spm->get_dx(c);
    spm->get_dy(c);
    spm->get_dz(c);
    spm->get_xloc(c);
    spm->get_yloc(c);
    spm->get_zloc(c);
    spm->get_xloc(f);
    spm->get_yloc(f);
    spm->get_zloc(f);
    spm->get_face_normals(fv);
    spm->get_face_areas(f);
    spm->get_face_lengths(f);
    spm->get_cell_volumes(c);
    spm->get_vertex_volumes(v);
    spm->get_node_volumes(n);

    cout << "t2: end\n";
}



// Test the MT scatters and gathers

void t3()
{
    cout << "t3: beginning.\n";

    typedef MT::FieldConstructor FieldConstructor;
    typedef MT::ccsf ccsf;
    typedef MT::fcdsf fcdsf;
    typedef MT::ncsf ncsf;
    typedef MT::vcsf vcsf;
    typedef MT::bssf bssf;

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );
    FieldConstructor FC = spm;

    ccsf c(FC);
    fcdsf f(FC);
    fcdsf f2(FC);
    vcsf v(FC);
    ncsf n(FC);
    bssf b(FC);

    // Cell centered to face centered scatter

    c = 1.0;
    f = 0.0;
    MT::scatter(f, c, MT::OpAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

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

    // Face centered to cell centered scatter

    f = 1.0;
    c = 0.0;
    MT::scatter(c, f, MT::OpAddAssign());
    for (ccsf::const_iterator citer = c.begin(); citer != c.end();
         ++citer)
      if (*citer != 6.0)
        passed = false;

    f = 3.0;
    c = 1.0;
    MT::scatter(c, f, MT::OpMultAssign());
    for (ccsf::const_iterator citer = c.begin(); citer != c.end();
         ++citer)
      if (*citer != 729.0)
        passed = false;

    // Vertex centered to face centered scatter

    v = 1.0;
    f = 0.0;
    MT::scatter(f, v, MT::OpAddAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 4.0)
        passed = false;

    v = 3.0;
    f = 1.0;
    MT::scatter(f, v, MT::OpMultAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 81.0)
        passed = false;

    // Face centered to vertex centered scatter

    f = 1.0;
    v = 0.0;
    MT::scatter(v, f, MT::OpAddAssign());
    for (vcsf::const_iterator citer = v.begin(); citer != v.end();
         ++citer)
      if (*citer != 3.0)
        passed = false;

    // Vertex centered to node centered scatter

    v = 1.0;
    n = 0.0;
    MT::scatter(n, v, MT::OpAssign());
    for (ncsf::const_iterator citer = n.begin(); citer != n.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    v = 1.0;
    n = 1.0;
    MT::scatter(n, v, MT::OpMultAssign());
    for (ncsf::const_iterator citer = n.begin(); citer != n.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    // Vertex centered to cell centered scatter

    v = 1.0;
    c = 0.0;
    MT::scatter(c, v, MT::OpAddAssign());
    for (ccsf::const_iterator citer = c.begin(); citer != c.end();
         ++citer)
      if (*citer != 8.0)
        passed = false;

    // Cell centered to face centered gather

    c = 1.0;
    f = 0.0;
    MT::gather(f, c, MT::OpAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    c = 5.0;
    f = 3.0;
    MT::gather(f, c, MT::OpMinAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 3.0)
        passed = false;

    c = 5.0;
    f = 6.0;
    MT::gather(f, c, MT::OpMinAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 5.0)
        passed = false;

    c = 5.0;
    f = 3.0;
    MT::gather(f, c, MT::OpMaxAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 5.0)
        passed = false;

    c = 5.0;
    f = 6.0;
    MT::gather(f, c, MT::OpMaxAssign());
    for (fcdsf::const_iterator citer = f.begin(); citer != f.end();
         ++citer)
      if (*citer != 6.0)
        passed = false;

    // Face centered to boundary specified gather

    f = 1.0;
    b = 0.0;
    MT::gather(b, f, MT::OpAssign());
    for (bssf::const_iterator citer = b.begin(); citer != b.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

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

    // Node centered to vertex centered gather

    n = 1.0;
    v = 0.0;
    MT::gather(v, n, MT::OpAssign());
    for (vcsf::const_iterator citer = v.begin(); citer != v.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    n = 1.0;
    v = 0.0;
    MT::gather(v, n, MT::OpAddAssign());
    for (vcsf::const_iterator citer = v.begin(); citer != v.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

    // Cell centered to vertex centered gather

    c = 1.0;
    v = 0.0;
    MT::gather(v, c, MT::OpAssign());
    for (vcsf::const_iterator citer = v.begin(); citer != v.end();
         ++citer)
      if (*citer != 1.0)
        passed = false;

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

    cout << "t3: end\n";
}



// Test the MT global reduction functions

void t4()
{
    cout << "t4: beginning.\n";

    typedef MT::FieldConstructor FieldConstructor;
    typedef MT::ccif ccif;
    typedef MT::fcdif fcdif;
    typedef MT::ncif ncif;
    typedef MT::vcif vcif;
    typedef MT::bsif bsif;

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );
    FieldConstructor FC = spm;

    ccif c(FC);
    fcdif f(FC);
    vcif v(FC);
    ncif n(FC);
    bsif b(FC);

    int value = 3;

    // Cell centered sum

    c = value;
    ccif::value_type csum = MT::sum(c);
    if (csum != value*spm->get_total_ncells())
        passed = false;

    // Face centered sum

    f = value;
    fcdif::value_type fsum = MT::sum(f);
    if (fsum != value*6*spm->get_total_ncells())
        passed = false;

    // Node centered sum

    n = value;
    ncif::value_type nsum = MT::sum(n);
    if (nsum !=
        value*((spm->get_ncx() + 1)*(spm->get_ncy() + 1)*(spm->get_ncz() + 1)))
        passed = false;

    // Vertex centered sum

    v = value;
    vcif::value_type vsum = MT::sum(v);
    if (vsum != value*8*spm->get_total_ncells())
        passed = false;

    // Boundary specified sum

    b = value;
    bsif::value_type bsum = MT::sum(b);
    if (bsum !=
        value*2*((spm->get_ncx()*spm->get_ncy()) +
             (spm->get_ncx()*spm->get_ncz()) +
             (spm->get_ncy()*spm->get_ncz())))
        passed = false;

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

    cout << "t4: end\n";
}



void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    C4::Finalize();
	    return 0;
	}
    }
    
    cout << "Initiating test of the Mesh_XYZ model of the MT concept.\n";

    try {
        t1();
        t2();
        t3();
        t4();
    }
    catch( dsxx::assertion& a )
    {
	cout << "Failed assertion: " << a.what() << endl;
    }

// Print the status of the test.

    cout << endl;
    cout <<     "************************************" << endl;
    if (passed) 
    {
        cout << "**** Mesh_XYZ Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** Mesh_XYZ Self Test: FAILED ****" << endl;
    }
    cout <<     "************************************" << endl;
    cout << endl;

    cout << "Done testing Mesh_XYZ.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of txyz.cc
//---------------------------------------------------------------------------//
