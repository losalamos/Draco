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
using namespace std;

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
    dsxx::SP<MT> spm = new MT( mdb );
    FieldConstructor FC = spm;

    {
        dsxx::SP<MT> spm2 = new MT( mdb );
        if (spm.bp() == spm2.bp())
            passed = false;
        if ((spm.bp() != spm2.bp()) != !(spm.bp() == spm2.bp()))
            passed = false;
        spm2 = spm;
        if (!(spm.bp() == spm2.bp()))
            passed = false;
        if ((spm.bp() != spm2.bp()) != !(spm.bp() == spm2.bp()))
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
