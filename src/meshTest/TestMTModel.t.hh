//----------------------------------*-C++-*----------------------------------//
// TestMTModel.t.hh
// Randy M. Roberts
// Fri Aug 20 09:11:42 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Test the MT Model for its compliance with the MT Concept.
//---------------------------------------------------------------------------//

#include "TestMTModel.hh"
#include "DoubleContainer.hh"

#include <iostream>

namespace rtt_meshTest
{

template<class MTFactory>
void TestMTModel<MTFactory>::run()
{
    // Run the tests in this test class.

    os_m << "Begin Running....... TestMTModel tests." << std::endl;

    passed_m = true;
    
    t1();

    t2();

    os_m << "Completed Running... TestMTModel tests." << std::endl;
}


// Test the required MT typedefs

template<class MTFactory>
void TestMTModel<MTFactory>::t1()
{
    os_m << "t1: beginning.\n";

    typedef typename MT::size_type size_type;
    typedef typename MT::FieldConstructor FieldConstructor;

    typedef typename MT::cctf<DoubleContainer> ccDCf;
    typedef typename MT::fcdtf<DoubleContainer> fcdDCf;
    typedef typename MT::nctf<DoubleContainer> ncDCf;
    typedef typename MT::vctf<DoubleContainer> vcDCf;
    typedef typename MT::bstf<DoubleContainer> bsDCf;

    typedef typename MT::ccsf ccsf;
    typedef typename MT::fcdsf fcdsf;
    typedef typename MT::ncsf ncsf;
    typedef typename MT::vcsf vcsf;
    typedef typename MT::bssf bssf;

    typedef typename MT::ccif ccif;
    typedef typename MT::fcdif fcdif;
    typedef typename MT::ncif ncif;
    typedef typename MT::vcif vcif;
    typedef typename MT::bsif bsif;

    typedef typename MT::ccvsf ccvsf;
    typedef typename MT::fcdvsf fcdvsf;
    typedef typename MT::ncvsf ncvsf;
    typedef typename MT::vcvsf vcvsf;
    typedef typename MT::bsvsf bsvsf;

    typedef typename MT::OpAssign OpAssign;
    typedef typename MT::OpAddAssign OpAddAssign;
    typedef typename MT::OpSubAssign OpSubAssign;
    typedef typename MT::OpMultAssign OpMultAssign;
    typedef typename MT::OpMinAssign OpMinAssign;
    typedef typename MT::OpMaxAssign OpMaxAssign;

    os_m << "t1: end\n";
}



// Test the general MT functions

template<class MTFactory>
void TestMTModel<MTFactory>::t2()
{
    os_m << "t2: beginning.\n";

    typedef typename MT::size_type size_type;
    typedef typename MT::FieldConstructor FieldConstructor;
    typedef typename MT::ccsf ccsf;
    typedef typename MT::fcdsf fcdsf;
    typedef typename MT::ncsf ncsf;
    typedef typename MT::vcsf vcsf;
    typedef typename MT::fcdvsf fcdvsf;

    // need a temporary passed variable for each subset of tests.

    bool passed = true;

    {
	// Create 3 meshes to test equivalence relations
	// and invariants.

	MTFactoryProduct meshProduct1 = meshFactory_m.create();
	MTFactoryProduct meshProduct2 = meshFactory_m.create();
	MTFactoryProduct meshProduct3 = meshFactory_m.create();

	MT &mesh1 = meshProduct1.mesh();
	MT &mesh2 = meshProduct2.mesh();
	MT &mesh3 = meshProduct3.mesh();

	// Insure that the MTFactory produced "distinct" meshes
	// for each call of create.

	passed = true;
	
        if (mesh1 == mesh2)
            passed = false;
        if (mesh1 == mesh3)
            passed = false;
        if (mesh3 == mesh2)
            passed = false;

	if (!passed)
	    os_m << "t2() failed testing MTFactory distinct products."
		 << std::endl;

	// update the object state.
	
	passed_m = passed && passed_m;
	passed = true;

	// Test equivalence relations.

        if ((mesh1 != mesh2) != !(mesh1 == mesh2))
            passed = false;

        if (!(mesh1 == mesh1))
            passed = false;
        if ((mesh1 != mesh1) != !(mesh1 == mesh1))
            passed = false;

	if (!passed)
	    os_m << "t2() failed testing MT equivalence relations."
		 << std::endl;
	
	// update the object state.
	
	passed_m = passed && passed_m;
	passed = true;
	
	// Invariants
	// Remember... !(A => B) <=> (A && !B)

	// Identity
        if ((&mesh1 == &mesh2) && !(mesh1 == mesh2))
            passed = false;
	// Reflexivity
        if (!(mesh1 == mesh1))
            passed = false;
	// Symmetry
        if ((mesh1 == mesh2) && !(mesh2 == mesh1))
            passed = false;
	// Transitivity
        if (((mesh1 == mesh2) && (mesh2 == mesh3)) && !(mesh1 == mesh3))
            passed = false;

	if (!passed)
	    os_m << "t2() failed testing MT Invariants."
		 << std::endl;
	
	// update the object state.
	
	passed_m = passed && passed_m;
	passed = true;
    }

    // Get the mesh and field constructor from the mesh factory.
    
    MTFactoryProduct meshProduct = meshFactory_m.create();
    MT &mesh = meshProduct.mesh();
    FieldConstructor &fCtor = meshProduct.fieldConstructor();
    
    size_type ncells = mesh.get_ncells();
    size_type total_ncells = mesh.get_total_ncells();
    size_type ncx = mesh.get_ncx();
    size_type ncy = mesh.get_ncy();
    size_type ncz = mesh.get_ncz();

    passed = true;
    
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

    if (!passed)
	os_m << "t2() failed testing MT mesh size relations."
	     << std::endl;
	
    // update the object state.
	
    passed_m = passed && passed_m;
    passed = true;

    ccsf c(fCtor);
    fcdsf f(fCtor);
    fcdvsf fv(fCtor);
    vcsf v(fCtor);
    ncsf n(fCtor);

    mesh.get_dx(c);
    mesh.get_dy(c);
    mesh.get_dz(c);
    mesh.get_xloc(c);
    mesh.get_yloc(c);
    mesh.get_zloc(c);
    mesh.get_xloc(f);
    mesh.get_yloc(f);
    mesh.get_zloc(f);
    mesh.get_face_normals(fv);
    mesh.get_face_areas(f);
    mesh.get_face_lengths(f);
    mesh.get_cell_volumes(c);
    mesh.get_vertex_volumes(v);
    mesh.get_node_volumes(n);

    os_m << "t2: end\n";
}

} // end namespace rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of TestMTModel.t.hh
//---------------------------------------------------------------------------//
