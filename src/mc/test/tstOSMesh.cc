//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstOSMesh.cc
 * \author Thomas M. Evans
 * \date   Thu Dec 20 16:34:09 2001
 * \brief  OS_Mesh test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../OS_Mesh.hh"
#include "../Layout.hh"
#include "../XYCoord_sys.hh"
#include "../XYZCoord_sys.hh"
#include "../OS_Builder.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using rtt_mc::XYCoord_sys;
using rtt_mc::XYZCoord_sys;
using rtt_mc::Layout;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc_test::Parser;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// mesh proxy class
class Mesh_Proxy
{
  private:
    SP<OS_Mesh> mesh;
  public:
    Mesh_Proxy(SP<OS_Mesh> m) : mesh(m) {}
    const OS_Mesh& get_Mesh() const { return *mesh; }
};

//---------------------------------------------------------------------------//

// test the builder

void Builder_2D()
{
    // >>> TEST THE MESH_BUILDER

    // make a builder from parser input 
    SP<Parser> parser(new Parser());
    OS_Builder builder(parser);

    // test cell region data
    {
	vector<int> regions(6,1);
	regions[3] = 2;
	regions[4] = 2;
	regions[5] = 2;

	if (builder.get_num_regions() != 2)   ITFAILS;
	if (builder.get_regions() != regions) ITFAILS;
    }

    // test zone mapping
    {
	vector<vector<int> > zone(2, vector<int>(3));
	zone[0][0] = 1;
	zone[0][1] = 2;
	zone[0][2] = 3;
	zone[1][0] = 4;
	zone[1][1] = 5;
	zone[1][2] = 6;

	if (builder.get_num_zones() != 2)            ITFAILS;
	if (builder.get_cells_in_zone(1) != zone[0]) ITFAILS;
	if (builder.get_cells_in_zone(2) != zone[1]) ITFAILS;
    }

    // build a mesh
    SP<OS_Mesh> mesh = builder.build_Mesh();

    // check defined surface source cells
    {
	vector<vector<int> > ss(2);
	ss[0].resize(3);
	ss[1].resize(2);
	ss[0][0] = 1;
	ss[0][1] = 2;
	ss[0][2] = 3;
	ss[1][0] = 5;
	ss[1][1] = 6;
	
	if (builder.get_defined_surcells() != ss) ITFAILS;

	vector<string> ssp(2);
	ssp[0] = "loy";
	ssp[1] = "hiy";
	
	if (builder.get_ss_pos() != ssp) ITFAILS;
    }

    // check zone mapper
    {
	vector<int> zone_field(2);
	zone_field[0] = 1000;
	zone_field[1] = 1001;
	vector<int> cell_field = builder.zone_cell_mapper(zone_field);

	if (cell_field.size() != mesh->num_cells()) ITFAILS;

	if (cell_field[0] != 1000) ITFAILS;
	if (cell_field[1] != 1000) ITFAILS;
	if (cell_field[2] != 1000) ITFAILS;
	if (cell_field[3] != 1001) ITFAILS;
	if (cell_field[4] != 1001) ITFAILS;
	if (cell_field[5] != 1001) ITFAILS;
    }
}

//---------------------------------------------------------------------------//

// 2D Mesh Tests
void Mesh_2D()
{    
    // make a builder from parser input 
    SP<Parser> parser(new Parser());
    OS_Builder builder_A(parser);
    OS_Builder builder_B(parser);

    // >>> BUILD AND TEST MESHES
   
    // build meshes using the builder
    SP<OS_Mesh> m1  = builder_A.build_Mesh();
    SP<OS_Mesh> m2  = builder_B.build_Mesh();
    SP<OS_Mesh> m1A = builder_A.get_Mesh();
    SP<OS_Mesh> m2B = builder_B.get_Mesh();

    OS_Mesh &m = *m1;

    // the SPs should not be equal
    if (m1 == m2)   ITFAILS; 
    if (*m1 != *m2) ITFAILS;

    // the m1A and m2B meshes should equal their progenitors
    if (m1 != m1A) ITFAILS;
    if (m2 != m2B) ITFAILS;
    
    // check equality of other mesh components
    if (m1->get_Layout() != m2->get_Layout())       ITFAILS;
    if (m1->get_Coord() != m2->get_Coord())         ITFAILS;
    if (m1->get_vertex() != m2->get_vertex())       ITFAILS;
    if (m1->get_cell_pair() != m2->get_cell_pair()) ITFAILS;
    
    // remove our second mesh
    m2 = SP<OS_Mesh>();
    if (m2) ITFAILS;

    // now lets use our proxy classes to ensure mesh equality
    {
	Mesh_Proxy p1(m1);
	Mesh_Proxy p2 = m1;
	if (p1.get_Mesh() != p2.get_Mesh()) ITFAILS;
    }

    // some mesh dimensionality tests
    if (m1->num_cells() != 6)                      ITFAILS;
    if (m1->begin(1) != -1.0 && m1->end(1) != 2.0) ITFAILS;
    if (m1->begin(2) != -1.0 && m1->end(2) != 3.0) ITFAILS;
    
    if (m.pos(1, 2) != 0.5 && m.pos(2, 2) != 0.0) ITFAILS;
    if (m.pos(1, 4) != -.5 && m.pos(2, 4) != 2.0) ITFAILS;
    if (m.pos(1, 6) != 1.5 && m.pos(2, 6) != 2.0) ITFAILS;

    for (int cell = 1; cell <= m.num_cells(); cell++)
    {
	if (m.dim(1, cell) != 1.0) ITFAILS;
	if (m.dim(2, cell) != 2.0) ITFAILS;
    }

    // next cell operations
    if (m.next_cell(1, 2) != 2) ITFAILS;
    if (m.next_cell(1, 1) != 1) ITFAILS;
    if (m.next_cell(2, 4) != 5) ITFAILS;
    if (m.next_cell(3, 3) != 0) ITFAILS;
    if (m.next_cell(4, 1) != 4) ITFAILS;
    if (m.next_cell(5, 2) != 6) ITFAILS;
    if (m.next_cell(6, 1) != 5) ITFAILS;

    // normals
    for (int i = 1; i <= 6; i++)
    {
	vector<double> px(3, 0.0); px[0] = 1.0;
	vector<double> py(3, 0.0); py[1] = 1.0;
	vector<double> mx(3, 0.0); mx[0] = -1.0;
	vector<double> my(3, 0.0); my[1] = -1.0;

	if (m.get_normal(i, 1) != mx)    ITFAILS;
	if (m.get_normal_in(i, 1) != px) ITFAILS;

	if (m.get_normal(i, 2) != px)    ITFAILS;
	if (m.get_normal_in(i, 2) != mx) ITFAILS;

	if (m.get_normal(i, 3) != my)    ITFAILS;
	if (m.get_normal_in(i, 3) != py) ITFAILS;

	if (m.get_normal(i, 4) != py)    ITFAILS;
	if (m.get_normal_in(i, 4) != my) ITFAILS;
    }

    // surface cells
    {
	vector<int> sc(2);
	sc[0] = 1;
	sc[1] = 4;
	if (m.get_surcells("lox") != sc) ITFAILS;
	sc[0] = 1;
	sc[1] = 2;
	sc.push_back(3);
	if (m.get_surcells("loy") != sc) ITFAILS;
    }

    // boundary faces
    for (int i = 1; i <= m.num_cells(); i++)
    {
	if (m.get_bndface("lox", i) != 1) ITFAILS;
	if (m.get_bndface("loy", i) != 3) ITFAILS;
	if (m.get_bndface("hix", i) != 2) ITFAILS;
	if (m.get_bndface("hiy", i) != 4) ITFAILS;
    }

    // vertices
    {
	vector<double> x(4);
	vector<double> y(4);
	x[0] = 0.0;
	x[1] = 1.0;
	x[2] = 1.0;
	x[3] = 0.0;
	y[0] = -1.0;
	y[1] = -1.0;
	y[2] = 1.0;
	y[3] = 1.0;
	
	vector<vector<double> > ref = m.get_vertices(2);
	if (ref[0] != x) ITFAILS;
	if (ref[1] != y) ITFAILS;

	x.resize(2);
	y.resize(2);
	x[0] = 2.0;
	x[1] = 2.0;
	y[0] = 1;
	y[1] = 3;
	ref = m.get_vertices(6, 2);
	if (ref[0] != x) ITFAILS;
	if (ref[1] != y) ITFAILS;
    }

    // face areas
    if (m.face_area(4, 2) != 2.0) ITFAILS;
    if (m.face_area(1, 1) != 2.0) ITFAILS;
    if (m.face_area(5, 4) != 1.0) ITFAILS;
    if (m.face_area(1, 3) != 1.0) ITFAILS;

    // volumes
    for (int i = 1; i <= m.num_cells(); i++)
	if (m.volume(i) != 2.0) ITFAILS;

    // distance-to-boundary calcs
    {
	int face;
	double d;
	vector<double> omega(3, 0.0);
	vector<double> r(2, 0.0);
	
	// case 1
	r[0]     = -.5; 
	r[1]     = -.9;
	omega[0] = .70711;
	omega[1] = .70711;
	omega[2] = 0.0;
	d        = m.get_db(r, omega, 1, face);
	if ((d - .707) >= .001) ITFAILS;
	if (face != 2)          ITFAILS;

	// case 2
	r[0]     = .8; 
	r[1]     = 1.6;
	omega[0] = -.32139;
	omega[1] = .38302;
	omega[2] = .86603;
	d        = m.get_db(r, omega, 5, face);
	if ((d - 2.489) >= .001) ITFAILS;
	if (face != 1)           ITFAILS;

	// case 3
	r[0]     = 1.6; 
	r[1]     = 0.8;
	omega[0] = -.75441;
	omega[1] = -.63302;
	omega[2] = -.17365;
	d        = m.get_db(r, omega, 3, face);
	if ((d - .795) >= .001) ITFAILS;
	if (face != 1)          ITFAILS;

	// case 4
	r[0]     = -.8; 
	r[1]     = 1.2;
	omega[0] = 0.25000;
	omega[1] = -.43301;
	omega[2] = 0.86603;
	d        = m.get_db(r, omega, 4, face);
	if ((d - .462) >= .001) ITFAILS;
	if (face != 3)          ITFAILS;
    }

    // minimum orthogonal distance
    {
	// choose a point in a cell and check the distance
	vector<double> r(2);
	double         ref = 0.0;
	double         d   = 0.0;

	// cell 1
	r[0] = -0.6;
	r[1] = 0.9;
	ref  = 0.1;
	d    = m.get_orthogonal_dist_to_bnd(r, 1);
	
	if (!soft_equiv(d, ref)) ITFAILS;

	r[0] = -0.6;
	r[1] = -0.2;
	ref  = 0.4;
	d    = m.get_orthogonal_dist_to_bnd(r, 1);
	
	if (!soft_equiv(d, ref)) ITFAILS;

	// cell 2
	r[0] = 0.5;
	r[1] = -0.5;
	ref  = 0.5;
	d    = m.get_orthogonal_dist_to_bnd(r, 2);
	
	if (!soft_equiv(d, ref)) ITFAILS;
	
	// cell 5
	r[0] = 0.11;
	r[1] = 2.8;
	ref  = 0.11;
	d    = m.get_orthogonal_dist_to_bnd(r, 5);
	
	if (!soft_equiv(d, ref)) ITFAILS;
    }

    // neighbors
    {
	vector<int> cn(4);

	// cell 1
	cn[0] = 1;
	cn[1] = 2;
	cn[2] = 0;
	cn[3] = 4;
	if (m.get_neighbors(1) != cn) ITFAILS;

	// cell 2
	cn[0] = 1; 
	cn[1] = 3;
	cn[2] = 0;
	cn[3] = 5;
	if (m.get_neighbors(2) != cn) ITFAILS;

	// cell 3
	cn[0] = 2; 
	cn[1] = 0;
	cn[2] = 0;
	cn[3] = 6;
	if (m.get_neighbors(3) != cn) ITFAILS;

	// cell 4
	cn[0] = 4; 
	cn[1] = 5;
	cn[2] = 1;
	cn[3] = 0;
	if (m.get_neighbors(4) != cn) ITFAILS;

	// cell 5
	cn[0] = 4; 
	cn[1] = 6;
	cn[2] = 2;
	cn[3] = 0;
	if (m.get_neighbors(5) != cn) ITFAILS;

	// cell 6
	cn[0] = 5; 
	cn[1] = 0;
	cn[2] = 3;
	cn[3] = 0;
	if (m.get_neighbors(6) != cn) ITFAILS;
    }

    // cell types
    vector<int> cell_types(m.num_cells(), 5);
    if (cell_types != m.get_cell_types()) ITFAILS;

    // point coordinates
    {
	vector<vector<double> > ptc(12, vector<double>(2));
	ptc[0][0]  = -1;
	ptc[1][0]  = 0;
	ptc[2][0]  = 1;
	ptc[3][0]  = 2;
	ptc[4][0]  = -1;
	ptc[5][0]  = 0;
	ptc[6][0]  = 1;
	ptc[7][0]  = 2;
	ptc[8][0]  = -1;
	ptc[9][0]  = 0;
	ptc[10][0] = 1;
	ptc[11][0] = 2;
	
	ptc[0][1]  = -1;
	ptc[1][1]  = -1;
	ptc[2][1]  = -1;
	ptc[3][1]  = -1;
	ptc[4][1]  = 1;
	ptc[5][1]  = 1;
	ptc[6][1]  = 1;
	ptc[7][1]  = 1;
	ptc[8][1]  = 3;
	ptc[9][1]  = 3;
	ptc[10][1] = 3;
	ptc[11][1] = 3;
	
	if (ptc != m.get_point_coord()) ITFAILS;
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is a serial test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// 2D Mesh tests
	Builder_2D();
	Mesh_2D();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstOSMesh, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstOSMesh Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing tstOSMesh on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstOSMesh.cc
//---------------------------------------------------------------------------//
