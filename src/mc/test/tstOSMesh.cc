//----------------------------------*-C++-*----------------------------------//
// tstOSMesh.cc
// Thomas M. Evans
// Wed Apr 21 18:50:01 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Test of OS_Mesh and OS_Builder
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../OS_Mesh.hh"
#include "../Layout.hh"
#include "../XYCoord_sys.hh"
#include "../XYZCoord_sys.hh"
#include "../OS_Builder.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

using rtt_mc::XYCoord_sys;
using rtt_mc::XYZCoord_sys;
using rtt_mc::Layout;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc_test::MC_Interface;
using dsxx::SP;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

// mesh proxy class
class Mesh_Proxy
{
  private:
    SP<OS_Mesh> mesh;
  public:
    Mesh_Proxy(SP<OS_Mesh> m) : mesh(m) {}
    const OS_Mesh& get_Mesh() const { return *mesh; }
};

// 2D Mesh Tests
void Test_2D()
{
    // build 2 identical meshes 
    SP<MC_Interface> interface(new MC_Interface());
    OS_Builder builder(interface);
    
    SP<OS_Mesh> m1 = builder.build_Mesh();
    SP<OS_Mesh> m2 = builder.build_Mesh();
    OS_Mesh &m = *m1;

    // the SPs should not be equal, also because of SPs to the Coord_sys
    // inside of mesh, the meshes are not identically equal
    if (m1 == m2)   ITFAILS; 
    if (*m1 == *m2) ITFAILS;
    
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
	x[2] = 0.0;
	x[3] = 1.0;
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
}

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
	    cout << argv[0] << ": version " << rtt_mc::release() << endl;
	    C4::Finalize();
	    return 0;
	}

    // 2D Mesh tests
    Test_2D();

    // status of test
    cout << endl;
    cout <<     "***********************************" << endl;
    if (passed) 
    {
        cout << "**** OS_Mesh Self Test: PASSED ****" << endl;
    }
    cout <<     "***********************************" << endl;
    cout << endl;

    cout << "Done testing OS_Mesh." << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                              end of tstOSMesh.cc
//---------------------------------------------------------------------------//
