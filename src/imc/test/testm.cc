//----------------------------------*-C++-*----------------------------------//
// testm.cc
// Thomas M. Evans
// Wed Feb  4 17:04:54 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "XYCoord_sys.hh"
#include "XYZCoord_sys.hh"
#include "Layout.hh"
#include "Global.hh"
#include "Random.hh"
#include "OS_Mesh.hh"
#include <vector>
#include <string>
#include "SP.hh"

IMCSPACE

using std::string;

class OSBuild
{
private:
    void Assign2D(vector<int> &dim_size, vector<string> &bnd_cond,
		  Layout &layout)
    {
      // 2D map of Mesh
	int num_xcells = dim_size[0];
	int num_ycells = dim_size[1];
	
      // loop over num_cells and assign cell across faces
	for (int cell = 1; cell <= layout.getNum_cell(); cell++)
	{
	    layout(cell,1) = cell - 1;
	    layout(cell,2) = cell + 1;
	    layout(cell,3) = cell - num_xcells;
	    layout(cell,4) = cell + num_xcells;
	}

      // take care of boundary conditions

      // low x boundary
	{
	    int bcell = 1;
	    for (int i = 1; i <= num_ycells; i++)
	    {
		if (bnd_cond[0] == "vacuum")
		    layout(bcell,1) = 0;
		else if (bnd_cond[0] == "reflect")
		    layout(bcell,1) = bcell;
		bcell += num_xcells;
	    }
	}
      // high x boundary
	{
	    int bcell = num_xcells;
	    for (int i = 1; i <= num_ycells; i++)
	    {
		if (bnd_cond[1] == "vacuum")
		    layout(bcell,2) = 0;
		else if (bnd_cond[1] == "reflect")
		    layout(bcell,2) = bcell;
		bcell += num_xcells;
	    }
	}	
      // low y boundary
	{
	    int bcell = 1;
	    for (int i = 1; i <= num_xcells; i++)
	    {
		if (bnd_cond[2] == "vacuum")
		    layout(bcell,3) = 0;
		else if (bnd_cond[2] == "reflect")
		    layout(bcell,3) = bcell;
		bcell += 1;
	    }
	}
      // high y boundary
	{
	    int bcell = num_xcells * (num_ycells-1) + 1;
	    for (int i = 1; i <= num_xcells; i++)
	    {
		if (bnd_cond[3] == "vacuum")
		    layout(bcell,4) = 0;
		else if (bnd_cond[3] == "reflect")
		    layout(bcell,4) = bcell;
		bcell += 1;
	    }
	}
    }
    void Assign3D(vector<int> &dim_size, vector<string> &bnd_cond,
		  Layout &layout)
    {}
    SP<Coord_sys> buildXY()
    {
	SP<XYCoord_sys> coord = new XYCoord_sys;
	SP<Coord_sys> base_coord = coord;
	return base_coord;
    }
    SP<Coord_sys> buildXYZ()
    {
	SP<XYZCoord_sys> coord = new XYZCoord_sys;
	SP<Coord_sys> base_coord = coord;
	return base_coord;
    }	   
    SP<OS_Mesh> build2DMesh(OS_Mesh::CCVF_a &fdata, SP<Coord_sys> coord,
			    const Layout &layout)
    {
      // variable declarations
	int num_cells  = layout.getNum_cell();
	int num_xsur   = fdata[0].size();
	int num_ysur   = fdata[1].size();
	int num_xcells = num_xsur - 1;
	int num_ycells = num_ysur - 1;
	int cell       = 0;

      // initialization variables for Mesh
	OS_Mesh::CCVF_a pos(2);
	OS_Mesh::CCVF_a dim(2);
	OS_Mesh::CCF_i  index(num_cells);

      // size pos and loc arrays
	for (int d = 1; d <= 2; d++)
	{
	    pos[d-1].resize(num_cells);
	    dim[d-1].resize(num_cells);
	}

      // set pos arrays
	for (int i = 1; i <= num_xcells; i++)
	    for (int j = 1; j <= num_ycells; j++)
	    {
		cell           = 1 + (i-1) + num_xcells * (j-1);
		pos[0][cell-1] = (fdata[0][i-1] + fdata[0][i]) / 2.0;
		pos[1][cell-1] = (fdata[1][j-1] + fdata[1][j]) / 2.0;
		dim[0][cell-1] = fdata[0][i] - fdata[0][i-1];
		dim[1][cell-1] = fdata[1][j] - fdata[1][j-1];
		index[cell-1]  = cell;
	    }
	
      // return mesh to builder
	SP<OS_Mesh> mesh_return = new OS_Mesh(coord, layout, pos, dim, 
					      fdata, index);
	return mesh_return;
    }
    SP<OS_Mesh> build3DMesh(OS_Mesh::CCVF_a &fdata, SP<Coord_sys> coord,
			    const Layout &layout)
    {
	SP<OS_Mesh> mesh_return;
	mesh_return = 0;
	return mesh_return;
    }
public:
    SP<Coord_sys> buildCoord(string &system)
    {
	SP<Coord_sys> coord;
	if (system == "xy" || system == "XY")
	    coord = buildXY();
	else if (system == "xyz" || system == "XYZ")
	    coord = buildXYZ();
	return coord;
    }
    SP<Layout> buildLayout(int size, vector<int> &dsize, 
			   vector<string> &bnd_cond, Coord_sys &coord)
    {
      // set size of new Layout
	SP<Layout> layout = new Layout(size);
	for (int i = 1; i <= size; i++)
	    layout->setSize(i, coord.getDim()*2);

      // assign cells and faces to Layout
	if (coord.getDim() == 2)
	    Assign2D(dsize, bnd_cond, *layout);
	else
	    Assign3D(dsize, bnd_cond, *layout);
	
	return layout;
    }
    SP<OS_Mesh> buildMesh(OS_Mesh::CCVF_a &face_data, SP<Coord_sys> coord,
			  const Layout &layout)
    {
	SP<OS_Mesh> mesh;
	int dim = coord->getDim();
	
      // call private member build functions
	if (dim == 2)
	    mesh = build2DMesh(face_data, coord, layout);
	else if (dim == 3)
	    mesh = build3DMesh(face_data, coord, layout);
	return mesh;
    }	    
};
	
void f(Random &random, const OS_Mesh &mesh)
{
    using std::cout;
    using std::endl;

    vector<double> r(2);
    r[0] = mesh.begin(1) + random.ran() * mesh.end(1);
    r[1] = mesh.begin(2) + random.ran() * mesh.end(2);

    int cell = mesh.getCell(r);

    cout << endl;
    cout << "Starting Cell : " << cell << endl;
    cout << "X             : " << r[0] << endl;
    cout << "Y             : " << r[1] << endl;
}

CSPACE

main()
{
    using namespace IMC;
    using namespace std;
    using Global::operator<<;
    
    SP<OS_Mesh> mesh;

    {
	SP<Coord_sys> coord;
	SP<Layout> layout;
	OSBuild build;

      // PARSER STUFF

	string system;
	cout << "Enter Coord_sys:" << endl;
	cin >> system;

	int mesh_size = 6;
	vector<int> dim_size(2);
	dim_size[0] = 3;
	dim_size[1] = 2;

	vector<string> bnd_cond(4);
	bnd_cond[0] = "vacuum";
	bnd_cond[1] = "reflect";
	bnd_cond[2] = "reflect";
	bnd_cond[3] = "reflect";

	OS_Mesh::CCVF_a edges(2);
	edges[0].resize(dim_size[0]+1);
	edges[1].resize(dim_size[1]+1);
	edges[0][0] = 0.0;
	edges[0][1] = 2.0;
	edges[0][2] = 4.0;
	edges[0][3] = 6.0;
	edges[1][0] = 0.0;
	edges[1][1] = 1.0;
	edges[1][2] = 2.0;

      // END PARSER STUFF

	coord  = build.buildCoord(system);
	layout = build.buildLayout(mesh_size, dim_size, bnd_cond, *coord);
	mesh = build.buildMesh(edges, coord, *layout);
    }

  // cell tests
    
    Random random(-142573);

    for (int i = 1; i <= 10; i++)
	f(random, *mesh);

}

//---------------------------------------------------------------------------//
//                              end of testm.cc
//---------------------------------------------------------------------------//
