//----------------------------------*-C++-*----------------------------------//
// OS_Mesh.cc
// Thomas M. Evans
// Tue Feb  3 16:50:13 1998
//---------------------------------------------------------------------------//
// @> OS_Mesh class implementation file
//---------------------------------------------------------------------------//

#include "imctest/OS_Mesh.hh"
#include "imctest/Constants.hh"
#include <iomanip>

IMCSPACE

using std::sort;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// default constructor
OS_Mesh::OS_Mesh(SP<Coord_sys> coord_, Layout &layout_, CCVF_a &vertex_, 
		 CCVF_i &cell_pair_) 
    : coord(coord_), layout(layout_), vertex(vertex_),
      cell_pair(cell_pair_), sur(coord->get_dim()) 
{
  // assertions to verify size of mesh and existence of a Layout and
  // Coord_sys  
    Check (coord);
	
  // variable initialization
    int ncells = num_cells();
    int dimension = coord->get_dim();
    
  // dimension assertions
    Check (dimension == vertex.size());
    Check (dimension == sur.size());
    
  // mesh size assertions
    Check (ncells == cell_pair.size());
      
  // calculate surface array
    calc_surface();
}

//---------------------------------------------------------------------------//
// private member functions
//---------------------------------------------------------------------------//
void OS_Mesh::calc_surface()
{
  // calculate an array of the dimensional surfaces which make up the OS_Mesh

  // initialize mesh_size for assertion at end of function
    int mesh_size = 1;

  // loop to calculate surface array
    for (int d = 0; d < coord->get_dim(); d++)
    {
      // define an array for dim which is sorted in ascending order
	vector<double> sorted = vertex[d];
	sort(sorted.begin(), sorted.end());

      // loop over sorted array, appending new surfaces onto sur array
	sur[d].push_back(sorted[0]);
	for (int i = 1; i < sorted.size(); i++)
	    if (sorted[i] > sorted[i-1])
		sur[d].push_back(sorted[i]);

      // calculate mesh_size by dimension
	mesh_size *= sur[d].size() - 1;
    }

  // assert mesh size
    Check (num_cells() == mesh_size);
}

//---------------------------------------------------------------------------//
// member functions
//---------------------------------------------------------------------------//
// do binary search on a cell

int OS_Mesh::get_cell(const vector<double> &r) const
{
  // used variables
    int dim         = coord->get_dim();
    int return_cell = 1;
    int subcells    = 1;
    
  // binary search of cells

    for (int i = 0; i < dim; i++)
    {
	int low_index  = 0;
	int high_index = sur[i].size() - 1;
	int index;
	while ((high_index - low_index) != 1)
	{
	    index = (high_index + low_index) / 2;
	    if (r[i] < sur[i][index])
		high_index = index;
	    else
		low_index  = index;
	}
	return_cell += subcells * (high_index - 1);
      // number of cells per dimension equals the number of surfaces along
      // that dimension minus one
	subcells    *= sur[i].size() - 1;
    }  
    
  // return cell index
    return return_cell;
}

//---------------------------------------------------------------------------//
// calculate the distance to boundary

double OS_Mesh::get_db(const vector<double> &r, const vector<double> &omega,
		       int cell, int &face) const
{
    using std::vector;
    using std::min_element;
    using Global::huge;
    
  // calculate distance to the vec(r) boundaries

  // boundary distances over each coordinate direction
    vector<double> dim_dist_boundary(coord->get_dim(), 0.0);
    
  // loop to get the distances to boundary in each coordinate direction
    for (int i = 0; i < coord->get_dim(); i++)
    {
      // define absolute dimension index
	int d = i + 1;

      // find the distances to boundary along each dimension
	if (omega[i] == 0.0)
	    dim_dist_boundary[i] = Global::huge;
	else if (omega[i] > 0.0)
	    dim_dist_boundary[i] = (max(d, cell) - r[i]) / omega[i];
	else
	    dim_dist_boundary[i] = (min(d, cell) - r[i]) / omega[i];
    }

  // calculate the distance to boundary
    vector<double>::iterator itor = min_element(dim_dist_boundary.begin(),
						dim_dist_boundary.end());
    double dist_boundary = *itor;

  // calculate the face that the boundary is on
    int index = itor - dim_dist_boundary.begin();
    if (omega[index] < 0.0)
	face = 1 + 2 * index;
    else
	face = 2 + 2 * index;

  // return the distance-to-boundary
    return dist_boundary;
}

//---------------------------------------------------------------------------//
// return the face number for a given cell boundary

int OS_Mesh::get_bndface(string boundary, int cell) const
{
  // return the face number for boundary on cell

  // return value
    int face;

    if (boundary == "lox")
	face = 1;
    if (boundary == "hix")
	face = 2;
    if (boundary == "loy")
	face = 3;
    if (boundary == "hiy")
	face = 4;
    if (boundary == "loz")
	face = 5;
    if (boundary == "hiz")
	face = 6;

  // return the face
    return face;
}

//---------------------------------------------------------------------------//
// return a list of cells along a specified boundary

vector<int> OS_Mesh::get_surcells(string boundary) const
{
  // return a list of cells along the specified boundary

  // make return vector containing a list of cells along specified boundary
    vector<int> return_list;

    int num_xcells = 1;
    int num_ycells = 1;
    int num_zcells = 1;

  // set dimensionality variables
    if (coord->get_dim() == 1)
    {
	num_xcells = sur[0].size() - 1;
	num_ycells = 1;
	num_zcells = 1;
    }
    else if (coord->get_dim() == 2)
    {
	num_xcells = sur[0].size() - 1;
	num_ycells = sur[1].size() - 1;
	num_zcells = 1;
    }
    else if (coord->get_dim() == 3)
    {
	num_xcells = sur[0].size() - 1;
	num_ycells = sur[1].size() - 1;
	num_zcells = sur[2].size() - 1;
    }

  // calculate the cells along the boundary
    if (boundary == "lox")
    {
	for (int k = 1; k <= num_zcells; k++)
	    for (int j = 1; j <= num_ycells; j++)
	    {
		int bcell = 1 + num_xcells * (j - 1) + 
		    num_xcells * num_ycells * (k - 1);
		return_list.push_back(bcell);
	    }
    }
    if (boundary == "hix")
    {
        for (int k = 1; k <= num_zcells; k++)
	    for (int j = 1; j <= num_ycells; j++)
	    {
		int bcell = 1 + (num_xcells - 1) + num_xcells * (j - 1) +
		    num_xcells * num_ycells * (k - 1);
		return_list.push_back(bcell);
	    }
    }
    if (boundary == "loy")
    {
	for (int k = 1; k <= num_zcells; k++)
	    for (int i = 1; i <= num_xcells; i++)
	    {
		int bcell = 1 + (i - 1) + num_xcells * num_ycells * (k - 1);
		return_list.push_back(bcell);
	    }	
    }
    if (boundary == "hiy")
    {
       	for (int k = 1; k <= num_zcells; k++)
	    for (int i = 1; i <= num_xcells; i++)
	    {
		int bcell = 1 + (i - 1) + num_xcells * (num_ycells - 1) +
		    num_xcells * num_ycells * (k - 1);
		return_list.push_back(bcell);
	    }
    }
    if (boundary == "loz")
    {
	for (int j = 1; j <= num_ycells; j++)
	    for (int i = 1; i <= num_xcells; i++)
	    {
		int bcell = 1 + (i - 1) + num_xcells * (j - 1);
		return_list.push_back(bcell);
	    }
    }
    if (boundary == "hiz")
    {
	for (int j = 1; j <= num_ycells; j++)
	    for (int i = 1; i <= num_xcells; i++)
	    {
		int bcell = 1 + (i - 1) + num_xcells * (j - 1) + 
		    num_xcells * num_ycells * (num_zcells - 1);
		return_list.push_back(bcell);
	    }
    }

  // return vector
    return return_list;
}

//---------------------------------------------------------------------------//
// Overloaded operators
//---------------------------------------------------------------------------//
// overloaded == for design-by-contract

bool OS_Mesh::operator==(const OS_Mesh &rhs) const
{
  // check to see that we have the same coordinate systems
    if (coord != rhs.coord)
	return false;

  // check to see that the Layouts are equal
    if (layout != rhs.layout)
	return false;

  // check the vertices
    if (vertex != rhs.vertex)
	return false;
    if (cell_pair != rhs.cell_pair)
	return false;

  // if we haven't returned, then the two meshes must be equal
    return true;
}

//---------------------------------------------------------------------------//
// public diagnostic member functions
//---------------------------------------------------------------------------//

void OS_Mesh::print(ostream &output, int cell) const
{
    using std::cout;
    using std::endl;
    using std::setw;

  // print out content info for 1 cell
    output << "+++++++++++++++" << endl;
    output << "---------------" << endl;
    output << "Cell : "         << cell << endl;
    output << "---------------" << endl;
    output << "Dimensions "     << endl;
    output << "---------------" << endl;
    if (coord->get_dim() == 2)
    {
	output << " x  : " << pos(1, cell) << endl;
	output << " y  : " << pos(2, cell) << endl;
    	output << " dx : " << dim(1, cell) << endl;
	output << " dy : " << dim(2, cell) << endl;
    }
    else
    {
	output << " x  : " << pos(1, cell) << endl;
	output << " y  : " << pos(2, cell) << endl;
	output << " z  : " << pos(3, cell) << endl;
    	output << " dx : " << dim(1, cell) << endl;
	output << " dy : " << dim(2, cell) << endl;
	output << " dz : " << dim(3, cell) << endl;
    }	
    output << "---------------" << endl;
    output << "Layout "         << endl;
    output << "---------------" << endl;
    layout.print(output, cell);
    output << "+++++++++++++++" << endl;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of OS_Mesh.cc
//---------------------------------------------------------------------------//
