//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Mesh.cc
// B.T. Adams (bta@lanl.gov)
// 18 May 99
//---------------------------------------------------------------------------//
// @> CAR_CU_Mesh class implementation file (developed from OS_Mesh.cc)
//---------------------------------------------------------------------------//

#include "CAR_CU_Mesh.hh"
#include "Constants.hh"
#include <iomanip>

namespace rtt_mc 
{

// std components
using std::sort;
using std::cout;
using std::endl;
using std::setw;
using std::ios;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// default constructor
CAR_CU_Mesh::CAR_CU_Mesh(SP<Coord_sys> coord_, Layout & layout_, 
   NCVF_d & vertex_, CCVF_i & cell_pair_,  CCSF_i & generation_,
   bool submesh_) : coord(coord_), layout(layout_), vertex(vertex_), 
   cell_pair(cell_pair_), generation(generation_), sur(coord->get_dim()), 
   submesh(submesh_)
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
    
  // calculate surface array
    if (!submesh)
	calc_surface();
}

//---------------------------------------------------------------------------//
// private member functions
//---------------------------------------------------------------------------//
bool CAR_CU_Mesh::compReal(const double & low_val, const double & high_val)
    const
{
    // require agreement to six significant figures for equality (set by
    // exponent of EPSILON).
    const double EPSILON = 1.0e-06;
    double epsilon;
    if (low_val != 0 && high_val != 0)
        epsilon = EPSILON * ((fabs(low_val) + fabs(high_val))/2.);
    else
        epsilon = EPSILON;

    return (fabs(high_val - low_val) > epsilon);
}

void CAR_CU_Mesh::calc_surface()
{
    // calculate an array of the dimensional surfaces which make up the 
    // CAR_CU_Mesh

    // loop to calculate surface array
    for (int d = 0; d < coord->get_dim(); d++)
    {
	// define an array for dim which is sorted in ascending order
	vector<double> sorted = vertex[d];
	sort(sorted.begin(), sorted.end());

	// loop over sorted array, appending new surfaces onto sur array
	// using the compReal function to check for machine level equality
	sur[d].push_back(sorted[0]);
	for (int i = 1; i < sorted.size(); i++)
	    if (compReal(sorted[i-1], sorted[i]))
		sur[d].push_back(sorted[i]);
    }
}

//---------------------------------------------------------------------------//
// member functions
//---------------------------------------------------------------------------//
// do binary search on a cell

int CAR_CU_Mesh::get_cell(const vector<double> &r) const
{
    Require (!submesh);

  // used variables
    int dim         = coord->get_dim();
    int return_cell;
    int index;

    // binary search of cell vertexes - vertex array is sorted with x varying 
    // fastest, followed by y, and then z. find the first vertex "below" this
    // point in space.
    int high_index = vertex[0].size() - 1;
    for (int i = dim; i >= 0; i--)
    {
        int low_index  = 0;
	while ((high_index - low_index) > 1)
	{
	    index = (high_index + low_index) / 2;
	    if (r[i] < vertex[i][index] && compReal(r[i],vertex[i][index]))
		high_index = index;
	    else
		low_index  = index;
	}
	high_index = index;
    }
    // have the vertex index - find the cell that has this vertex at its
    // negative x,y,z corner (first vertex), again using a binary search.

    int low_index  = 0;
    high_index = num_cells() - 1;
    while ((high_index - low_index) > 1)
    {
        return_cell = (high_index + low_index) / 2;
	if (cell_pair[return_cell][0] - 1 < index)
	    high_index = return_cell;
	else
	    low_index  = return_cell;
    }
  // return cell index
    return return_cell;
}

//---------------------------------------------------------------------------//
// calculate the distance to boundary

double CAR_CU_Mesh::get_db(const vector<double> &r, 
   const vector<double> &omega, int cell, int &face) const
{
    using std::vector;
    using std::min_element;
    using global::huge;
    
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
	    dim_dist_boundary[i] = global::huge;
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
	face = 3 - index;
    else
	face = 4 + index;

  // return the distance-to-boundary
    return dist_boundary;
}

//---------------------------------------------------------------------------//
// return the face number for a given cell boundary

int CAR_CU_Mesh::get_bndface(string boundary, int cell) const
{
  // return the face number for boundary on cell

  // return value
    int face;

    if (boundary == "loz")
	face = 1;
    else if (boundary == "loy")
	face = 2;
    else if (boundary == "lox")
	face = 3;
    else if (boundary == "hix")
	face = 4;
    else if (boundary == "hiy")
	face = 5;
    else if (boundary == "hiz")
	face = 6;
    else
        Insist(0, "Illegal cell boundary face!");

  // return the face
    return face;
}

//---------------------------------------------------------------------------//
// return a list of cells along a specified boundary

vector<int> CAR_CU_Mesh::get_surcells(string boundary) const
{
    Require (!submesh);

    // return a list of cells along the specified boundary

    // make return vector containing a list of cells along specified boundary
    vector<int> return_list;

    // index to the end of the boundary string. Allowed values for the 
    // boundary string are loz, loy, lox, hix, hiy, hiz.
    int end_ind = boundary.size() - 1;
    // the vertex coordinate dimension index associated with the boundary 
    // surface (-z = 2, -y = 1, -x = 0, +x = 0, +y = 1, +z = 2);
    int dim_ind = 0;
    // the value of the fixed coordinate along the face
    double value;
    // starting cell number index (zero for negative faces, ncells - 1 for
    // the positive faces)
    int cell_ind = 0;
    // node checked in each cell to determine if it lies on this surface
    int node_ind = 0;

    // reset the starting cell number index and node index for positive faces.
    if (boundary[0] == 'h')
    {
        cell_ind = num_cells() - 1;
	if (boundary[end_ind] == 'x')
	    node_ind = 1;
	else if (boundary[end_ind] == 'y')
	    node_ind = 2;
	else if (boundary[end_ind] == 'z')
	    node_ind = 4;
    }

    // reset the dimension index for y and z faces
    if (boundary[end_ind] == 'y')
        dim_ind = 1;
    else if (boundary[end_ind] == 'z')
        dim_ind = 2;

    // calculate the cells along the boundary. Have to do a linear search for
    // cells on the -y, -x, x, and y boundaries for the time being. In the
    // future it might be smart to create static vectors that mark the changes
    // in z and y to speed up the process. The need for such depends on how
    // often this routine is called.
    if (boundary[end_ind] != 'z')
    {
	value = vertex[dim_ind][cell_pair[cell_ind][node_ind] - 1];
	for (int k = 0; k < num_cells(); k++)
	{
	    if (!compReal(vertex[dim_ind][cell_pair[k][node_ind] - 1],value))
		return_list.push_back(k+1);
	}
    }
    // only have to search the beginning (-z) and end (+z) of the last 
    // column of the vertex vector since they are in ascending order.
    else
    {
        if (boundary[0] == 'l')
	{
	    value = vertex[dim_ind][cell_pair[cell_ind][node_ind] - 1];
	    while (!compReal(vertex[dim_ind][cell_pair[cell_ind][node_ind] -1],
		   value))
	    {
	        return_list.push_back(cell_ind+1);
		++cell_ind;
	    }
	}
	else
	{
	    value = vertex[dim_ind][cell_pair[cell_ind][node_ind] - 1];
	    while (!compReal(vertex[dim_ind][cell_pair[cell_ind][node_ind]-1],
		   value))
	    {
	        return_list.push_back(cell_ind+1);
		--cell_ind;
	    }
	}
    }

  // return vector
    return return_list;
}

//---------------------------------------------------------------------------//
// Overloaded operators
//---------------------------------------------------------------------------//
// overloaded == for design-by-contract

bool CAR_CU_Mesh::operator==(const CAR_CU_Mesh &rhs) const
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
// print out the whole mesh

void CAR_CU_Mesh::print(ostream &out) const
{
    out << endl;
    out << ">>> MESH <<<" << endl;
    out << "============" << endl;

    for (int cell = 1; cell <= num_cells(); cell++)
	print(out, cell);
}

//---------------------------------------------------------------------------//
// print individual cells

void CAR_CU_Mesh::print(ostream &output, int cell) const
{
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

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of CAR_CU_Mesh.cc
//---------------------------------------------------------------------------//
