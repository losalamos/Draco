//----------------------------------*-C++-*----------------------------------//
// Mesh.cc
// B.T. Adams (bta@lanl.gov)
// 18 May 99
/*! 
 * \file   amr_mesh/Mesh.cc
 * \author B.T. Adams
 * \date   Tue May 18 10:33:26 1999
 * \brief  Implementation file for CAR_CU_Mesh class library.
 */
//---------------------------------------------------------------------------//
// @> CAR_CU_Mesh class implementation file (developed from OS_Mesh.cc)
//---------------------------------------------------------------------------//

#include "Mesh.hh"
#include <iomanip>

namespace rtt_amr 
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
/*!
 * \brief Constructs a CAR_CU_Mesh class object (typically using objects that
 *        were created by a CAR_CU_Builder class object).
 * \param layout_ An existing, initialized Layout class object.
 * \param node_coords_ An existing, initialized ncvf_d type object containing
 *                     the mesh node coordinates data.
 * \param cell_nodes_ An existing, initialized ccvf_i type object containing 
 *                   the mesh cell-nodes data.
 * \param generation_ An existing, initialized ccsf_i type object containing 
 *                   the mesh cell generation data.
 * \param submesh_ Status flag indicating that this is a submesh for parallel
 *                 execution (defaults to false). 
 */
CAR_CU_Mesh::CAR_CU_Mesh(Layout & layout_, ncvf_d & node_coords_, 
			 ccvf_i & cell_nodes_,  ccsf_i & generation_,
			 bool submesh_) : layout(layout_), 
                         node_coords(node_coords_), cell_nodes(cell_nodes_), 
                         generation(generation_), submesh(submesh_)
{
  // assertions to verify size of mesh

  // variable initialization
    int ncells = get_num_cells();
    int dimension = get_ndim();
    
  // dimension assertions
    Check (dimension == node_coords.size());

  // mesh size assertions
    Check (ncells == cell_nodes.size());
      

}

//---------------------------------------------------------------------------//
// private member functions
//---------------------------------------------------------------------------//
bool CAR_CU_Mesh::compReal(const double & low_val, const double & high_val)
    const
{
    // require agreement to six significant figures for equality (set by
    // the exponent of EPSILON).
    const double EPSILON = 1.0e-06;
    double epsilon;
    if (low_val != 0 && high_val != 0)
        epsilon = EPSILON * ((fabs(low_val) + fabs(high_val))/2.);
    else
        epsilon = EPSILON;

    return (fabs(high_val - low_val) > epsilon);
}

//---------------------------------------------------------------------------//
// Public member functions
//---------------------------------------------------------------------------//
// do binary search on a cell

/*!
 * \brief Returns the cell number that contains the specified point in space.
 * \param r Coordinate values of point in space.
 * \return Cell number. 
 */
int CAR_CU_Mesh::get_cell_from_coords(const vector<double> &r) const
{
    Require (!submesh);

  // used variables
    int dim         = get_ndim();
    int return_cell;
    int index;

    // binary search of cell node coordinates - node_coords array is sorted 
    // with x varying fastest, followed by y, and then z. find the first 
    // node "below" this point in space.
    int high_index = node_coords[0].size() - 1;
    for (int i = dim; i >= 0; i--)
    {
        int low_index  = 0;
	while ((high_index - low_index) > 1)
	{
	    index = (high_index + low_index) / 2;
	    if (r[i] < node_coords[i][index] && 
		compReal(r[i],node_coords[i][index]))
		high_index = index;
	    else
		low_index  = index;
	}
	high_index = index;
    }
    // have the node_coords index - find the cell that has this node at its
    // negative x,y,z corner (first node_coords), again using a binary search.

    int low_index  = 0;
    high_index = get_num_cells() - 1;
    while ((high_index - low_index) > 1)
    {
        return_cell = (high_index + low_index) / 2;
	if (cell_nodes[return_cell][0] - 1 < index)
	    high_index = return_cell;
	else
	    low_index  = return_cell;
    }
  // return cell index
    return return_cell;
}

//---------------------------------------------------------------------------//
// calculate the distance to boundary

/*!
 * \brief Returns the distance to the nearest boundary in the specified cell 
 *        and direction for the specified point in space.
 * \param r Coordinate values of point in space.
 * \param omega Direction.
 * \param cell Cell number. 
 * \param face Cell face number at the nearest boundary (calculated). 
 * \return Distance to the nearest boundary. 
 */
double CAR_CU_Mesh::get_dist_2_bndry(const vector<double> & r, 
   const vector<double> & omega, int cell, int & face) const
{
    using std::vector;
    using std::min_element;
    const double huge = DBL_MAX;
    
  // calculate distance to the vec(r) boundaries

  // boundary distances over each coordinate direction
    vector<double> dim_dist_boundary(get_ndim(), 0.0);
    
  // loop to get the distances to boundary in each coordinate direction
    for (int i = 0; i < get_ndim(); i++)
    {
      // define absolute dimension index
	int d = i + 1;

      // find the distances to boundary along each dimension
	if (omega[i] == 0.0)
	    dim_dist_boundary[i] = huge;
	else if (omega[i] > 0.0)
	    dim_dist_boundary[i] = 
	        (get_cell_max_coord(d, cell) - r[i]) / omega[i];
	else
	    dim_dist_boundary[i] = 
	        (get_cell_min_coord(d, cell) - r[i]) / omega[i];
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

/*!
 * \brief Returns the face number that corresponds to the specified boundary
 *        side.
 * \param boundary Boundary side (either lox, loy, loz, hix, hiy, or hiz).
 * \param cell Cell (not used).
 * \return Face number. 
 */
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

/*!
 * \brief Returns a list of cells along the specified boundary.
 * \param boundary Boundary (either lox, loy, loz, hix, hiy, or hiz).
 * \return Surface cell numbers. 
 */
vector<int> CAR_CU_Mesh::get_surcells(string boundary) const
{
    Require (!submesh);

    // return a list of cells along the specified boundary

    // make return vector containing a list of cells along specified boundary
    vector<int> return_list;

    // index to the end of the boundary string. Allowed values for the 
    // boundary string are loz, loy, lox, hix, hiy, hiz.
    int end_ind = boundary.size() - 1;
    // the dimension index associated with the boundary surface 
    // (-z = 2, -y = 1, -x = 0, +x = 0, +y = 1, +z = 2);
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
        cell_ind = get_num_cells() - 1;
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
	value = node_coords[dim_ind][cell_nodes[cell_ind][node_ind] - 1];
	for (int k = 0; k < get_num_cells(); k++)
	{
	    if (!compReal(node_coords[dim_ind][cell_nodes[k][node_ind] - 1],
			  value))
		return_list.push_back(k+1);
	}
    }
    // only have to search the beginning (-z) and end (+z) of the last 
    // column of the node_coords vector since they are in ascending order.
    else
    {
        if (boundary[0] == 'l')
	{
	    value = node_coords[dim_ind][cell_nodes[cell_ind][node_ind] - 1];
	    while (!compReal(node_coords[dim_ind]
			     [cell_nodes[cell_ind][node_ind] -1], value))
	    {
	        return_list.push_back(cell_ind+1);
		++cell_ind;
	    }
	}
	else
	{
	    value = node_coords[dim_ind][cell_nodes[cell_ind][node_ind] - 1];
	    while (!compReal(node_coords[dim_ind]
			     [cell_nodes[cell_ind][node_ind]-1], value))
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
// check that a user-/host-defined set of surface source cells actually
// resides on the surface of the system (requires eitehr a vacuum or 
// reflection bnd).

/*!
 * \brief Checks to insure that user-defined surface source cells actually
 *        reside on the specified face and are on the outer boundary of the 
 *        problem.
 * \param ss_face Boundary side (either lox, loy, loz, hix, hiy, or hiz).
 * \param ss_list Surface source cells.
 */
void CAR_CU_Mesh::check_defined_surcells(const string ss_face, 
					 const vector<int> &ss_list) const
{
    // a weak check on number of surface cells
    Check (ss_list.size() <= get_num_cells());

    for (int ss_indx = 0; ss_indx < ss_list.size(); ss_indx++)
    {
        // convert face on which ss resides from string to int.
        // despite its args, get_bndface actually has no cell dependence
	int ss_face_num = get_bndface(ss_face, ss_list[ss_indx]);

        // get bnd condition on ss face; had better be  either vacuum (0) or
	// reflection (cell number).
	int bc = layout(ss_list[ss_indx], ss_face_num, 1);
	Check (bc == 0 || bc == ss_list[ss_indx]);
    }
}

//---------------------------------------------------------------------------//
// Overloaded operators
//---------------------------------------------------------------------------//
// overloaded == for design-by-contract

bool CAR_CU_Mesh::operator==(const CAR_CU_Mesh & rhs) const
{
  // check to see that the Layouts are equal
    if (layout != rhs.layout)
	return false;

  // check the vertices
    if (node_coords != rhs.node_coords)
	return false;
    if (cell_nodes != rhs.cell_nodes)
	return false;

  // if we haven't returned, then the two meshes must be equal
    return true;
}

//---------------------------------------------------------------------------//
// public diagnostic member functions
//---------------------------------------------------------------------------//
// print out the whole mesh

/*!
 * \brief Diagnostic member function used to print out the cell center-point 
 *        coordinate values and widths for the entire mesh.
 * \param output Stream-output class object.
 */
void CAR_CU_Mesh::print(ostream & output) const
{
    output << endl;
    output << ">>> MESH <<<" << endl;
    output << "============" << endl;

    for (int cell = 1; cell <= get_num_cells(); cell++)
	print(output, cell);
}

//---------------------------------------------------------------------------//
// print individual cells

/*!
 * \brief Diagnostic member function used to print out the center-point 
 *        coordinate values and width for the specified cell.
 * \param output Stream-output class object.
 * \param cell Cell number.
 */
void CAR_CU_Mesh::print(ostream & output, int cell) const
{
  // print out content info for 1 cell
    output << "+++++++++++++++" << endl;
    output << "---------------" << endl;
    output << "Cell : "         << cell << endl;
    output << "---------------" << endl;
    output << "Dimensions "     << endl;
    output << "---------------" << endl;
    if (get_ndim() == 2)
    {
	output << " x  : " << get_cell_center_coord(1, cell) << endl;
	output << " y  : " << get_cell_center_coord(2, cell) << endl;
    	output << " dx : " << get_cell_width(1, cell) << endl;
	output << " dy : " << get_cell_width(2, cell) << endl;
    }
    else
    {
	output << " x  : " << get_cell_center_coord(1, cell) << endl;
	output << " y  : " << get_cell_center_coord(2, cell) << endl;
	output << " z  : " << get_cell_center_coord(3, cell) << endl;
    	output << " dx : " << get_cell_width(1, cell) << endl;
	output << " dy : " << get_cell_width(2, cell) << endl;
	output << " dz : " << get_cell_width(3, cell) << endl;
    }	
    output << "---------------" << endl;
    output << "Layout "         << endl;
    output << "---------------" << endl;
    layout.print(output, cell);
    output << "+++++++++++++++" << endl;
}

} // end namespace rtt_amr

//---------------------------------------------------------------------------//
//                              end of Mesh.cc
//---------------------------------------------------------------------------//
