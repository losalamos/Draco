//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/AMR_Layout.cc
 * \author Thomas M. Evans
 * \date   Tue Jul 18 16:07:35 2000
 * \brief  AMR_Layout member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "AMR_Layout.hh"
#include <iomanip>

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS 
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.

 * This constructor allows the user to set up the data inside of AMR_Layout
 * based on the number of cells in the layout. The number of faces per cell
 * and the number of cells across a face are as yet undetermined.

 * \param num_cells[=0] number of cells in the problem, defaults to zero

 */
AMR_Layout::AMR_Layout(int num_cells)
    : cell_face(num_cells)
{
    Ensure (cell_face.size() == num_cells);
}

//---------------------------------------------------------------------------//

/*!
 * \brief Constructor when the problem dimension is known up front.

 * This constructor allows the user to enter the number of coarse faces in
 * each cell of the mesh for which this layout is designed.  The internal
 * data is then appropriately constructed based on the size of the mesh and
 * the number of coarse faces per cell.  A drawback when using this
 * constructor is that each cell must have the same number of coarse faces.
 * However, this is not a limitation of the layout. But, if each cell has a
 * different number of coarse faces then they will have to be sized using the
 * set_size service.

 * Each coarse face is given one fine face by default.

 * \param num_cells number of cells in the layout

 * \param num_faces number of coarse faces per cell

 */
AMR_Layout::AMR_Layout(int num_cells, int num_faces)
    : cell_face(num_cells, vf_int(num_faces, sf_int(1)))
{
    Ensure (cell_face.size() == num_cells);
    Ensure (cell_face[num_cells-1].size() == num_faces);
    Ensure (cell_face[num_cells-1][num_faces-1].size() == 1);
}

//---------------------------------------------------------------------------//
// PUBLIC MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Pack up a Layout object.

 * This function returns a Layout::Pack object containing a packed Layout.
 * The client specifies a list of cells to include in the packed layout.
 * Each cell in the current layout must be addressed.  For example, if there
 * are 4 cells in the layout, but only 2 are desired in the packed layout,
 * the map should have 10 entries.  Boundary cells must be addressed in a
 * consistent manner.

 * If the layout contains boundary cells (indicated with a negative integer
 * index), the mapping must be one-to-one.  The layout packer has no way to
 * remap boundary cells.

 * \param current_layout_to_new_layout map from the current layout to the new
 * requested, packed layout.

 */
AMR_Layout::SP_Pack 
AMR_Layout::pack(const sf_int &current_layout_to_new_layout) const
{
    Require (current_layout_to_new_layout.size() == cell_face.size());

    // number of cells in the layout
    int num_cells = cell_face.size();

    // count number of cells in new mesh
    int num_packed_cells = 0;
    for (int new_cell, cell = 0; cell < num_cells; cell++)
    {
	new_cell = current_layout_to_new_layout[cell];
	
	if (new_cell > 0)
	{
	    Check (new_cell <= num_cells);
	    num_packed_cells++;
	}
    }

    // allocate space for the new layout
    tf_int layout(num_packed_cells);

    // build the layout
    int size                = 0;
    int new_cell            = 0;
    int next_cell           = 0;
    int next_n_cell         = 0;
    int num_of_coarse_faces = 0;
    for (int num_coarse, cell = 1; cell <= num_cells; cell++)
    {
	// get cell index for the CURRENT layout
	new_cell = current_layout_to_new_layout[cell-1];

	if (new_cell > 0)
	{
	    Check (new_cell <= num_cells);

	    // set the number of coarse faces
	    num_coarse = num_faces(cell);
	    layout[new_cell-1].resize(num_coarse);
	    num_of_coarse_faces += num_coarse;

	    // loop over coarse faces
	    for (int num_fine, cf = 1; cf <= num_coarse; cf++)
	    {
		// how many cells are across the coarse face
		num_fine = num_cells_across(cell, cf);
		layout[new_cell-1][cf-1].resize(num_fine);
		Check (num_fine >= 1);

		// loop over cells across the face
		for (int ff = 1; ff <= num_fine; ff++)
		{
		    // get the cell across the face in the CURRENT mesh
		    next_cell = cell_face[cell-1][cf-1][ff-1];

		    // if the next cell is less than zero then this layout
		    // has boundary cells and we must check that the mapping
		    // is one to one
		    if (next_cell < 0)
		    {
			Insist (new_cell == cell,
				"Tried to remap a layout with boundary cells!"); 
		    }

		    // if this is a vacuum or boundary cell store it
		    if (next_cell <= 0)
		    {
			layout[new_cell-1][cf-1][ff-1] = next_cell;
			size++;
		    }

		    // else do the conversion to the new mesh (which maybe
		    // 1-1)
		    else
		    {
			next_n_cell =
			    current_layout_to_new_layout[next_cell-1];

			// add the next new cell to the layout
			layout[new_cell-1][cf-1][ff-1] = next_n_cell;
			size++;
		    }
		}
	    }
	}
    }
    
    // each cell must have at least 2 faces
    Check (size >= num_packed_cells * 2);

    // add the number of packed cells + a number of faces/cell array to size
    size += 1 + num_packed_cells + num_of_coarse_faces;

    // make the packed data and add the number of packed cells to the first
    // element; no need to delete the data, Pack will take care of it
    int *data = new int[size];
    data[0]   = num_packed_cells;

    // pack the new layout
    int ctr = 1;
    for (int i = 0; i < num_packed_cells; i++)
    {
	// add the number of coarse faces for this cell
	data[ctr++] = layout[i].size();
	
	// loop over coarse faces
	for (int j = 0; j < layout[i].size(); j++)
	{
	    // add the number of fine faces for this coarse face
	    data[ctr++] = layout[i][j].size();

	    // pack the fine face cells
	    for(int k = 0; k < layout[i][j].size(); k++)
		data[ctr++] = layout[i][j][k];
	}
    }
    Check (ctr == size);

    // build the pack object
    SP_Pack pack(new Pack(size, data));

    Ensure (pack->get_size() == ctr);

    // return the pack object
    return pack;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print out the AMR_Layout.
 * \param cell print out layout for given cell index
 */
void AMR_Layout::print(std::ostream &out, int cell) const
{
    using std::endl;
    using std::setw;

    Require (cell > 0 && cell <= cell_face.size());

    int faces = cell_face[cell-1].size();

    // print out the cell info
    out << "====================" << endl;
    out << "Cell         : " << cell << endl;
    if (faces)
	out << "Coarse Faces : " << faces << endl;
    out << "--------------------" << endl;
    for (int i = 0; i < faces; i++)
    {
	out << setw(4) << i+1 << endl;
	for (int j = 0; j < cell_face[cell-1][i].size(); j++)
	    out << setw(8) << j+1 << setw(5) 
		<< cell_face[cell-1][i][j] << endl;
    }
    out << "====================" << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Test for equality.
 */
bool AMR_Layout::operator==(const AMR_Layout &rhs) const
{
    // if the cell data is equal then return true
    if (cell_face == rhs.cell_face)
	return true;

    return false;
}

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out all cells using operator<<.
 */
std::ostream& operator<<(std::ostream &out, const AMR_Layout &object)
{
    int num_cells = object.num_cells();
    for (int i = 1; i <= num_cells; i++)
	object.print(out, i);
    return out;
}

//===========================================================================//
// AMR_LAYOUT::PACK CLASS DEFINITIONS
//===========================================================================//
/*!
 * \brief Constructor.

 * Construct an AMR_Layout::Pack instance.  Once allocated integer data is
 * given to the AMR_Layout::Pack constructor in the form of an int*, the Pack
 * object owns it.  When the Pack object goes out of scope it will clean up
 * the memory.  In general, Pack objects are only created by calling the
 * AMR_Layout::pack() function.

 * \param s size of integer data stream
 * \param d pointer to integer data stream

 */
AMR_Layout::Pack::Pack(int s, int *d)
    : data(d),
      size(s)
{
    // nothing more to do here
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.

 * Do copy construction while preserving memory.  This is not a reference
 * counted class so data is copied from one class to the other during
 * function calls and the like (wherever a copy constructor is called).

 */
AMR_Layout::Pack::Pack(const Pack &rhs)
    : data(new int[rhs.size]),
      size(rhs.size)
{
    // fill up new data array
    for (int i = 0; i < size; i++)
	data[i] = rhs.data[i];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.

 * Cleans up memory when the Pack object goes out of scope.  Once allocated
 * pointers are given to the Pack object the Pack object takes control of
 * them.

 */
AMR_Layout::Pack::~Pack()
{
    delete [] data;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack the AMR_Layout

 * Unpacks and returns a smart pointer to the new layout.

 * \return smart pointer to the unpacked layout

 */
AMR_Layout::SP_Layout AMR_Layout::Pack::unpack() const
{
    Require (size >= 1);

    SP_Layout new_layout(new AMR_Layout(data[0]));
    
    // make a reference to the layout for easy operator overloading
    AMR_Layout &reflay = *new_layout;

    // index counter
    int ic = 1;

    // unpack the layout
    for (int ncoarse, cell = 1; cell <= reflay.num_cells(); cell++)
    {
	// determine the number of coarse faces for this cell
	ncoarse = data[ic++];
	Check (ncoarse >= 2);

	// resize the layout
	reflay.set_size(cell, ncoarse);

	// loop over coarse faces
	for (int nfine, coarse = 1; coarse <= ncoarse; coarse++)
	{
	    // determine the number of fine faces for this cell
	    nfine = data[ic++];
	    Check (nfine >= 1);

	    // resize the layout
	    reflay.set_size(cell, coarse, nfine);

	    // loop over fine faces
	    for (int fine = 1; fine <= nfine; fine++)
		reflay(cell, coarse, fine) = data[ic++];
	}
    }

    Ensure (ic == size);
    Ensure (new_layout->num_cells() == data[0]);

    // return the layout
    return new_layout;
}

} // end namespace rtt_mc


//---------------------------------------------------------------------------//
//                              end of AMR_Layout.cc
//---------------------------------------------------------------------------//
