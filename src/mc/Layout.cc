//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Layout.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 30 15:53:52 1998
 * \brief  Layout class implementation file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Layout.hh"
#include "ds++/Assert.hh"

#include <iomanip>

namespace rtt_mc 
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
  
 * \brief Construct a layout object.

 * This constructor allows the user to set up the data inside of Layout based
 * on the number of cells in the layout. The number of faces per cell are as
 * yet undetermined.

 * \param num_cells the number of cells used to build the layout.

 */
Layout::Layout(int num_cells)
    : face_cell(num_cells)
{
    Ensure (face_cell.size() == num_cells);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor when the problem dimension is known up front.

 * This constructor allows the user to enter the number of faces in
 * each cell of the mesh for which this layout is designed.  The internal
 * data is then appropriately constructed based on the size of the mesh and
 * the number of coarse faces per cell.  A drawback when using this
 * constructor is that each cell must have the same number of coarse faces.
 * However, this is not a limitation of the layout. But, if each cell has a
 * different number of coarse faces then they will have to be sized using the
 * set_size service.

 * Each coarse face is given one fine face by default.

 * \param num_cells number of cells in the layout

 * \param num_faces number of faces per cell

 */
Layout::Layout(int num_cells, int num_faces)
    : face_cell(num_cells, sf_int(num_faces))
{
    Ensure (face_cell.size() == num_cells);
    Ensure (face_cell[num_cells-1].size() == num_faces);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
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
Layout::SP_Pack Layout::pack(const sf_int &current_layout_to_new_layout) const
{
    Require (current_layout_to_new_layout.size() == face_cell.size());

    // number of cells in the layout
    int num_cells = face_cell.size();

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
    vf_int cell_layout(num_packed_cells, sf_int(0));
	
    // build the layout
    int size = 0;
    for (int new_cell, cell = 1; cell <= num_cells; cell++)
    {
	// get cell index for the CURRENT layout
	new_cell = current_layout_to_new_layout[cell-1];

	if (new_cell > 0)
	{
	    Check (new_cell <= num_cells);

	    // set the faces
	    cell_layout[new_cell-1].resize(num_faces(cell));

	    for (int next_n_cell,next_cell,f = 1; f <= num_faces(cell); f++) 
	    {
		// get the cell across the face in the CURRENT mesh
		next_cell = face_cell[cell-1][f-1];

		// if the next cell is less than zero then this layout has
		// boundary cells and we must check that the mapping is one
		// to one
		if (next_cell < 0)
		{
		    Insist (new_cell == cell,
			    "Tried to remap a layout with boundary cells!"); 
		}

		// if this is a vacuum or boundary cell store it
		if (next_cell <= 0)
		{
		    cell_layout[new_cell-1][f-1] = next_cell;
		    size++;
		}

		// else do the conversion to the new mesh (which maybe 1-1)
		else
		{
		    next_n_cell = current_layout_to_new_layout[next_cell-1]; 

		    // add the next new cell to the layout
		    cell_layout[new_cell-1][f-1] = next_n_cell;

		    size++;
		}
	    }
	}
    }
    
    // each cell must have at least 2 faces
    Check (size >= num_packed_cells * 2);

    // add the number of packed cells + a number of faces/cell array to size
    size += 1 + num_packed_cells;

    // make the packed data and add the number of packed cells to the first
    // element; no need to delete the data, Pack will take care of it
    int *data = new int[size];
    data[0]   = num_packed_cells;

    // pack the new layout
    int ctr = 1;
    for (int i = 0; i < num_packed_cells; i++)
    {
	// add the number of faces for this cell
	data[ctr++] = cell_layout[i].size();
	
	// pack the faces
	for(int j = 0; j < cell_layout[i].size(); j++)
	    data[ctr++] = cell_layout[i][j];
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
 * \brief A diagnostic print functions.
 */
void Layout::print(std::ostream &output, int cell_index) const
{
    using std::endl;
    using std::setw;
      
    int faces = num_faces(cell_index);

    // print the face information for cell indexed by cell_index
    output << "=============" << endl;
    output << "Cell   : "     << cell_index << endl;
    output << "# faces: "     << faces << endl;
    output << "-------------" << endl;
    output << setw(4) << "Face" << setw(9) << " Neighbor" << endl;
    output << "-------------" << endl;
    for (int i = 0; i < faces; i++)
        output << setw(4) << i+1 << setw(9)
	       << face_cell[cell_index-1][i] << endl;
    output << "=============" << endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check two Layout objects for equality.
 */
bool Layout::operator==(const Layout &rhs) const
{
    // if the data is equal, the Layouts are equal
    if (face_cell == rhs.face_cell)
	return true;
    
    // if we haven't returned then the Layouts aren't equal
    return false;
}

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

std::ostream& operator<<(std::ostream &output, const Layout &object)
{
    int num_cells = object.num_cells();
    for (int i = 1; i <= num_cells; i++)
        object.print(output, i);
    return output;
}

//===========================================================================//
// LAYOUT::PACK CLASS DEFINITIONS
//===========================================================================//
/*!
 * \brief Constructor.

 * Construct a Layout::Pack instance.  Once allocated integer data is given
 * to the Layout::Pack constructor in the form of an int*, the Pack object
 * owns it.  When the Pack object goes out of scope it will clean up the
 * memory.  In general, Pack objects are only created by calling the
 * Layout::pack() function.

 * \param s size of integer data stream
 * \param d pointer to integer data stream

 */
Layout::Pack::Pack(int s, int *d)
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
Layout::Pack::Pack(const Pack &rhs)
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
Layout::Pack::~Pack()
{
    delete [] data;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack the Layout

 * Unpacks and returns a smart pointer to the new layout.

 * \return smart pointer to the unpacked layout

 */
Layout::SP_Layout Layout::Pack::unpack() const
{
    Require (size >= 1);

    // make a new layout with the appropriate number of cells
    SP_Layout new_layout(new Layout(data[0]));

    // make a reference to the layout for easy operator overloading
    Layout &reflay = *new_layout;

    // index counter
    int ic = 1;

    // unpack the layout
    for (int nf,cell = 1; cell <= reflay.num_cells(); cell++)
    {
	// determine number of face for this cell
	nf = data[ic++];
	Check (nf >= 2);
	
	// resize the layout
	reflay.set_size(cell, nf);

	// assign the faces
	for (int f = 1; f <= nf; f++)
	    reflay(cell,f) = data[ic++];
    }
    
    Ensure (ic == size);
    Ensure (new_layout->num_cells() == data[0]);

    // return the layout
    return new_layout;
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of Layout.cc
//---------------------------------------------------------------------------//
