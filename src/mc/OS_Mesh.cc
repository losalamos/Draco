//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/OS_Mesh.cc
 * \author Thomas M. Evans
 * \date   Tue Feb  3 16:50:13 1998
 * \brief  OS_Mesh class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "XYCoord_sys.hh"
#include "XYZCoord_sys.hh"
#include "Constants.hh"
#include "viz/Ensight_Translator.hh"
#include <iomanip>

namespace rtt_mc 
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!

 * \brief OS_Mesh constructor.

 * \param coord_  coordinate system smart pointer
 * \param layout_ Layout object describing mesh layout
 * \param vertex_ list of mesh vertices
 * \param cell_pair_ list of cells paired with vertices
 * \param submesh_ boolean indicating if this is a submesh

 */
OS_Mesh::OS_Mesh(SP_Coord_sys coord_, 
		 Layout &layout_, 
		 vf_double &vertex_, 
		 vf_int &cell_pair_, 
		 bool submesh_) 
    : coord(coord_), layout(layout_), vertex(vertex_),
      cell_pair(cell_pair_), sur(coord->get_dim()), submesh(submesh_)
{
    // assertions to verify size of mesh and existence of a Layout and
    // Coord_sys  
    Require (coord);
	
    // variable initialization
    int ncells    = num_cells();
    int dimension = coord->get_dim();
    
    // dimension assertions
    Check (dimension == vertex.size()); 
    Check (dimension == sur.size()); 
    
    // mesh size assertions
    Check (ncells == cell_pair.size());
      
    // calculate surface array
    if (!submesh)
	calc_surface();
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
// calculate an array of the dimensional surfaces which make up the OS_Mesh

void OS_Mesh::calc_surface()
{
    using std::vector;
    using std::sort;

    // initialize mesh_size for assertion at end of function
    int mesh_size = 1;

    // loop to calculate surface array
    for (int d = 0; d < coord->get_dim(); d++)
    {
	// define an array for dim which is sorted in ascending order
	vector<double> sorted = vertex[d];
	sort(sorted.begin(), sorted.end());

	// loop over sorted array, appending new surfaces onto sur array, watch 
	// out for possible machine error (especially when merging host codes)
	// in the sorted[i] > sorted[i-1] comparison!!!
	sur[d].push_back(sorted[0]);
	for (int i = 1; i < sorted.size(); i++)
	    if (sorted[i] > sorted[i-1])
		sur[d].push_back(sorted[i]);

	// calculate mesh_size by dimension
	mesh_size *= (sur[d].size() - 1);
    }
    
    // assert mesh size
    Require (num_cells() == mesh_size);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FOR IMC
//---------------------------------------------------------------------------//
// do binary search on a cell

int OS_Mesh::get_cell(const sf_double &r) const
{
    Require (!submesh);

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

double OS_Mesh::get_db(const sf_double &r, const sf_double &omega, int cell, 
		       int &face) const
{
    using std::min_element;
    using std::vector;
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
	face = 1 + 2 * index;
    else
	face = 2 + 2 * index;

    // return the distance-to-boundary
    return dist_boundary;
}

//---------------------------------------------------------------------------//
// return the face number for a given cell boundary

int OS_Mesh::get_bndface(std_string boundary, int cell) const
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

OS_Mesh::sf_int OS_Mesh::get_surcells(std::string boundary) const
{
    using std::vector;

    Require (!submesh);

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
// check that a user-/host-defined set of surface source cells actually
// resides on the surface of the system (requires a vacuum bnd).

bool OS_Mesh::check_defined_surcells(const std_string ss_face, 
				     const sf_int &ss_list) const
{
    // a weak check on number of surface cells
    Check (ss_list.size() <= num_cells());

    for (int ss_indx = 0; ss_indx < ss_list.size(); ss_indx++)
    {
        // convert face on which ss resides from string to int.
        // despite its args, get_bndface actually has no cell dependence
	int ss_face_num = get_bndface(ss_face, ss_list[ss_indx]);

        // get bnd condition on ss face; had better be vacuum (0)
	int bc = layout(ss_list[ss_indx], ss_face_num);
	if (bc != 0) 
	    return false;
    }

    return true;
}

//---------------------------------------------------------------------------//
// MESH PACKING INTERFACE
//---------------------------------------------------------------------------//
/*!
  
 * \brief Pack up a mesh into a Pack struct for communication and
 * persistence.

 * The cell list provides the cells to pack up.  It also is a map from a full
 * mesh to a spatially decomposed mesh.  Thus, the packed mesh will only
 * contain cells in the cell list with the provided mappings.

 * The packer will not produce an "exact" copy of the mesh even if the
 * current_mesh_to_new_mesh mapping is one to one.  It will produce an
 * equivalent copy (the internal data will be organized differently), and
 * operator== will fail on such a comparison.  To produce an exact copy, call
 * pack without any arguments.
 
 * \param current_mesh_to_new_mesh list of cells to include in the packed
 * mesh, set this to NULL to produce an exact copy

 */
OS_Mesh::SP_Pack OS_Mesh::pack(const sf_int &current_mesh_to_new_mesh) const
{
    Require (current_mesh_to_new_mesh.size() <= layout.num_cells() ||
	     current_mesh_to_new_mesh.size() == 0);

    // determine if the mesh is replicated or not
    sf_int current_to_new_replicate;
    bool   replicate;
    if (current_mesh_to_new_mesh.size() == 0)
    {
	replicate = true;
	current_to_new_replicate.resize(layout.num_cells());

	for (int cell = 1; cell <= layout.num_cells(); cell++)
	    current_to_new_replicate[cell-1] = cell;
    }
    else
    {
	replicate = false;
    }

    // determine coord system
    int coord_indicator = 0;
    if (coord->get_Coord() == "xy")
	coord_indicator = 1;
    else if (coord->get_Coord() == "xyz")
	coord_indicator = 2;
    else
    {
	Insist (0, "Improper coordinate system specified for OS_Mesh.");
    }

    // packed layout
    Layout::SP_Pack packed_layout;
    int             layout_size;
    int             num_packed_cells;

    // packup the vertices and cell-pairs
    char *mesh_data      = 0;
    int   mesh_data_size = 0;
    if (replicate)
    {
	// packup the layout
	packed_layout = layout.pack(current_to_new_replicate);
	layout_size   = packed_layout->get_size();
	Check (layout_size >= 1);

	// number of packed cells in this mesh
	num_packed_cells = packed_layout->get_num_packed_cells();
	Check (num_packed_cells == layout.num_cells());

	// pack up the mesh data
	mesh_data = pack_mesh_data(mesh_data_size, vertex, cell_pair);
    }
    else
    {
	// packup the layout
	packed_layout = layout.pack(current_mesh_to_new_mesh);
	layout_size   = packed_layout->get_size();
	Check (layout_size >= 1);

	// number of packed cells in this mesh
	num_packed_cells = packed_layout->get_num_packed_cells();
	Check (num_packed_cells <= layout.num_cells());

	// make local vertex and cellpair data for the packed array
	vf_double local_vertex(vertex.size());
	vf_int    local_cellpair(num_packed_cells);
	pack_compressed(current_mesh_to_new_mesh, local_vertex,
			local_cellpair);

	// pack up the mesh data
	mesh_data = pack_mesh_data(mesh_data_size, local_vertex,
				   local_cellpair);
    }
    Check (num_packed_cells >= 0 && num_packed_cells <= layout.num_cells());
    Check (mesh_data_size >= 2 * sizeof(int));
    Check (mesh_data != 0);

    // now pack up the mesh

    // ints (1-coordsys; 1-size of packed layout; 1-size of mesh data; 
    // 1-num packed cells; packed layout)
    int total_ints  = (4 + packed_layout->get_size()) * sizeof(int);
    int total_chars = mesh_data_size;

    int   size = total_ints + total_chars;
    char *data = new char[size];

    // iterator for packing int data
    const char *itor = 0;
    int          ctr = 0;

    // pack up the mesh

    // pack up the number of packed cells
    itor = reinterpret_cast<const char *>(&num_packed_cells);
    for (int i = 0; i < sizeof(int); i++)
	data[ctr++] = itor[i];

    // pack up the coord indicator
    itor = reinterpret_cast<const char *>(&coord_indicator);
    for (int i = 0; i < sizeof(int); i++)
	data[ctr++] = itor[i];

    // pack up the layout size
    itor = reinterpret_cast<const char *>(&layout_size);
    for (int i = 0; i < sizeof(int); i++)
	data[ctr++] = itor[i];

    // pack up the layout
    itor = reinterpret_cast<const char *>(packed_layout->begin());
    for (int i = 0; i < packed_layout->get_size() * sizeof(int); i++)
	data[ctr++] = itor[i];

    // pack up the mesh size
    itor = reinterpret_cast<const char *>(&mesh_data_size);
    for (int i = 0; i < sizeof(int); i++)
	data[ctr++] = itor[i];

    // pack up the mesh data
    for (int i = 0; i < mesh_data_size; i++)
	data[ctr++] = mesh_data[i];

    Ensure (ctr == size);

    // clean up some memory
    delete [] mesh_data;
    
    // make a packed mesh
    SP_Pack packed_mesh(new OS_Mesh::Pack(size, data));

    Ensure (packed_mesh->get_num_packed_cells() == num_packed_cells);
    return packed_mesh;
}

//---------------------------------------------------------------------------//
// MESH PACKING IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Pack up mesh data and return a pointer to it.

 * This function packs up the vertex and cell pair data given to it.

 */
char* OS_Mesh::pack_mesh_data(int &size,
			      const vf_double &local_vertex,
			      const vf_int    &local_cellpair) const
{
    Require (size == 0);
    Require (local_vertex.size() >= 1);
    Require (local_cellpair.size() <= layout.num_cells());

    // get sizes needed for mesh
    int num_vtcs = local_vertex.size() * local_vertex[0].size();
    int num_cp   = 0;

    for (int i = 0; i < local_cellpair.size(); i++)
	num_cp += (1 + local_cellpair[i].size());

    // total size (add 2 ints for sizes)
    int size_v  = num_vtcs * sizeof(double);
    int size_cp = num_cp * sizeof(int);
    size        = 2 * sizeof(int) + size_v + size_cp;
    char *data  = new char[size];

    // write out the vertices
    int vctr     = 0;
    double *vtcs = new double[num_vtcs];
    for (int i = 0; i < local_vertex.size(); i++)
	for (int j = 0; j < local_vertex[i].size(); j++)
	    vtcs[vctr++] = local_vertex[i][j];
    Check (vctr == num_vtcs);

    // write out the cell pairs
    int cpctr = 0;
    int *cp   = new int[num_cp];
    for (int i = 0; i < local_cellpair.size(); i++)
    {
	cp[cpctr++] = local_cellpair[i].size();
	for (int j = 0; j < local_cellpair[i].size(); j++)
	    cp[cpctr++] = local_cellpair[i][j];
    }
    Check (cpctr == num_cp);

    // collapse the vertices and cell pairs into a char array
    const char *itor = 0;
    int ctr          = 0;

    // first add the vertices

    // add the vertices size
    itor = reinterpret_cast<const char *>(&size_v);
    for (int i = 0; i < sizeof(int); i++)
	data[ctr++] = itor[i];

    // add the vertices array
    itor = reinterpret_cast<const char *>(vtcs);
    for (int i = 0; i < size_v; i++)
	data[ctr++] = itor[i];

    // add the cell pair size
    itor = reinterpret_cast<const char *>(&size_cp);
    for (int i = 0; i < sizeof(int); i++)
	data[ctr++] = itor[i];

    // add the cell pairs
    itor = reinterpret_cast<const char *>(cp);
    for (int i = 0; i < size_cp; i++)
	data[ctr++] = itor[i];
    
    // clean up memory
    delete [] vtcs;
    delete [] cp;
    
    Ensure (ctr == size);
    return data;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack up mesh data for compressed mesh.

 * Note that this will not produce an "exact" copy of the mesh even if the
 * mapping is one to one.

 */
void OS_Mesh::pack_compressed(const sf_int &current_new, 
			      vf_double &local_vertex,
			      vf_int &local_cellpair) const
{
    Require (local_vertex.size() >= 1);
    Require (current_new.size() == layout.num_cells());
    Require (local_vertex.size() == vertex.size());
    Require (local_cellpair.size() <= layout.num_cells());

    // sizes
    int num_packed     = local_cellpair.size();
    int num_total_vert = 0;

    // vertex map
    sf_int vert_map(vertex[0].size(), 0);

    // loop through the cells and build the new vertex and cell pair data
    for (int nvert, ncell, cell = 0; cell < current_new.size(); cell++)
    {
	// determine the new cell index
	ncell = current_new[cell];

	// if the new cell index is > 0 then we need to make a cell-pair and
	// vertices entry for it
	if (ncell > 0)
	{
	    Check (ncell <= num_packed);
	    Check (local_cellpair[ncell-1].size() == 0);

	    // determine the number of vertices this cell has
	    nvert = cell_pair[cell].size();

	    // allocate this in the new cell pair vector
	    local_cellpair[ncell-1].resize(nvert);

	    // add these vertices to the new cell pair and to the new
	    // vertices array
	    for (int v, i = 0; i < nvert; i++)
	    {
		// get the vertex index
		v = cell_pair[cell][i];
		Check (v > 0 && v <= vert_map.size());

		// calculate a new vertex index 
		if (vert_map[v-1] == 0)
		{
		    num_total_vert++;
		    vert_map[v-1] = num_total_vert;

		    // add the new vertex to the vertex array
		    for (int d = 0; d < vertex.size(); d++)
		    {
			local_vertex[d].push_back(vertex[d][v-1]);
			Check (local_vertex[d].size() == num_total_vert);;
		    }
		} 

		// add this vertex to the new cell_pair array
		local_cellpair[ncell-1][i] = vert_map[v-1];
	    }
	}
    }
    Check (local_vertex.size() == vertex.size());
    Check (local_vertex[0].size() == num_total_vert); 
}

//---------------------------------------------------------------------------//
// MEMBER FUNCTION OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
// overloaded == for design-by-contract

bool OS_Mesh::operator==(const OS_Mesh &rhs) const
{
    // check to see that we have the same type of coordinate systems
    if (*coord != *rhs.coord)
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
// INTERFACE FOR GRAPHIC DUMPS
//---------------------------------------------------------------------------//
// return the cell type for each cell in the mesh

OS_Mesh::sf_int OS_Mesh::get_cell_types() const
{
    using std::vector;

    vector<int> cell_type(layout.num_cells());

    if (coord->get_dim() == 2)
	std::fill(cell_type.begin(), cell_type.end(),
		  rtt_viz::four_node_quadrangle);
	
    else if (coord->get_dim() == 3)
	std::fill(cell_type.begin(), cell_type.end(),
		  rtt_viz::eight_node_hexahedron);

    return cell_type;
}

//---------------------------------------------------------------------------//
// get point coordinates [0:npoints-1, 0:ndim-1]

OS_Mesh::vf_double OS_Mesh::get_point_coord() const
{
    using std::vector;

    int npoints = vertex[0].size();
    vector<vector<double> > return_coord(npoints);
    for (int i = 0; i < return_coord.size(); i++)
    {
	return_coord[i].resize(coord->get_dim());
	for (int j = 0; j < return_coord[i].size(); j++)
	{
	    Check (return_coord[i].size() == vertex.size());
	    Check (vertex[j].size() == return_coord.size());

	    return_coord[i][j] = vertex[j][i];
	}
    }

    return return_coord;
}

//---------------------------------------------------------------------------//
// DIAGNOSTIC FUNCTIONS (PRINTS)
//---------------------------------------------------------------------------//
// print out the whole mesh

void OS_Mesh::print(std::ostream &out) const
{
    using std::endl;

    out << endl;
    out << ">>> MESH <<<" << endl;
    out << "============" << endl;

    for (int cell = 1; cell <= num_cells(); cell++)
	print(out, cell);
}

//---------------------------------------------------------------------------//
// print individual cells

void OS_Mesh::print(std::ostream &output, int cell) const
{
    using std::endl;

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

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

std::ostream& operator<<(std::ostream &output, const OS_Mesh &object)
{
    object.print(output);
    return output;
}

//===========================================================================//
// OS_MESH::PACK DEFINITIONS
//===========================================================================//
/*!
 * \brief Constructor.

 * Construct a OS_Mesh::Pack instance.  Once allocated mesh data is given to
 * the OS_Mesh::Pack constructor in the form of a char*, the Pack object owns
 * it.  When the Pack object goes out of scope it will clean up the memory.
 * In general, Pack objects are only created by calling the OS_Mesh::pack()
 * function.

 * \param s size of char data stream
 * \param d pointer to char data stream

 */
OS_Mesh::Pack::Pack(int s, char *d)
    : data(d),
      size(s)
{
    Require (size >= 4 * sizeof(int));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.

 * Do copy construction while preserving memory.  This is not a reference
 * counted class so data is copied from one class to the other during
 * function calls and the like (wherever a copy constructor is called).

 */
OS_Mesh::Pack::Pack(const Pack &rhs)
    : data(new char[rhs.size]),
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
OS_Mesh::Pack::~Pack()
{
    delete [] data;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get number of cells in the packed mesh.
 */
int OS_Mesh::Pack::get_num_packed_cells() const
{
    Require (size >= 4 * sizeof(int));

    int   num_cells = 0;
    char *itor      = reinterpret_cast<char*>(&num_cells);

    for (int i = 0; i < sizeof(int); i++)
	itor[i] = data[i];
    
    Ensure (num_cells >= 0);
    return num_cells;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack the OS_Mesh.

 * Unpacks and returns a smart pointer to the new OS_Mesh.

 * \return smart pointer to the unpacked mesh

 */
OS_Mesh::SP_Mesh OS_Mesh::Pack::unpack() const
{
    using rtt_dsxx::SP;

    Require (size >= 4 * sizeof(int));

    // counter for unpacking
    int ctr = 0;

    // iterator for unpacking
    char *itor = 0;

    // determine the number of packed cells
    int num_packed_cells = 0;
    itor                 = reinterpret_cast<char *>(&num_packed_cells);
    for (int i = 0; i < sizeof(int); i++)
	itor[i] = data[ctr++];
    Check (num_packed_cells >= 0);

    // >>> UNPACK THE COORD SYS
    int coord_indicator = 0;
    itor                = reinterpret_cast<char *>(&coord_indicator);
    SP<Coord_sys> coord;
    for (int i = 0; i < sizeof(int); i++)
	itor[i] = data[ctr++];

    if (coord_indicator == 1)
	coord = new XYCoord_sys();
    else if (coord_indicator == 2)
	coord = new XYZCoord_sys();
    else
    {
	Insist (0, "Improper coordinate system specified for OS_Mesh.");
    }

    // UNPACK THE LAYOUT
    int layout_size  = 0;
    itor             = reinterpret_cast<char *>(&layout_size);
    for (int i = 0; i < sizeof(int); i++)
	itor[i] = data[ctr++];

    // don't need to reclaim this memory because we are giving it to the
    // layout packer
    int *layout_data = new int[layout_size];
    itor             = reinterpret_cast<char *>(layout_data);
    for (int i = 0; i < layout_size * sizeof(int); i++)
	itor[i] = data[ctr++];

    Layout::Pack packed_layout(layout_size, layout_data);
    SP<Layout> layout = packed_layout.unpack();
    Check (layout->num_cells() == num_packed_cells);

    // UNPACK THE MESH DATA
    
    // get the size of the mesh data
    int mesh_data_size = 0;
    itor               = reinterpret_cast<char *>(&mesh_data_size);
    for (int i = 0; i < sizeof(int); i++)
	itor[i] = data[ctr++];
    Check (size >= 2 * sizeof(int));

    // get the size of the vertices in bytes
    int vertx_size = 0;
    itor           = reinterpret_cast<char *>(&vertx_size);
    for (int i = 0; i < sizeof(int); i++)
	itor[i] = data[ctr++];
    Check (vertx_size >= 0);
    Check (vertx_size ? vertx_size % sizeof(double) == 0 : true);
    Check (vertx_size ? vertx_size % coord->get_dim() == 0 : true);

    // get the vertices (we can't use iterators to a vector here because we
    // cannot assume that they can be cast to a char *)
    int num_v_dbl = vertx_size / sizeof(double);
    double *vertx = new double[num_v_dbl];
    itor          = reinterpret_cast<char *>(vertx);
    for (int i = 0; i < vertx_size; i++)
	itor[i] = data[ctr++];
    
    // make the vertex data
    int num_verts = num_v_dbl / coord->get_dim();
    int vctr      = 0;
    vf_double vertex(coord->get_dim(), sf_double(num_verts));
    for (int i = 0; i < vertex.size(); i++)
	for (int j = 0; j < vertex[i].size(); j++)
	    vertex[i][j] = vertx[vctr++];
    Check (vctr == num_v_dbl);

    // reclaim memory
    delete [] vertx;
    
    // get the size of the cell pairs in bytes
    int cp_size = 0;
    itor        = reinterpret_cast<char *>(&cp_size);
    for (int i = 0; i < sizeof(int); i++)
	itor[i] = data[ctr++];
    Check (cp_size >= num_packed_cells);
    Check (cp_size ? cp_size % sizeof(int) == 0 : true);

    // get the cell pairs
    int num_cp_int = cp_size / sizeof(int);
    int *cp        = new int[num_cp_int];
    itor           = reinterpret_cast<char *>(cp);
    for (int i = 0; i < cp_size; i++)
	itor[i] = data[ctr++];

    // make the cell pair data
    int cpctr      = 0;
    int num_v_cell = 0;
    vf_int cell_pair(num_packed_cells);
    for (int i = 0; i < cell_pair.size(); i++)
    {
	// resize for this cell
	num_v_cell = cp[cpctr++];
	cell_pair[i].resize(num_v_cell);
	Check (num_v_cell == std::pow(static_cast<double>(2),
				      coord->get_dim())); 
	for (int j = 0; j < cell_pair[i].size(); j++)
	    cell_pair[i][j] = cp[cpctr++];
    }
    Check (cpctr == num_cp_int);

    // reclaim memory
    delete [] cp;

    // make sure the count is accurate
    Ensure (ctr == size);

    // build the new mesh
    SP_Mesh unpacked_mesh(new OS_Mesh(coord, *layout, vertex, cell_pair,
				      true));
    
    Ensure (unpacked_mesh->num_cells() == num_packed_cells);
    Ensure (unpacked_mesh->get_spatial_dimension() == coord->get_dim());
    Ensure (!unpacked_mesh->full_Mesh());

    return unpacked_mesh;
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of OS_Mesh.cc
//---------------------------------------------------------------------------//
