//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Ensight_Translator.t.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:36:10 2000
 * \brief  Ensight_Translator template definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

namespace rtt_viz
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for Ensight_Translator.
 *
 * This constructor builds an Ensight_Translator.  The behavior of the
 * (existing) ensight dump files is controlled by the overwrite parameter.
 * If this is true, any existing ensight dumps (with the same problem name)
 * will be overwritten.  If overwrite is false then the ensight dumps are
 * appended.  The ensight case files are parsed to get the dump times if
 * overwrite is false.
 *
 * \param prefix std_string giving the name of the problem
 * \param gd_wpath directory where dumps are stored
 * \param ens_vdata_names string field containing vertex data names
 * \param ens_cdata_names string field containing cell data names
 * \param overwrite bool that controls whether an existing ensight
 * directory is to be appended to or overwritten.  If true, overwrites the
 * existing ensight directory.  If false, and the ensight directory exists,
 * the case file is appended to.  In either case, if the ensight directory
 * does not exist it is created.  The default for overwrite is false.
 * \param static_geom optional input that if true, geometry is assumed
 * the same across all calls to Ensight_Translator::ensight_dump.
 * \param binary If true, geometry and variable data files are output in
 * binary format.
 *
 * NOTE: If appending data (\a overwrite is false), then \a binary must
 * be the same value as the first ensight dump.  This class does NOT check
 * for this potential error (yes, it's possible to check and is left for a
 * future exercise).
 */
template<class SSF>
Ensight_Translator::Ensight_Translator(const std_string &prefix,
				       const std_string &gd_wpath,
				       const SSF        &ens_vdata_names,
				       const SSF        &ens_cdata_names,
				       const bool        overwrite,
				       const bool        static_geom,
				       const bool        binary)
    : d_static_geom(static_geom)
    , d_binary(binary)
    , d_ens_vdata_names(ens_vdata_names)
    , d_ens_cdata_names(ens_cdata_names)
      
{
    Require (d_dump_times.empty());
    createFilenames(prefix, gd_wpath);
    
    bool graphics_continue = false; // default behavior
    
    if ( ! overwrite ) {
	// then try to parse the case file.  Case files are always ascii.

	std::ifstream casefile(d_case_filename.c_str());

	if ( casefile ) {
	    // then case file exists, so parse the dump times
	    std_string key("number of steps:");
	    std_string line;
	    int num_steps = 0;

	    for (;;) {
		std::getline(casefile, line);
		Insist(casefile.good(),
		       "Error getting number of steps from case file!");
		if ( line.find(key) == 0 ) {
		    std::istringstream ss(line.substr(key.size()));
		    ss >> num_steps;
		    //std::cout << "FOUND " << num_steps
		    //      << " STEPS " << std::endl;
		    break;
		}
	    }

	    // read next three lines and discard
	    std::getline(casefile, line);
	    std::getline(casefile, line);
	    std::getline(casefile, line);
	    Insist(casefile.good(), "Error reading case file!");

	    // read the dump_times

	    d_dump_times.resize(num_steps);
	    
	    for ( int i = 0; i < num_steps; i++ ) {
		casefile >> d_dump_times[i];
		Insist(casefile.good(),
		       "Error reading dump_times from case file!");
		//std::cout << "   STEP " << i
		//	  << " TIME " << d_dump_times[i] << std::endl;
	    }
	    
	    casefile.close();
	    graphics_continue = true;
	}
    }

    initialize(graphics_continue);
}

//---------------------------------------------------------------------------//
// ENSIGHT DUMP PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Do an Ensight dump to disk.
 *
 * Performs an ensight dump in the directory specified in
 * Ensight_Translator::Ensight_Translator().
 *
 * \param icycle time cycle number associated with this dump
 *
 * \param time elapsed problem time
 *
 * \param dt current problem timestep
 *
 * \param ipar IVF field of pointers to vertices.  Dimensioned
 * [0:ncells-1, 0:n_local_vertices_per_cell-1], where
 * n_local_vertices_per_cell is the number of vertices that make up the cell.
 * ipar(i,j) maps the jth+1 vertex number, in the ith+1 cell, to Ensight's
 * "vertex number."  The "vertex number" is in [1:nvertices], so that for
 * example, the corresponding x-coordinate is pt_coor(ipar(i,j)-1, 0).
 *
 * \param iel_type ISF field of Ensight_Cell_Types.  Dimensioned
 * [0:ncells-1].  Each cell in the problem must be associated with a
 * Ensight_Cell_Types enumeration object.
 *
 * \param cell_rgn_index ISF field of region identifiers.  Dimensioned
 * [0:ncells-1].  This matches a region index to each cell in the problem.
 *
 * \param pt_coor FVF field of vertex coordinates. pt_coor is
 * dimensioned [0:nvertices-1, 0:ndim-1].  For each vertex point give the
 * value in the appropriate dimension.
 *
 * \param ens_vrtx_data FVF field of vertex data.  ens_vrtx_data is
 * dimensioned [0:nvertices-1, 0:number of vertex data fields - 1].  The
 * ordering of the second index must match the ens_vdata_names field input
 * argument to Ensight_Translator::Ensight_Translator().  The ordering of the
 * first index must match the vertex ordering from pt_coor.
 *
 * \param ens_cell_data FVF field of cell data.  ens_cell_data is
 * dimensioned [0:ncells-1, 0:number of cell data fields - 1].  The ordering
 * of the second index must match the ens_cdata_names field input argument
 * to Ensight_Translator::Ensight_Translator().  The ordering of the first
 * index must match the cell ordering from ipar.
 *
 * \param rgn_numbers ISF field of unique region ids.  This has dimensions of
 * the number of unique values found in the cell_rgn_index field.
 *
 * \param rgn_name SSF field of unique region names.  This has the same
 * dimensions and ordering as rgn_numbers.  In summary, rgn_numbers gives a
 * list of the unique region ids in the problem and rgn_name gives a list of
 * the names associated with each region id.
 *
 * \sa \ref Ensight_Translator_strings "Ensight_Translator class" for
 * restrictions on name strings.
 *
 * \sa \ref Ensight_Translator_description "Ensight_Translator class" for
 * information on templated field types.
 *
 * \sa Examples page for more details about how to do Ensight dumps.
 */
template<class ISF, class IVF, class SSF, class FVF>
void Ensight_Translator::ensight_dump(int        icycle,
				      double     time,
				      double     dt, 
				      const IVF &ipar_in,
				      const ISF &iel_type,
				      const ISF &cell_rgn_index,
				      const FVF &pt_coor_in,
				      const FVF &ens_vrtx_data_in,
				      const FVF &ens_cell_data_in,
				      const ISF &rgn_numbers,
				      const SSF &rgn_name)
{
    using rtt_traits::Viz_Traits;
    using std::ostringstream;
    using std::string;
    using std::vector;
    using std::find;

    // >>> PREPARE DATA TO SET ENSIGHT OUTPUT

    // Increment local dump counter and add dump time
    d_dump_times.push_back(time);
    int igrdump_num = d_dump_times.size();
    Check (igrdump_num < 10000);

    // load traits for vector field types
    Viz_Traits<IVF> ipar(ipar_in);
    Viz_Traits<FVF> pt_coor(pt_coor_in);
    Viz_Traits<FVF> ens_vrtx_data(ens_vrtx_data_in);
    Viz_Traits<FVF> ens_cell_data(ens_cell_data_in);

    // define sizes used throughout
    int ncells  = ipar.nrows();
    int npoints = pt_coor.nrows();
    int nrgn    = rgn_name.size();

    // Check sizes of all data.
    Require (iel_type.size() == ncells);
    Require (cell_rgn_index.size() == ncells);
    Require (ens_cell_data.nrows() == ncells  || ens_cell_data.nrows() == 0);
    Require (ens_vrtx_data.nrows() == npoints || ens_vrtx_data.nrows() == 0);
    Require (rgn_numbers.size() == nrgn);

    // create ensight postfix indicators
    ostringstream ens_postfix_build;
    ostringstream post_number;

    if (igrdump_num < 10)
	post_number << "000" << igrdump_num;
    else if (igrdump_num < 100) 
	post_number << "00" << igrdump_num;
    else if (igrdump_num < 1000) 
	post_number << "0" << igrdump_num;
    else 
	post_number << igrdump_num;

    ens_postfix_build << "data." << post_number.str();
    string ens_postfix = ens_postfix_build.str();

    // announce the graphics dump
    std::cout << ">>> ENSIGHT GRAPHICS DUMP: icycle= " << icycle 
	      << " time= " << time << " dt= " << dt << std::endl
	      << "dir= " << d_ens_prefix << ", dump_number= " 
	      << igrdump_num << std::endl;

    // create the parts list
    vector<int>::const_iterator find_location_c;
    vector<int>::iterator       find_location;
    vector<int>                 parts_list;

    for (int i = 0; i < ncells; i++)
    {
	find_location = find(parts_list.begin(), parts_list.end(),
			     cell_rgn_index[i]); 

	if (find_location == parts_list.end())
	    parts_list.push_back(cell_rgn_index[i]);
    }
    
    // store the number of parts
    int nparts = parts_list.size();
    
    // create the parts names
    vector<string> part_names;
    
    for (int i = 0; i < nparts; i++)
    {
	find_location_c = find(rgn_numbers.begin(), rgn_numbers.end(),
			       parts_list[i]);

	if (find_location_c != rgn_numbers.end())
	{
	    int index = find_location_c - rgn_numbers.begin();
	    part_names.push_back(rgn_name[index]);
	}
	else if (find_location_c == rgn_numbers.end())
	{
	    Insist (0, "Didn't supply a region name!");
	}
    }

    Insist (parts_list.size() == part_names.size(), "Mismatch on part size!");
    Insist (rgn_name.size() == parts_list.size(), "Mismatch on region size!");

    // create the cells that make up each part

    // vertices_of_part[ipart] is the set of vertex indices that make up part
    // ipart.
    vec_set_int vertices_of_part(nparts);

    // cells_of_type[ipart][itype][i] is the cell index of the i'th cell of
    // type itype in part ipart.
    sf3_int cells_of_type(nparts);
    for ( int i = 0; i < nparts; i++ )
	cells_of_type[i].resize(d_num_ensight_cell_types);

    // Initialize cells_of_type and vertices_of_part.

    for ( int i = 0; i < ncells; i++ )
    {
	find_location = find(parts_list.begin(), parts_list.end(),
			     cell_rgn_index[i]);

	Check(find_location != parts_list.end());
	Check(iel_type[i] < d_num_ensight_cell_types);

	int ipart = find_location - parts_list.begin();

	cells_of_type[ipart][iel_type[i]].push_back(i);

	int nvertices = d_vrtx_cnt[iel_type[i]];
	
	for ( int iv = 0; iv < nvertices; ++iv )
	    vertices_of_part[ipart].insert(ipar(i, iv)-1);
    }

    // >>> WRITE OUT DATA TO DIRECTORIES
    
    // write time to case file
    ensight_case();

    // WRITE THE GEOMETRY FILE
    if ( (! d_static_geom ) ||
	 (d_dump_times.size() == 1) ) {
	ensight_geom(ens_postfix, icycle, time, dt, ipar, pt_coor,
		     part_names, cells_of_type, vertices_of_part);
    }

    // write the vertex data
    if (ens_vrtx_data.nrows() > 0)
	ensight_vrtx_data(ens_postfix, ens_vrtx_data, vertices_of_part);

    // write out the cell data
    if (ens_cell_data.nrows() > 0)
	ensight_cell_data(ens_postfix, ens_cell_data, cells_of_type,
			  part_names); 
}

//---------------------------------------------------------------------------//
// ENSIGHT DATA OUTPUT FUNCTIONS (PRIVATE)
//---------------------------------------------------------------------------//
/*!
 * \brief Write out data to ensight geometry file.
 */
template<class IVF, class FVF> void 
Ensight_Translator::ensight_geom(
    const std_string                  &ens_postfix,
    const int                          icycle,
    const double                       time,
    const double                       dt,
    const rtt_traits::Viz_Traits<IVF> &ipar, 
    const rtt_traits::Viz_Traits<FVF> &pt_coor,
    const sf_string                   &part_names,
    const sf3_int                     &cells_of_type, 
    const vec_set_int                 &vertices_of_part)
{
    using rtt_viz::endl;
    using std::string;
    using std::ostringstream;

    // make output file for this timestep
    string filename  = d_geo_dir + "/";
    if ( d_static_geom ) {
	filename += "data";
    }
    else {
	filename += ens_postfix;
    }

    Ensight_Stream geomout(filename, d_binary, true);

    // write the header
    geomout << "Description line 1" << endl;
    
    ostringstream s;
    s << "probtime " << time << " cycleno " << icycle;
    geomout << s.str() << endl;
    
    geomout << "node id given" << endl;
    geomout << "element id given" << endl;
    
    int ndim = pt_coor.ncols(0);

    // write the data for each part

    for (int ipart = 0; ipart < part_names.size(); ipart++)
    {
	// output part number and names
 	geomout << "part" << endl;
 	geomout << ipart+1 << endl;
	geomout << part_names[ipart] << endl;

	// output the global vertex indices and form local_vertex.
	// local_vertex maps the global vertex index to the local vertex
	// index, for this particular part.
	geomout << "coordinates" << endl;
	const set_int &v = vertices_of_part[ipart];
	geomout << v.size() << endl;
	std::map<int, int> local_vertex;
	int count = 1;
	for ( set_const_iterator iv = v.begin(); iv != v.end(); ++iv )
	{
	    geomout << *iv+1 << endl;
	    local_vertex[*iv+1] = count;
	    ++count;
	}

	// output the coordinates
	for ( int idim = 0; idim < ndim; idim++ )
	    for ( set_const_iterator iv = v.begin(); iv != v.end(); ++iv )
		geomout << pt_coor(*iv, idim) << endl;

	// ensight expects coordinates for three dimensions, so fill any
	// remaining dimensions with zeroes
	double zero = 0.0;
	for ( int idim = ndim; idim < 3; idim++ )
	    for ( set_const_iterator iv = v.begin(); iv != v.end(); ++iv )
		geomout << zero << endl;

	// for each cell type, dump the local vertex indices for each cell.
	for (int type = 0; type < d_num_ensight_cell_types; type++)
	{
	    const sf_int &c = cells_of_type[ipart][type];
	    const int num_elem = c.size();
	    
	    if (num_elem > 0)
	    {
		geomout << d_ensight_cell_names[type] << endl;
		geomout << num_elem << endl;

		for ( int i = 0; i < num_elem; ++i )
		    geomout << c[i] << endl;
		
		for ( int i = 0; i < num_elem; ++i )
		{
		    for (int j = 0; j < d_vrtx_cnt[type]; j++)
			geomout << local_vertex[ipar(c[i],j)];
		    geomout << endl;
		}
	    }
	} // done looping over cell types
    } // done looping over parts
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write out data to ensight vertex data.
 */
template<class FVF> 
void Ensight_Translator::ensight_vrtx_data(
    const std_string                  &ens_postfix,
    const rtt_traits::Viz_Traits<FVF> &ens_vrtx_data,
    const vec_set_int                 &vertices_of_part)
{
    using rtt_viz::endl;
    using std::string;

    const int nparts = vertices_of_part.size();

    // loop over all vertex data fields and write out data for each field
    for (int nvd = 0; nvd < ens_vrtx_data.ncols(0); nvd++)
    {
	// open file for this data
	string filename  = d_vdata_dirs[nvd] + "/" + ens_postfix;
	Ensight_Stream vout(filename, d_binary);

	vout << d_ens_vdata_names[nvd] << endl;

	// loop over parts and output the data for this particular field
	for ( int ipart = 0; ipart < nparts; ++ipart )
	{
	    vout << "part" << endl;
	    vout << ipart+1 << endl;
	    vout << "coordinates" << endl;

	    const set_int &v = vertices_of_part[ipart];
	    for ( set_const_iterator iv = v.begin(); iv != v.end(); ++iv )
		vout << ens_vrtx_data(*iv, nvd) << endl;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write out data to ensight cell data.
 */
template<class FVF> 
void Ensight_Translator::ensight_cell_data(
    const std_string                  &ens_postfix,
    const rtt_traits::Viz_Traits<FVF> &ens_cell_data,
    const sf3_int                     &cells_of_type,
    const sf_string                   &part_names)
{
    using rtt_viz::endl;
    using std::string;

    // loop over all vertex data fields and write out data for each field
    for (int ncd = 0; ncd < ens_cell_data.ncols(0); ncd++)
    {
	// open file for this data
	string filename  = d_cdata_dirs[ncd] + "/" + ens_postfix;
	Ensight_Stream cellout(filename, d_binary);

	cellout << d_ens_cdata_names[ncd] << endl;

	for (int ipart = 0; ipart < part_names.size(); ipart++)
	{
	    cellout << "part" << endl;
	    cellout << ipart+1 << endl;
	    
	    // loop over ensight cell types
	    for (int type = 0; type < d_num_ensight_cell_types; type++)
	    {
		const sf_int &c = cells_of_type[ipart][type];
		
		int num_elem = c.size();

		// print out data if there are cells of this type
		if (num_elem > 0)
		{
		    // printout cell-type name
		    cellout << d_ensight_cell_names[type] << endl;

		    // print out data
		    for (int i = 0; i < num_elem; i++)
			cellout << ens_cell_data(c[i], ncd) << endl;
		}
	    } 
	}
    }
}

} // end of rtt_viz

//---------------------------------------------------------------------------//
//                        end of viz/Ensight_Translator.t.hh
//---------------------------------------------------------------------------//
