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
 * The constructor automatically knows the number of Ensight cell types (15),
 * and it sets data dependent upon this variable appropriately.
 *
 * Ensight data directories are created in the constructor based upon the
 * values of prefix_in and gd_wpath and the data names in ens_vdata_names_in
 * and ens_cdata_names_in.  Files are overwritten at runtime.  However,
 * directories are not purged.
 *
 * \param prefix_in std_string giving the name of the problem
 * \param gd_wpath directory where dumps are stored
 * \param ens_vdata_names_in string field containing vertex data names
 * \param ens_cdata_names_in string field containing cell data names 
 *
 * \sa \ref Ensight_Translator_description "Ensight_Translator class" 
 * for details about SSF templated field types.  
 *
 * \sa \ref Ensight_Translator_strings "Ensight_Translator class" for
 * details about restrictions on names (strings).
 */
template<class SSF>
Ensight_Translator::Ensight_Translator(const std_string &prefix_in,
				       const std_string &gd_wpath,
				       const SSF &ens_vdata_names_in,
				       const SSF &ens_cdata_names_in)
    : num_ensight_cell_types(15),
      ensight_cell_names(num_ensight_cell_types),
      vrtx_cnt(num_ensight_cell_types),
      cell_type_index(num_ensight_cell_types),
      ntime_current(0),
      dump_times(0),
      igrdump_num(0),
      prefix(prefix_in),
      ens_vdata_names(ens_vdata_names_in.size()),
      ens_cdata_names(ens_cdata_names_in.size()),
      vdata_dirs(ens_vdata_names.size()),
      cdata_dirs(ens_cdata_names.size())
{
    // assign values to ensight_cell_names
    ensight_cell_names[0]  = "point";
    ensight_cell_names[1]  = "bar2";
    ensight_cell_names[2]  = "bar3";
    ensight_cell_names[3]  = "tria3";
    ensight_cell_names[4]  = "tria6";
    ensight_cell_names[5]  = "quad4";
    ensight_cell_names[6]  = "quad8";
    ensight_cell_names[7]  = "tetra4";
    ensight_cell_names[8]  = "tetra10";
    ensight_cell_names[9]  = "pyramid5";
    ensight_cell_names[10] = "pyramid13";
    ensight_cell_names[11] = "hexa8";
    ensight_cell_names[12] = "hexa20";
    ensight_cell_names[13] = "penta6";
    ensight_cell_names[14] = "penta15";

    // assign values to vrtx_center
    vrtx_cnt[0]  = 1;
    vrtx_cnt[1]  = 2;
    vrtx_cnt[2]  = 3;
    vrtx_cnt[3]  = 3;
    vrtx_cnt[4]  = 6;
    vrtx_cnt[5]  = 4;
    vrtx_cnt[6]  = 8;
    vrtx_cnt[7]  = 4;
    vrtx_cnt[8]  = 10;
    vrtx_cnt[9]  = 5;
    vrtx_cnt[10] = 13;
    vrtx_cnt[11] = 8;
    vrtx_cnt[12] = 20;
    vrtx_cnt[13] = 6;
    vrtx_cnt[14] = 15;

    // assign values to cell_type_index
    cell_type_index[0]  = point;
    cell_type_index[1]  = two_node_bar;
    cell_type_index[2]  = three_node_bar;
    cell_type_index[3]  = three_node_triangle;
    cell_type_index[4]  = six_node_triangle;
    cell_type_index[5]  = four_node_quadrangle;
    cell_type_index[6]  = eight_node_quadrangle;
    cell_type_index[7]  = four_node_tetrahedron;
    cell_type_index[8]  = ten_node_tetrahedron;
    cell_type_index[9]  = five_node_pyramid;
    cell_type_index[10] = thirteen_node_pyramid;
    cell_type_index[11] = eight_node_hexahedron;
    cell_type_index[12] = twenty_node_hexahedron;
    cell_type_index[13] = six_node_wedge;
    cell_type_index[14] = fifteen_node_wedge;

    // calculate file prefixes and mkdir
   
    // create the ensight directory name (ens_prefix)
    ens_prefix = gd_wpath + "/" + prefix + "_ensight"; 

    // remove old ensight directory
    std::ostringstream rm_ensight;
    rm_ensight << "rm -rf " << ens_prefix;
    system(rm_ensight.str().c_str());

    // build the ensight directory
    mkdir(ens_prefix.c_str(), ENSIGHT_DIR_MODE);
    if (errno != 0)
    {
	std::ostringstream dir_error;
	dir_error << "Error opening ensight directory: " << strerror(errno);
	Insist (0,  dir_error.str().c_str());
    }

    // fill cell centered and vertex data names from input fields
    int nens_vdata = ens_vdata_names.size();
    int nens_cdata = ens_cdata_names.size();
    for (int i = 0; i < ens_vdata_names_in.size(); i++)
	ens_vdata_names[i] = ens_vdata_names_in[i];
    for (int i = 0; i < ens_cdata_names_in.size(); i++)
	ens_cdata_names[i] = ens_cdata_names_in[i];
       
    // CREATE THE VARIABLE NAMES, AND CHECK TO MAKE SURE THE 
    // VARIABLE NAMES ARE UNIQUE AND OF LENGTH GREATER THAN ZERO	      
    
    // here we check for size and spaces in names; however the following
    // characters are not valid for enight names and should be not be in the
    // names even though we don't check for them (yet):
    //    ( ) [ ] + - @ ! # * ^ $ / space

    // check for blanks in names
    for (int i = 0; i < nens_vdata; i++)
    {
	std_string test = ens_vdata_names[i];
	
	if (test.size() == 0 || test.size() > 19)
	    Insist (0, "Variable name inappropriately sized!");

	if (test.find(' ') != std_string::npos) 
	    Insist (0, "Spaces found in the vertex data name!");

	for (int j = i+1; j < nens_vdata; j++)
	    if (ens_vdata_names[j] == test) 
		Insist (0, "Non-unique names found in vdata!");
    }
    for (int i = 0; i < nens_cdata; i++)
    {
	std_string test = ens_cdata_names[i];
	
	if (test.size() == 0 || test.size() > 19)
	    Insist (0, "Variable name inappropriately sized!");

	if (test.find(' ') != std_string::npos) 
	    Insist (0, "Spaces found in the cell data name!");

	for (int j = i+1; j < nens_cdata; j++)
	    if (ens_cdata_names[j] == test) 
		Insist (0, "Non-unique names found in cdata!");
    }

    // calculate case file filename
    case_filename = ens_prefix + "/" + prefix + ".case";

    // calculate and make the geometry directory
    geo_dir = ens_prefix + "/geo";
    mkdir(geo_dir.c_str(), ENSIGHT_DIR_MODE);

    // make data directory names and directories
    Check (ens_vdata_names.size() == vdata_dirs.size());
    Check (ens_cdata_names.size() == cdata_dirs.size());
    for (int i = 0; i < ens_vdata_names.size(); i++)
    {
	vdata_dirs[i] = ens_prefix + "/" + ens_vdata_names[i];
	mkdir(vdata_dirs[i].c_str(), ENSIGHT_DIR_MODE);
    }
    for (int i = 0; i < ens_cdata_names.size(); i++)
    {
	cdata_dirs[i] = ens_prefix + "/" + ens_cdata_names[i];
	mkdir(cdata_dirs[i].c_str(), ENSIGHT_DIR_MODE);
    }   
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
 * \param ipar_in IVF field of pointers to cell vertices.  ipar_in is
 * dimensioned [0:ncells-1, 0:nvertices/cell-1].  Given ipar[i][j], is the
 * jth+1 vertex number for the ith+1 cell.  The vertex number is the i+1
 * entry from the pt_coor_in field.
 *
 * \param iel_type ISF field of Ensight_Cell_Types.  Each cell in the problem
 * must be associated with a Ensight_Cell_Types enumeration object.
 *
 * \param cell_rgn_index ISF field of region identifiers for each cell.  This
 * matches a region index to each cell in the problem.
 *
 * \param pt_coor_in FVF field of vertex coordinates. pt_coor_in is
 * dimensioned [0:nvertices-1, 0:ndim-1].  For each vertex point give the
 * value in the appropriate dimension.
 *
 * \param ens_vrtx_data_in FVF field of vertex data.  ens_vrtx_data_in is
 * dimensioned [0:nvertices-1, 0:number of vertex data fields - 1].  The
 * ordering of the second index must match the ens_vdata_names_in field input
 * argument to Ensight_Translator::Ensight_Translator().  The ordering of the
 * first index must match the vertex ordering from pt_coor_in.
 *
 * \param ens_cell_data_in FVF field of cell data.  ens_cell_data_in is
 * dimensioned [0:ncells-1, 0:number of cell data fields - 1].  The ordering
 * of the second index must match the ens_cdata_names_in field input argument
 * to Ensight_Translator::Ensight_Translator().  The ordering of the first
 * index must match the cell ordering from ipar_in.
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
void Ensight_Translator::ensight_dump(int icycle,
				      double time,
				      double dt, 
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

    // Increment dump counter adn add dump time
    igrdump_num++;
    dump_times.push_back(time);
    Check (igrdump_num == dump_times.size());
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
	      << " time= " << time << " dt= " << dt << endl
	      << "dir= " << ens_prefix << ", dump_number= " 
	      << igrdump_num << endl;

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
    int index;
    
    for (int i = 0; i < nparts; i++)
    {
	find_location_c = find(rgn_numbers.begin(), rgn_numbers.end(),
			       parts_list[i]);

	if (find_location_c != rgn_numbers.end())
	{
	    index = find_location_c - rgn_numbers.begin();
	    part_names.push_back(rgn_name[index]);
	}
	else if (find_location_c == rgn_numbers.end())
	{
	    Insist (0, "Didn't supply a region name!");
	}
    }
    Insist (parts_list.size() == part_names.size(), "Mismatch on part size!");
    Insist (rgn_name.size() == parts_list.size(), "Mismatch on region size!");

    // CREATE THE INDEXES FOR DUMPING CELL DATA.  ELEMENT (CELL) DATA IS
    // DUMPED BY PART NUMBER AS THE PRIMARY INDEX, AND BY CELL-TYPE AS THE
    // SECONDARY INDEX.  INDEX_CELL CONTAINS THE ORDER IN WHICH TO DUMP THE
    // CELL DATA.  IPTR_INDEX_CELL(i+(j-1)*ncell_types_ensight) POINTS TO THE
    // ELEMENT OF INDEX_CELL WHICH IS THE FIRST ELEMENT OF CELL-TYPE i, AND
    // PART j. 

    // cell index and integer pointer into index_cell
    int kk = num_ensight_cell_types * nparts + 1;
    vector<int> index_cell(ncells);
    vector<int> ptr_index_cell(kk, 0);

    {
	// define some temporary indices
	vector<int> itmp(ncells);
	vector<int> jtmp(kk, 0);

	for (int i = 0; i < ncells; i++)
	{
	    // find cell type
	    find_location = find(cell_type_index.begin(),
				 cell_type_index.end(), iel_type[i]);

	    if (find_location != cell_type_index.end())
		itmp[i] = find_location - cell_type_index.begin();
	    else
		Insist(0, "Element type not in cell_type_index!");

	    // find cell part number
	    find_location = find(parts_list.begin(), parts_list.end(),
				 cell_rgn_index[i]);

	    if (find_location != parts_list.end())
		itmp[i] += (find_location - parts_list.begin()) *
		    num_ensight_cell_types;
	    else
		Insist(0, "Part index not in parts_list!");

	    jtmp[itmp[i]] += 1;
	}
	
	for (int i = 1; i < kk; i++)
	{
	    ptr_index_cell[i] = ptr_index_cell[i-1] + jtmp[i-1];
	}
	for (int i = 0; i < kk; i++)
	    jtmp[i] = ptr_index_cell[i];

	for (int i = 0; i < ncells; i++)
	{
	    index_cell[jtmp[itmp[i]]] = i;
	    jtmp[itmp[i]] += 1;
	}
    }

    // >>> WRITE OUT DATA TO DIRECTORIES
    
    // write time to case file
    ensight_case(time);

    // WRITE THE GEOMETRY FILE
    ensight_geom(ens_postfix, icycle, time, dt, ipar, pt_coor,
		 part_names, index_cell, ptr_index_cell);

    // write the vertex data
    if (ens_vrtx_data.nrows() > 0)
	ensight_vrtx_data(ens_postfix, ens_vrtx_data);

    // write out the cell data
    if (ens_cell_data.nrows() > 0)
	ensight_cell_data(ens_postfix, ens_cell_data, index_cell,
			  ptr_index_cell, part_names); 
}

//---------------------------------------------------------------------------//
// ENSIGHT DATA OUTPUT FUNCTIONS (PRIVATE)
//---------------------------------------------------------------------------//
/*!
 * \brief Write out data to ensight geometry file.
 */
template<class IVF, class FVF> void 
Ensight_Translator::ensight_geom(const std_string &ens_postfix,
				 const int icycle,
				 const double time,
				 const double dt,
				 const rtt_traits::Viz_Traits<IVF> &ipar, 
				 const rtt_traits::Viz_Traits<FVF> &pt_coor,
				 const sf_string &part_names,
				 const sf_int &index_cell, 
				 const sf_int &ptr_index_cell)
{
    using std::ofstream;
    using std::setw;
    using std::ios;
    using std::endl;
    using std::string;

    Require (ipar.nrows() == index_cell.size());

    // make output file for this timestep
    string filename  = geo_dir + "/" + ens_postfix;
    const char *file = filename.c_str();
    ofstream geomout(file);

    // write the header
    geomout << "Description line 1" << endl;
    geomout << "probtime " << time << " cycleno " << icycle << endl;
    geomout << "node id given" << endl;
    geomout << "element id given" << endl;
    geomout << "coordinates" << endl;
    geomout << pt_coor.nrows() << endl;

    // WRITE THE UNSTRUCTURED COORDINATE DATA. DATA IS ALWAYS
    // 3D, 2D DATA IS CONFINED TO A PLANE, 1D TO A LINE IN 3-SPACE
    
    int ndim = pt_coor.ncols(0);

    for (int i = 0; i < pt_coor.nrows(); i++)
    {
	Check (pt_coor.ncols(i) == ndim);

	// print out node index: REMEMBER, node indices run from [1,N]
	// whereas the pt_coor array runs from [0,N-1]
	geomout << setw(8) << i+1;
	
	// set precision
	geomout.precision(5);
	geomout.setf(ios::scientific, ios::floatfield);

	// do 3-D output
	if (ndim == 3)
	    geomout << setw(12) << pt_coor(i, 0) << setw(12) << pt_coor(i, 1)
		    << setw(12) << pt_coor(i, 2) << endl;

	// do 2-D output
	if (ndim == 2)
	{
	    double zero = 0.0;
	    geomout << setw(12) << pt_coor(i, 0) << setw(12) << pt_coor(i, 1) 
		    << setw(12) << zero << endl;
	}
	   
	// do 1-D output
	if (ndim == 1)
	{
	    double zero = 0.0;
	    geomout << setw(12) << pt_coor(i, 0) << setw(12) << zero
		    << setw(12) << zero << endl;
	}
    }

    // write the cell data
    int counter  = 0;
    int num_elem = 0;
    int high;
    int low;
    for (int ipart = 0; ipart < part_names.size(); ipart++)
    {
	// output part number and names
	geomout << "part " << ipart+1 << endl;
	geomout << part_names[ipart] << endl;

	// find ensight cell type
	for (int type = 0; type < num_ensight_cell_types; type++)
	{
	    // calculate number of elements of type
	    low  = ptr_index_cell[counter];
	    high = ptr_index_cell[counter+1];
	    num_elem = high - low;
	    
	    if (num_elem > 0)
	    {
		geomout << ensight_cell_names[type] << endl;
		geomout << setw(8) << num_elem << endl;
		
		for (int i = low; i < high; i++)
		{
		    geomout << setw(8) << index_cell[i] + 1;
		    for (int j = 0; j < vrtx_cnt[type]; j++)
			geomout << setw(8) << ipar(index_cell[i], j);
		    geomout << endl;
		}
	    }

	    // increment type counter in parts list
	    counter++;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write out data to ensight vertex data.
 */
template<class FVF> void Ensight_Translator::ensight_vrtx_data
(const std_string &ens_postfix,
 const rtt_traits::Viz_Traits<FVF> &ens_vrtx_data)
{
    using std::ofstream;
    using std::endl;
    using std::ios;
    using std::setw;
    using std::string;

    // loop over all vertex data fields and write out data for each field
    for (int nvd = 0; nvd < ens_vrtx_data.ncols(0); nvd++)
    {
	// open file for this data
	string filename  = vdata_dirs[nvd] + "/" + ens_postfix;
	const char *file = filename.c_str();
	ofstream vout(file);

	vout << ens_vdata_names[nvd] << endl;
	
	// set precisions
	vout.precision(5);
	vout.setf(ios::scientific, ios::floatfield);
	
	// write 6 data entries per line
	int count   = 0;
	int columns = 0;

	if (ens_vrtx_data.nrows() > 6)
	    columns = 6;
	else 
	    columns = ens_vrtx_data.nrows();

	while (count < ens_vrtx_data.nrows())
	{
	    for (int j = 0; j < columns; j++)
	    {
		vout << setw(12) << ens_vrtx_data(count, nvd);
		count++;
	    }
	    vout << endl;

	    if (ens_vrtx_data.nrows() - count < columns) 
		columns = ens_vrtx_data.nrows() - count;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write out data to ensight cell data.
 */
template<class FVF> void Ensight_Translator::ensight_cell_data
(const std_string &ens_postfix,
 const rtt_traits::Viz_Traits<FVF> &ens_cell_data,
 const sf_int &index_cell,
 const sf_int &ptr_index_cell,
 const sf_string &part_names)
{
    using std::ofstream;
    using std::endl;
    using std::ios;
    using std::setw;
    using std::string;

    // loop over all vertex data fields and write out data for each field
    for (int ncd = 0; ncd < ens_cell_data.ncols(0); ncd++)
    {
	// open file for this data
	string filename  = cdata_dirs[ncd] + "/" + ens_postfix;
	const char *file = filename.c_str();
	ofstream cellout(file);

	cellout << ens_cdata_names[ncd] << endl;

	// set precision
	cellout.precision(5);
	cellout.setf(ios::scientific, ios::floatfield);

	// write the cell data
	int counter  = 0;
	int num_elem = 0;
	int high;
	int low;

	// count runs from 0 to ncells-1 across types
	int cell_count = 0; 

	for (int ipart = 0; ipart < part_names.size(); ipart++)
	{
	    // output part number
	    cellout << "part " << ipart+1 << endl;
	    
	    // loop over ensight cell types
	    for (int type = 0; type < num_ensight_cell_types; type++)
	    {
		// calculate number of elements of type
		low  = ptr_index_cell[counter];
		high = ptr_index_cell[counter+1];
		num_elem = high - low;

		// print out if number of elements is greater than zero
		if (num_elem > 0)
		{
		    int columns    = 0;
		    int elem_count = 0;

		    // determine groups of six for printing
		    if (num_elem > 6) 
			columns = 6;
		    else
			columns = num_elem;

		    // printout cell-type name
		    cellout << ensight_cell_names[type] << endl;

		    // print out data
		    for (int i = low; i < high; i++)
		    {
			while (elem_count < num_elem)
			{
			    for (int j = 0; j < columns; j++)
			    {
				cellout << setw(12) <<
				    ens_cell_data(index_cell[cell_count],
						  ncd); 
				cell_count++;
				elem_count++;
			    }
			    cellout << endl;

			    if (num_elem - elem_count < columns)
				columns = num_elem - elem_count;
			}
		    }
		}
		    
		// increment counter of death
		counter++;		
	    } 
	}
    }
}

} // end of rtt_viz

//---------------------------------------------------------------------------//
//                        end of viz/Ensight_Translator.t.hh
//---------------------------------------------------------------------------//
