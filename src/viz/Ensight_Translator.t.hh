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

//---------------------------------------------------------------------------//
// ENSIGHT DUMP PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Do an Ensight dump.
 *
 * Performs an ensight dump in the appropriate directory locations.
 */
template<class ISF, class IVF, class SSF, class FVF>
void Ensight_Translator::ensight_dump(const std_string &prefix,
				      int icycle,
				      double time,
				      double dt, 
				      const std_string &gd_wpath,
				      const IVF &ipar_in,
				      const ISF &iel_type,
				      const ISF &rgn_index,
				      const FVF &pt_coor_in,
				      const FVF &ens_vrtx_data_in,
				      const FVF &ens_cell_data_in,
				      const SSF &ens_vdata_names,
				      const SSF &ens_cdata_names,
				      const ISF &rgn_data,
				      const SSF &rgn_name,
				      const ISF &iproc)
{
    using rtt_traits::Viz_Traits;
    using std::ostringstream;
    using std::string;
    using std::vector;
    using std::find;

    // Increment dump counter
    igrdump_num++;

    // load traits for vector field types
    Viz_Traits<IVF> ipar(ipar_in);
    Viz_Traits<FVF> pt_coor(pt_coor_in);
    Viz_Traits<FVF> ens_vrtx_data(ens_vrtx_data_in);
    Viz_Traits<FVF> ens_cell_data(ens_cell_data_in);

    // define sizes used throughout
    int ncells     = ipar.nrows();
    int npoints    = pt_coor.nrows();
    int nens_vdata = ens_vdata_names.size();
    int nens_cdata = ens_cdata_names.size();
    int nrgn       = rgn_name.size();

    // Check sizes of all data.
    Require (iel_type.size() == ncells);
    Require (rgn_index.size() == ncells);
    Require (ens_cell_data.nrows() == ncells);
    Require (ens_vrtx_data.nrows() == npoints);
    Require (rgn_data.size() == nrgn);

    // total problem size
    int nens_data = nens_vdata + nens_cdata;

    // create the ensight directory name (ens_prefix)
    string ens_prefix;
    ens_prefix = gd_wpath + "/" + prefix + "_ensight"; 

    // build the ensight directory
    mkdir(ens_prefix.c_str(), ENSIGHT_DIR_MODE);

    // create ensight postfix indicators
    ostringstream ens_postfix_build;
    ens_postfix_build << "data." << igrdump_num;
    string ens_postfix = ens_postfix_build.str();

    // announce the graphics dump
    std::cout << "ENSIGHT GRAPHICS DUMP: icycle= " << icycle 
	      << " time= " << time << " dt= " << dt << endl
	      << "dir= " << ens_prefix << ", dump_number= " 
	      << igrdump_num << endl;

    // create the parts list
    vector<int>::const_iterator find_location_c;
    vector<int>::iterator       find_location;
    vector<int> parts_list;

    for (int i = 0; i < ncells; i++)
    {
	find_location = find(parts_list.begin(), parts_list.end(),
			     rgn_index[i]); 

	if (find_location == parts_list.end())
	    parts_list.push_back(rgn_index[i]);
    }
    
    // store the number of parts
    int nparts = parts_list.size();
    
    // create the parts names
    vector<string> part_names(nparts);
    int index;
    
    for (int i = 0; i < nparts; i++)
    {
	find_location_c = find(rgn_data.begin(), rgn_data.end(),
			       parts_list[i]);

	if (find_location_c != rgn_data.end())
	{
	    index = find_location_c - rgn_data.begin();
	    part_names.push_back(rgn_name[index]);
	}
	else if (find_location_c == rgn_data.end())
	{
	    Insist (0, "Didn't supply a region name!");
	}
    }

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
				 rgn_index[i]);

	    if (find_location != parts_list.end())
		itmp[i] += (find_location - parts_list.begin()) *
		    num_ensight_cell_types;
	    else
		Insist(0, "Part index not in parts_list!");

	    jtmp[itmp[i]] += 1;
	}

	for (int i = 1; i < kk; i++)
	    ptr_index_cell[i] = ptr_index_cell[i-1] + jtmp[i-1];
	for (int i = 0; i < kk; i++)
	    jtmp[i] = ptr_index_cell[i];

	for (int i = 0; i < ncells; i++)
	{
	    index_cell[jtmp[itmp[i]]] = i;
	    jtmp[itmp[i]] += 1;
	}
    }    

    // CREATE THE VARIABLE NAMES, AND CHECK TO MAKE SURE THE 
    // VARIABLE NAMES ARE UNIQUE AND OF LENGTH GREATER THAN ZERO	      
    
    // check for blanks in names
    for (int i = 0; i < nens_vdata; i++)
    {
	string test = ens_vdata_names[i];
	if (test.find(' ') != string::npos) 
	    Insist (0, "Spaces found in the vertex data name!");

	for (int j = i+1; j < nens_vdata; j++)
	    if (ens_vdata_names[j] == test) 
		Insist (0, "Non-unique names found in vdata!");
    }
    for (int i = 0; i < nens_cdata; i++)
    {
	string test = ens_cdata_names[i];
	if (test.find(' ') != string::npos) 
	    Insist (0, "Spaces found in the cell data name!");

	for (int j = i+1; j < nens_cdata; j++)
	    if (ens_cdata_names[j] == test) 
		Insist (0, "Non-unique names found in cdata!");
    }
}
//---------------------------------------------------------------------------//
//                        end of viz/Ensight_Translator.t.hh
//---------------------------------------------------------------------------//
