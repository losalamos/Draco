//----------------------------------*-C++-*--------------------------------//
// CellFlags.hh
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/CellFlags.cc
 * \author B.T. Adams
 * \date   Mon Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/CellFlags class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "CellFlags.hh"

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Parses the cell_flags data block of the mesh file via calls to 
 *        private member functions.
 * \param meshfile Mesh file name.
 */
void CellFlags::readCellFlags(ifstream & meshfile)
{
    readKeyword(meshfile);
    readFlagTypes(meshfile);
    readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the cell_flags block keyword.
 * \param meshfile Mesh file name.
 */
void CellFlags::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "cell_flags",
	   "Invalid mesh file: cell_flags block missing");
    std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the cell_flags block data.
 * \param meshfile Mesh file name.
 */
void CellFlags::readFlagTypes(ifstream & meshfile)
{
    int flagTypeNum;
    string dummyString;

    for (int i = 0; i < dims.get_ncell_flag_types(); ++i)
    {
	meshfile >> flagTypeNum >> dummyString;
	// Ignore plurals in cell flag definitions
	if (dummyString[dummyString.size()-1] == 's' || 
	    dummyString[dummyString.size()-1] == 'S')
	    dummyString.resize(dummyString.size()-1);
	Insist(flagTypeNum == i+1,
	       "Invalid mesh file: cell flag type out of order");
	flagTypes[i] = new Flags(dims.get_ncell_flags(i), dummyString);
	std::getline(meshfile, dummyString);
	flagTypes[i]->readFlags(meshfile);
    }
}
/*!
 * \brief Reads and validates the end_cell_flags block keyword.
 * \param meshfile Mesh file name.
 */
void CellFlags::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_cell_flags",
	   "Invalid mesh file: cell_flags block missing end");
    std::getline(meshfile, dummyString);       // read and discard blank line.
}
/*!
 * \brief Returns the index to the cell flag type that contains the cell 
 *        materials.
 * \return The material cell flag type index.
 */
int CellFlags::get_material_flag_number() const 
{
    // Allow any combination of the phrase "material". 
    string matl("materilMATERIL");
    int material = -1;
    int length = 0;
    for (int f = 0; f < dims.get_ncell_flag_types(); f++)
    {
        string flag = flagTypes[f]->getFlagType();
        if ((flag[0] == 'm' || flag[0] == 'M') &&
      	    flag.find_first_not_of(matl) == string::npos &&  
       	    flag.find_first_not_of(matl) >= length)
       	{
	    length = flag.size();
      	    material = f;
       	}
    }
    return material;
}
/*!
 * \brief Returns the index to the optional cell flag type that contains the 
 *        cell volumetric sources.
 * \return The volumetric source cell flag type index.
 */
int CellFlags::get_volume_src_flag_number() const 
{
    // Allow any combination of the phrase "volume_source". 
    string source("volumesrcVOLUMESRC_");
    int vol_src = -1;
    int length = 0;
    for (int f = 0; f < dims.get_ncell_flag_types(); f++)
    {
        string flag = flagTypes[f]->getFlagType();
        if ((flag[0] == 'v' || flag[0] == 'V') &&
	    flag.find_first_not_of(source) == string::npos &&  
	    flag.find_first_not_of(source) >= length)
	{
	    length = flag.size();
	    vol_src = f;
	}
    }
    return vol_src;
}    
/*!
 * \brief Returns the index to the optional cell flag type that contains the 
 *        cell radiation sources.
 * \return The radiation source cell flag type index.
 */
int CellFlags::get_radiation_src_flag_number() const 
{
    // Allow any combination of the phrase "raditiation_source". 
    string source("raditonsuceRADITONSUCE_");
    int rad_src = -1;
    int length = 0;
    for (int f = 0; f < dims.get_ncell_flag_types(); f++)
    {
        string flag = flagTypes[f]->getFlagType();
        if ((flag[0] == 'r' || flag[0] == 'R') &&
	    flag.find_first_not_of(source) == string::npos &&  
	    flag.find_first_not_of(source) >= length)
	{
	    length = flag.size();
	    rad_src = f;
	}
    }
    return rad_src;
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                   end of RTT_Format_Reader/CellFlags.cc
//---------------------------------------------------------------------------//
