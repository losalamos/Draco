//----------------------------------*-C++-*--------------------------------//
// SideFlags.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/SideFlags.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/SideFlags class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "SideFlags.hh"

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Parses the side_flags data block of the mesh file via calls to 
 *        private member functions.
 * \param meshfile Mesh file name.
 */
void SideFlags::readSideFlags(ifstream & meshfile)
{
    readKeyword(meshfile);
    readFlagTypes(meshfile);
    readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the side_flags block keyword.
 * \param meshfile Mesh file name.
 */
void SideFlags::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "side_flags",
	   "Invalid mesh file: side_flags block missing");
    std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the side_flags block data.
 * \param meshfile Mesh file name.
 */
void SideFlags::readFlagTypes(ifstream & meshfile)
{
    int flagTypeNum;
    string dummyString;

    for (int i = 0; i < dims.get_nside_flag_types(); ++i)
    {
	meshfile >> flagTypeNum >> dummyString;
	// Ignore plurals in side flag definitions
	if (dummyString[dummyString.size()-1] == 's' || 
	    dummyString[dummyString.size()-1] == 'S')
	    dummyString.resize(dummyString.size()-1);
	Insist(flagTypeNum == i+1,
	       "Invalid mesh file: side flag type out of order");
	flagTypes[i] = new Flags(dims.get_nside_flags(i), dummyString);
	std::getline(meshfile, dummyString);
	flagTypes[i]->readFlags(meshfile);
    }
}
/*!
 * \brief Reads and validates the end_side_flags block keyword.
 * \param meshfile Mesh file name.
 */
void SideFlags::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_side_flags",
	   "Invalid mesh file: side_flags block missing end");
    std::getline(meshfile, dummyString);       // read and discard blank line.
}
/*!
 * \brief Returns the index to the side flag type that contains the problem 
 *        boundary conditions.
 * \return The boundary conditions side flag type index.
 */
int SideFlags::get_boundary_flag_number() const 
{
    // Allow any combination of the phrase "boundary_conditions". 
    string bc("boundarycditBOUNDARYCDIT_");
    int boundary = -1;
    int length = 0;
    for (int f = 0; f < dims.get_nside_flag_types(); f++)
    {
        string flag = flagTypes[f]->getFlagType();
        if ((flag[0] == 'b' || flag[0] == 'B') &&
       	    flag.find_first_not_of(bc) == string::npos &&  
       	    flag.find_first_not_of(bc) >= length)
       	{
	  length = flag.size();
	    boundary = f;
	}
    }
    return boundary;
}
/*!
 * \brief Returns the index to the optional side flag type that contains the 
 *        problem external sources.
 * \return The external source side flag type index.
 */
int SideFlags::get_surface_src_flag_number() const 
{
    // Allow any combination of the phrase "surface_source". 
    string surface("surfaceoSURFACEO_");
    int source = -1;
    int length = 0;
    for (int f = 0; f < dims.get_nside_flag_types(); f++)
    {
        string flag = flagTypes[f]->getFlagType();
	if ((flag[0] == 's' || flag[0] == 'S') &&
	    flag.find_first_not_of(surface) == string::npos &&  
	    flag.find_first_not_of(surface) >= length)
	{
	    length = flag.size();
	    source = f;
	}
    }
    return source;
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                 end of RTT_Format_Reader/SideFlags.cc
//---------------------------------------------------------------------------//
