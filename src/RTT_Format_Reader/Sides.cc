//----------------------------------*-C++-*--------------------------------//
// Sides.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/Sides.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/Sides class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Sides.hh"

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Parses the sides block data from the mesh file via calls to 
 *        private member functions.
 * \param meshfile Mesh file name.
 */
void Sides::readSides(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the sides block keyword.
 * \param meshfile Mesh file name.
 */
void Sides::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "sides",
	   "Invalid mesh file: sides block missing");
    std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the sides block data.
 * \param meshfile Mesh file name.
 */
void Sides::readData(ifstream & meshfile)
{
    string dummyString;
    int sideNum;

    for (int i = 0; i < dims.get_nsides(); ++i)
    {
	meshfile >> sideNum;
	Insist(sideNum == i+1,
	       "Invalid mesh file: side index out of order");
	meshfile >> sideType[i];
	--sideType[i];
	Insist(dims.allowed_side_type(sideType[i]),
	       "Invalid mesh file: illegal side type");
	nodes[i].resize(cellDefs.get_nnodes(sideType[i]));
	for (int j = 0; j < cellDefs.get_nnodes(sideType[i]); ++j)
	{
	    meshfile >> nodes[i][j];
	    --nodes[i][j];
	}
	for (int j = 0; j < dims.get_nside_flag_types(); ++j)
	{
	    meshfile >> flags[i][j];
	    Insist(sideFlags.allowed_flag(j, flags[i][j]),
		   "Invalid mesh file: illegal side flag");
	}
	std::getline(meshfile, dummyString);
    }
}
/*!
 * \brief Reads and validates the end_sides block keyword.
 * \param meshfile Mesh file name.
 */
void Sides::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_sides",
	   "Invalid mesh file: sides block missing end");
    std::getline(meshfile, dummyString);       // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                     end of RTT_Format_Reader/Sides.cc
//---------------------------------------------------------------------------//
