//----------------------------------*-C++-*--------------------------------//
// Nodes.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/Nodes.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/Nodes class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Nodes.hh"

namespace rtt_RTT_Format_Reader
{
/*! 
 * \brief Parses the nodes block data from the mesh file via calls to private 
 *        member functions.
 * \param meshfile Mesh file name.
 */
void Nodes::readNodes(ifstream & meshfile)
{
    readKeyword(meshfile);
    readData(meshfile);
    readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the nodes block keyword.
 * \param meshfile Mesh file name.
 */
void Nodes::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "nodes",
	   "Invalid mesh file: nodes block missing");
    std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the nodes block data.
 * \param meshfile Mesh file name.
 */
void Nodes::readData(ifstream & meshfile)
{
    string dummyString;
    int nodeNum;

    for (int i = 0; i < dims.get_nnodes(); ++i)
    {
	meshfile >> nodeNum;
	Insist(nodeNum == i+1,
	       "Invalid mesh file: node index out of order");
	for (int j = 0; j < dims.get_ndim(); ++j)
	    meshfile >> coords[i][j];
	meshfile >> parents[i];
	--parents[i];
	for (int j = 0; j < dims.get_nnode_flag_types(); ++j)
	{
	    meshfile >> flags[i][j];
	    Insist(nodeFlags.allowed_flag(j, flags[i][j]),
		   "Invalid mesh file: illegal node flag");
	}
	std::getline(meshfile, dummyString);
    }
}
/*!
 * \brief Reads and validates the end_nodes block keyword.
 * \param meshfile Mesh file name.
 */
void Nodes::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_nodes",
	   "Invalid mesh file: nodes block missing end");
    std::getline(meshfile, dummyString);       // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                     end of RTT_Format_Reader/Nodes.cc
//---------------------------------------------------------------------------//
