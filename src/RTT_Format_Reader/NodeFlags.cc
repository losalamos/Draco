//----------------------------------*-C++-*--------------------------------//
// NodeFlags.hh
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/NodeFlags.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Format_Reader/NodeFlags class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "NodeFlags.hh"

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Parses the node_flags data block from the mesh file via calls to 
 *        private member functions.
 * \param meshfile Mesh file name.
 */
void NodeFlags::readNodeFlags(ifstream & meshfile)
{
    readKeyword(meshfile);
    readFlagTypes(meshfile);
    readEndKeyword(meshfile);
}
/*!
 * \brief Reads and validates the node_flags block keyword.
 * \param meshfile Mesh file name.
 */
void NodeFlags::readKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "node_flags",
	   "Invalid mesh file: node_flags block missing");
    std::getline(meshfile, dummyString);
}
/*!
 * \brief Reads and validates the node_flags block data.
 * \param meshfile Mesh file name.
 */
void NodeFlags::readFlagTypes(ifstream & meshfile)
{
    int flagTypeNum;
    string dummyString;

    for (int i = 0; i < dims.get_nnode_flag_types(); ++i)
    {
	meshfile >> flagTypeNum >> dummyString;
	// Ignore plurals in node flag definitions
	if (dummyString[dummyString.size()-1] == 's' || 
	    dummyString[dummyString.size()-1] == 'S')
	    dummyString.resize(dummyString.size()-1);
	Insist(flagTypeNum == i+1,
	       "Invalid mesh file: node flag type out of order");
	flagTypes[i] = new Flags(dims.get_nnode_flags(i), dummyString);
	std::getline(meshfile, dummyString);
	flagTypes[i]->readFlags(meshfile);
    }
}
/*!
 * \brief Reads and validates the end_node_flags block keyword.
 * \param meshfile Mesh file name.
 */
void NodeFlags::readEndKeyword(ifstream & meshfile)
{
    string dummyString;

    meshfile >> dummyString;
    Insist(dummyString == "end_node_flags",
	   "Invalid mesh file: node_flags block missing end");
    std::getline(meshfile, dummyString);       // read and discard blank line.
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                       end of RTT_Format_Reader/NodeFlags.cc
//---------------------------------------------------------------------------//
