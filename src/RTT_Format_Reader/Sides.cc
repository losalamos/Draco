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
/*!
 * \brief Compares two integer vectors for equality. 
 * \param low_val Integer vector that is assumed to be the lower value.
 * \param high_val Integer vector that is assumed to be the higher value.
 * \return high_val > low_val.
 */
bool Sides::compareVectorInt(vector_int low_val, vector_int high_val)
{
    // compare two integer vectors to determine which contains the lowest 
    // numbers. First sort the vectors into ascending order.
    std::sort(low_val.begin(),low_val.end());
    std::sort(high_val.begin(),high_val.end());

    // compare the vectors
    return (low_val < high_val);
}
/*!
 * \brief Renumbers the sides in ascending order based upon their node numbers.
 */
void Sides::sortData()
{
    vector_vector_int sort_vector(dims.get_nsides(),1);
    vector_int original(1);
    vector_int temp_sideType = sideType ;
    vector_vector_int temp_flags = flags;
    sort_map.resize(dims.get_nsides());

    for (int i = 0; i < dims.get_nsides(); ++i)
    {
        sort_vector[i].resize(cellDefs.get_nnodes(sideType[i]));
	for (int j = 0; j < cellDefs.get_nnodes(sideType[i]); ++j)
	    // map the user-input node numbers to the sorted node numbers.
	    nodes[i][j] = nodesClass.get_map(nodes[i][j]);
	sort_vector[i] = nodes[i];
    }
    std::sort(sort_vector.begin(),sort_vector.end(),Sides::compareVectorInt);

    // establish the mapping between the old and new side numbers, and assign
    // the nodes, sidetypes, and flags with the new numbering.
    for (int i = 0; i < dims.get_nsides(); ++i)
    {
        if (original.size() != cellDefs.get_nnodes(sideType[i]))
	    original.resize(cellDefs.get_nnodes(sideType[i]));
	original = nodes[i];

        bool mapped = false;
        int low_index  = 0;
        int high_index = dims.get_nsides() - 1;
        int k = (high_index + low_index) / 2;
	while (!mapped)
	{
	    // Correlate the original and sorted side numbers with a binary 
	    // search
            while ((high_index - low_index) > 1)
            {
	        if (compareVectorInt(original,sort_vector[k]))
	            high_index = k;
	        else
	            low_index  = k;
                k = (high_index + low_index) / 2;
            }
	    Insist(k < dims.get_nsides(), "Overflow in side sort while loop!");

	    k = low_index;
	    // Check to make sure these sides have the same nodes
	    if (original == sort_vector[k])
	    {
	        mapped = true;
		sort_map[i] = k;
		nodes[i] = sort_vector[i];
		sideType[k] = temp_sideType[i];
		flags[k] = temp_flags[i];
	    }
	    else if (k < high_index)
	        ++low_index;
	    else
	        Insist(0, "Error in side sort loop!");
	}
    }

    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;

    if (debugging)
    {
        for (int i = 0; i < dims.get_nsides(); ++i)
	{
	    std::cout << i << " " << sort_map[i] << " " << sideType[i] << " ";
	    for (int j = 0; j < cellDefs.get_nnodes(sideType[i]); ++j)
	        std::cout << nodes[i][j] << " ";
	    for (int f = 0; f < dims.get_nside_flag_types(); ++f)
	        std::cout << flags[i][f];
	    std::cout << std::endl;
	}
    }
    // free memory
    sort_vector.resize(0);
    original.resize(0);
    temp_sideType.resize(0);
    temp_flags.resize(0);
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                     end of RTT_Format_Reader/Sides.cc
//---------------------------------------------------------------------------//
