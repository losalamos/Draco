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
/*!
 * \brief Compares real coordinate values for equality while considering
 *        machine precision. 
 * \param low_val (Assumed) lower value coordinates.
 * \param high_val (Assumed) higher value coordinates.
 * \return high_val > low_val && 
 *         abs(low_val- high_val) > epsilon (desired precision).
 */
bool Nodes::compareXYZ(const vector_dbl & low_val, const vector_dbl & high_val)
{
    // require agreement to six significant figures for equality. Note that
    // Shawn's tet mesh only agrees to four significant digits.
    const double EPSILON = 1.0e-06;
    vector_dbl epsilon(low_val.size());
    bool sorted = true;

    Insist(low_val.size() == high_val.size(),"Improper sort arguments!");
    int dim = low_val.size();
    for (int d = dim - 1; d >= 0; d--)
    {
        if (low_val[d] != 0 && high_val[d] != 0)
	    epsilon[d] = EPSILON * ((std::fabs(low_val[d]) + 
				     std::fabs(high_val[d]))/2.);
	else
	    epsilon[d] = EPSILON;
        // this strange looking logical operator will sort x,y,(z) coordinates
        // with x varying fastest, followed by y, and lastly z.
        if (high_val[d] < low_val[d] && 
	     std::fabs(low_val[d] - high_val[d]) > epsilon[d] && 
	       (d == dim-1 || 
	           (d == dim-2 && 
		     std::fabs(high_val[d+1]-low_val[d+1]) < epsilon[d+1]) ||
	                 (std::fabs(high_val[d+1]-low_val[d+1]) < epsilon[d+1]
			  &&
		          std::fabs(high_val[d+2]-low_val[d+2])<epsilon[d+2])))
	{
	    sorted = false;
	    d = -1;
	}
    }
    return sorted;
}
/*!
 * \brief Renumbers the nodes in ascending order based upon their coordinates
 *        (x, y, and then z).
 */
void Nodes::sortData()
{
    vector_vector_dbl sort_vector = coords;
    vector_dbl original(dims.get_ndim());
    vector_int temp_parents = parents;
    vector_vector_int temp_flags = flags;
    sort_map.resize(dims.get_nnodes());

    std::sort(sort_vector.begin(), sort_vector.end(), Nodes::compareXYZ);

    // establish the mapping between the old and new node numbers, and assign
    // the coordinates, parents, and flags with the new numbering.
    for (int i = 0; i < dims.get_nnodes(); ++i)
    {
	original = coords[i];
        bool mapped = false;
        int low_index  = 0;
        int high_index = dims.get_nnodes() - 1;
        int k = (high_index + low_index) / 2;
	while (!mapped)
	{
	    // Correlate the original and sorted node numbers/coordinates with
	    // a binary search
            while ((high_index - low_index) > 1)
            {
	        if (compareXYZ(original,sort_vector[k]))
	            high_index = k;
	        else
	            low_index  = k;
                k = (high_index + low_index) / 2;
            }
	    Insist(k < dims.get_nnodes(), "Overflow in node sort while loop!");

	    k = low_index;
	    // Check to make sure these nodes have the same coordinates
	    if (original == sort_vector[k])
	    {
	        mapped = true;
	        sort_map[i] = k;
		coords[i] = sort_vector[i];
		parents[k] = temp_parents[i];
		flags[k] = temp_flags[i];
	    }
	    else if (k < high_index)
	        ++low_index;
	    else
	        Insist(0, "Error in node sort loop!");
	}
    }
    // finally, now that the nodes are presumably sorted, update the parent
    // nodes to the new numbering system. Note that this step cannot be done
    // in the previous loop because the node and parent numbers are not 
    // necessarily the same.
    for (int i = 0; i < dims.get_nnodes(); ++i)
        parents[i] =  sort_map[parents[i]];

    // this flag is for debugging use only (not user-selectable).
    bool debugging = false;
    if (debugging)
    {
        for (int i = 0; i < dims.get_nnodes(); ++i)
	{
	  std::cout << i << " " << sort_map[i] << " " << parents[i] << " ";
	  for (int j = 0; j < dims.get_ndim(); ++j)
	      std::cout << coords[i][j] << " ";
	  for (int f = 0; f < dims.get_nnode_flag_types(); ++f)
	      std::cout << flags[i][f];
	  std::cout << std::endl;
	}
    }
    // free memory
    sort_vector.resize(0);
    original.resize(0);
    temp_parents.resize(0);
    temp_flags.resize(0);
}
/*!
 * \brief Determines a node number based upon the specified coordinate values.
 * \param node_coords Coordinate values.
 * \return The node number.
 */
int Nodes::get_node(vector_dbl node_coords) const
{
    const double EPSILON = 1.0e-06;
    int node_number = 0;
    bool found = false;
    int dim = dims.get_ndim();

    // Find the desired coordinates with a binary search if the nodes are
    // sorted, otherwise we are stuck with a linear search
    if (dims.get_renumber())
    {
        vector_dbl local_coords(dim);
        int low_index  = 0;
        int high_index = dims.get_nnodes() - 1;
        node_number = (high_index + low_index) / 2;
        while ((high_index - low_index) > 1)
        {
	    local_coords = get_coords(node_number);
            if (compareXYZ(node_coords,local_coords))
	        high_index = node_number;
	    else
	        low_index  = node_number;
            node_number = (high_index + low_index) / 2;
        }
        Insist(node_number < dims.get_nnodes(), 
	       "Overflow in get_node binary search routine!");
        node_number = low_index;
    }

    // A full linear search if the nodes are not sorted or just a couple
    // of nodes otherwise.
    vector_dbl epsilon(dim);
    while (node_number < dims.get_nnodes() && !found)
    {
	int true_count = 0;
        for (int d = 0; d < dim; d++)
	{
            if (node_coords[d] != 0 && get_coords(node_number,d) != 0)
	        epsilon[d] = EPSILON * ((std::fabs(node_coords[d]) + 
				         std::fabs(get_coords(node_number,d)))
					 /2.);
	    else
	        epsilon[d] = EPSILON;

	    if (std::fabs(get_coords(node_number,d) - node_coords[d]) <= 
		epsilon[d])
	        ++true_count;
	}
        if (true_count == dim)
	    found = true;
        else
      	    ++node_number;
    }
    Insist(node_number < dims.get_nnodes(),
     	   "Node number could not be found from its coordinates!");
    return node_number;
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                     end of RTT_Format_Reader/Nodes.cc
//---------------------------------------------------------------------------//
