//----------------------------------*-C++-*----------------------------------//
// Connect.hh
// B.T. Adams
// 7 June 00
/*! 
 * \file   meshReaders_Services/Connect.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for meshReaders_Services/Connect class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __meshReaders_Services_Connect_hh__
#define __meshReaders_Services_Connect_hh__

#include "ds++/isSorted.hh"
#include "ds++/SortPermutation.hh"
#include "meshReaders/RTT_Mesh_Reader.hh"
#include <vector>
#include <set>
#include <map>

namespace rtt_meshReaders_Services
{
/*!
 * \brief Determines the mesh connectivity from the input mesh file data.
 */
class Connect
{
    // typedefs
    typedef std::set<int> set_int;
    typedef std::vector<int> vector_int;
    typedef std::vector<std::vector<int> > vector_vector_int;
    typedef std::vector<double> vector_dbl;
    typedef std::vector<std::vector<double> > vector_vector_dbl;
    typedef rtt_meshReaders::RTT_Mesh_Reader RTT_Mesh_Reader;
    typedef rtt_meshReaders::Element_Definition Element_Definition;
    typedef Element_Definition::Element_Type Element_Type;

    // data needed for all cells.
    // Mesh reader data
    rtt_dsxx::SP<RTT_Mesh_Reader> mesh_reader;
    // Cell Definitions data
    std::map<Element_Type, rtt_dsxx::SP<Element_Definition> > elem_defs;
    std::map<std::string, std::set<int> > elem_sets;
    // Nodes data
    vector_vector_dbl nodes_coords;
    int nnodes;
    // Sides data
    std::vector<std::string> bndry_flags;
    vector_int flag_numbs;
    std::vector<Element_Type> sides_types;
    vector_vector_int sides_nodes;
    vector_int sides_bndry;
    int nsides;
    // Cells data
     std::vector<Element_Type> cells_types;
    vector_vector_int cells_nodes;
    int ncells;
    // data needed specific to mutual faces.
    std::multimap<set_int, int> sideNodes;
    std::multimap<set_int, int> cellFaceNodes;
    // data needed specific to hanging nodes.
    std::multimap<int, vector_int > cmplxNodes;
    std::multimap<vector_int, vector_int > cmplxCells;
    std::multimap<double, vector_int, std::greater<double> > circum;
    vector_int nodes_map;

    // data derived from the calculation.
    std::vector<std::vector<std::vector<int> > > adjCell;
    std::multimap<int, int> bndryFaces;
    std::multimap<int, vector_int > Side_to_Cell_Face;

  public:
    Connect() { /* empty */ } 
/*!
 * \brief Constructs a Connect class object which determines the 
 *        mesh connectivity via a call to a private memember function.
 * \param mesh_reader_ meshReaders class object.
 * \param bndry_flags_ Flags used to indicate boundary conditions.
 * \param flag_numbs_ Flag numbers associated with the bndry_flags.
 * \param comp Functor used to compare node coordinates for sorting if 
 *             hanging nodes exist.
 */
    template<class StrictWeakOrdering>
    Connect(rtt_dsxx::SP<RTT_Mesh_Reader> mesh_reader_,
		       const std::vector<std::string> & bndry_flags_,
		       const std::vector<int> & flag_numbs_,
		       const StrictWeakOrdering & comp) 
        : mesh_reader(mesh_reader_), bndry_flags(bndry_flags_), 
	  flag_numbs(flag_numbs_)
{
    organizeData(mesh_reader);
    connectMutualFaces();
    if (!cmplxCells.empty())
        connectHangingNodes(comp);
}
    // Defaulted destructor.
    ~Connect() {}

  private:
    void organizeData(rtt_dsxx::SP<RTT_Mesh_Reader> mesh_reader);
    std::multimap<set_int, int> generateSideNodes();
    std::multimap<set_int, int> generateCellFaceNodes();
    void treatComplexFace(const vector_int & cell_faces);
    void treatSimpleFace(const int & cell, const int & faceNum,
			 std::multimap<set_int, int>::iterator & mmiter);
    void connectMutualFaces();
    template<class StrictWeakOrdering>
    int get_node(std::vector<double> node_coords, 
		 const StrictWeakOrdering & comp) const;
    template<class StrictWeakOrdering>
    void connectHangingNodes(const StrictWeakOrdering & comp);

  public:
    void printAdjacentCells() const;

    int get_adjCell_size(int cell, int face) const
    { return adjCell[cell][face].size(); }
/*!
 * \brief Returns the adjacent cell array.
 * \return The adjacent cell number array.
 */
    std::vector<std::vector<std::vector<int> > > get_adjCell() const
    { return adjCell; }
/*!
 * \brief Returns the number of the cell adjacent to the specified cell, face,
 *        and optional adjacent cell index.
 * \param cell Cell number.
 * \param face Face number.
 * \param adjcell Adjacent cell number (defaults to 0).
 * \return The adjacent cell number.
 */
    int get_adjCell(int cell, int face, int adjcell = 0) const
    { return adjCell[cell][face][adjcell]; }
/*!
 * \brief Returns the number of boundary faces (i.e., faces that are either
 *        on the outer boundary of the problem geometry or a connection between
 *        cells with different refinement levels in an AMR mesh) with the 
 *        specified face number
 * \param face Face number.
 * \return The number of boundary faces.
 */
    int get_bndryFaces_count(int face) const { return bndryFaces.count(face); }

    set_int get_bndryCells(int face) const;
    bool check_bndryFace(int cell, int face) const;
    int get_Cell_from_Side(int side) const;
    int get_Cell_Face_from_Side(int side) const;
};
/*!
 * \brief Determines a node number based upon the specified coordinate values.
 * \param node_coords Coordinate values.
 * \param comp StrictWeakOrdering used in the comparison.
 * \return The node number.
 */
template<class StrictWeakOrdering>
int Connect::get_node(std::vector<double> coords, 
		      const StrictWeakOrdering & comp) const
{
    const double EPSILON = 1.0e-06;
    int node_number = 0;
    std::vector<double> node_coords;
    bool found = false;
    int dim = coords.size();

    // Find the desired coordinates with a binary search.
    int low_index  = 0;
    int high_index = nnodes - 1;
    node_number = (high_index + low_index) / 2;
    while ((high_index - low_index) > 1)
    {
        node_coords = nodes_coords[nodes_map[node_number]];
        if (comp(coords, node_coords))
            high_index = node_number;
        else
            low_index  = node_number;
        node_number = (high_index + low_index) / 2;
    }
    Insist(node_number < nnodes, 
           "Overflow in get_node binary search routine!");
    node_number = low_index;

    // A full linear search if the nodes are not sorted or just a couple
    // of nodes otherwise.
    vector_dbl epsilon(dim);
    while (node_number < nnodes && !found)
    {
	int true_count = 0;
        for (int d = 0; d < dim; d++)
	{
            if (coords[d] != 0 && nodes_coords[nodes_map[node_number]][d] != 0)
	        epsilon[d] = EPSILON*((std::fabs(coords[d]) + 
		    std::fabs(nodes_coords[nodes_map[node_number]][d]))/2.);
	    else
	        epsilon[d] = EPSILON;

	    if (std::fabs(nodes_coords[nodes_map[node_number]][d] - coords[d]) 
		<= epsilon[d])
	        ++true_count;
	}
        if (true_count == dim)
	    found = true;
        else
      	    ++node_number;
    }
    Insist(node_number < nnodes,
     	   "Node number could not be found from its coordinates!");
    return nodes_map[node_number];
}
/*!
 * \brief Determines the mesh connectivity for cells that have hanging nodes.
 * \param comp Functor used to compare node coordinates for sorting.
 */
template<class StrictWeakOrdering>
void Connect::connectHangingNodes(const StrictWeakOrdering & comp)
{
    vector_int cell_faces(2);

    // Done with all of the structured cells. This essentially leaves a three-
    // dimensional jigsaw puzzle that is composed of blocks formed by cells 
    // that share a common geometry. The unassigned faces represent the
    // interlocks between these various pieces. Some of these faces will be 
    // "whole" in that groups of cells need to be combined into "subfaces" to
    // form the corresponding mating face. The whole faces can be distinguished
    // from the subfaces by size. The circum multimap contains the cell face 
    // circumferences in descending order. For the case of a continuous AMR
    // mesh, the center of each face corresponds to the location of a single 
    // node that is common to each of the subgrouped cells (and this node will
    // not occur anywhere else). Use the data in the cmplxCells  multimap to 
    // calculate the coordinates of the center of the face, find the node
    // that lives there, and assign each of the cell faces that contain this
    // node (as identified in the cmplxNode multimap) to this whole face. Then
    // remove all of these cell faces from the cmplxCells multimap. Note that 
    // the trick of locating the grouped subfaces by finding the whole face
    // center is unique to the continuous adaptive refinement mesh. Another
    // technique would probably be required for another mesh topology. The fact
    // that these "junction nodes" appear less frequently (because they do not
    // actually exist on the whole face) might be of use. The fact that the
    // "junction lines" between two grouped subfaces appear exactly twice could
    // also be valuable.

    // Create a mapping of the nodes based upon their coordinates to allow a 
    // binary search in the get_node function.
    nodes_map.resize(nodes_coords.size());
    if (rtt_dsxx::isSorted(nodes_coords.begin(), nodes_coords.end(), 
			   comp))
    {
        for (int n = 0; n < nodes_coords.size(); n++)
	    nodes_map[n] = n;
    }
    else
    {
        rtt_dsxx::SortPermutation nodeSorter(nodes_coords.begin(), 
					     nodes_coords.end(), comp);
	for (int n = 0; n < nodes_coords.size(); n++)
	    nodes_map[nodeSorter.inv(n)] = n;
    }

    std::multimap<double, vector_int, std::greater<double> >::iterator cfiter= 
        circum.begin();
    // this is a no-op if the mesh is fully structured - a countdown otherwise.
    while (!cmplxCells.empty())
    {
        // skip over grouped cell faces in circum that have already been 
        // assigned to whole faces and erased from the cmplxCells multimap.
	while (cmplxCells.count(cfiter->second) == 0)
	    ++cfiter;

	cell_faces = cfiter->second;	    
	int cell = cell_faces[0];
	int faceNum = cell_faces[1];

	// calculate the coordinates of the face center point.
	vector_int & faceNodes = cmplxCells.find(cell_faces)->second;
	std::vector<double> center(nodes_coords[0].size());
	for (int n = 0; n < faceNodes.size(); n++)
	{
	    for (int d = 0; d < nodes_coords[0].size(); d++)
	        center[d] += nodes_coords[faceNodes[n]][d];
	}
	for (int d = 0; d < nodes_coords[0].size(); d++)
	    center[d] /= static_cast<double>(faceNodes.size());

	// find a node with these coordinates, and set up a corresponding 
	// iterator in the cmplxNodes multimap.
	int centerNode = get_node(center, comp);
	std::multimap<int,vector_int >::iterator fiter = 
	    cmplxNodes.find(centerNode);

	// assign all of the cell faces that contain this node to this whole
	// face, and assign this whole face to each of the subfaces.
	while (fiter !=  cmplxNodes.end() && fiter->first == centerNode)
	{
	    vector_int other_cell_faces = fiter->second;
	    int otherCell = other_cell_faces[0];
	    int otherFaceNum = other_cell_faces[1];

	    adjCell[cell][faceNum].push_back(otherCell);
	    adjCell[otherCell][otherFaceNum].push_back(cell);
	    cmplxCells.erase(other_cell_faces);
	    ++fiter;
	}
	cmplxCells.erase(cell_faces);
	++cfiter;
    }
    // clean house.
    circum.clear();
    cmplxNodes.clear();
    cell_faces.resize(0);
}

} // end namespace rtt_meshReaders_Services

#endif                      // __meshReaders_Services_Connect_hh__


//---------------------------------------------------------------------------//
//                end of meshReaders_Services/Connect.hh
//---------------------------------------------------------------------------//
