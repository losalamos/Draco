//----------------------------------*-C++-*----------------------------------//
// Builder.hh
// B.T. Adams (bta@lanl.gov)
// 18 May 99
/*! 
 * \file   amr_mesh/Builder.hh
 * \author B.T. Adams
 * \date   Tue May 18 10:33:26 1999
 * \brief  Header file for CAR_CU_Builder class library.
 */
//---------------------------------------------------------------------------//
// @> CAR_CU_Builder class header file
//---------------------------------------------------------------------------//

#ifndef __amr_CAR_CU_Builder_hh__
#define __amr_CAR_CU_Builder_hh__

//===========================================================================//
// class CAR_CU_Builder - 
//
// Purpose : builds an CAR_CU_Mesh object
//
// revision history:
// -----------------
//  0) original (developed from OS_Builder.hh)
//===========================================================================//

#include "mc/Coord_sys.hh"
#include "Layout.hh"
#include "Mesh.hh"
#include "ds++/SP.hh"
#include "meshReaders/RTT_Format.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>

namespace rtt_amr 
{

// stl components
using std::vector;
using std::string;
using std::set;
using std::multimap;

// draco components
using rtt_dsxx::SP;

// RTT_Format input file class
using rtt_format::RTT_Format;

/*!
 * \brief  The continuous adaptive refinement (CAR) Cartesion unstructured 
 *         (CU) Builder class provides member functions that can be used to
 *         construct a CAR_CU_Mesh class object. The build_Mesh public member
 *         function takes a smart pointer to an rtt_format::RTT_Format class 
 *         object as its single input and returns a smart pointer to the new 
 *         CAR_CU_Mesh class object. This member function calls the private 
 *         member functions assign_Generations, build_Coord, and build_Layout
 *         to construct CAR_CU_Mesh::CCSF_i type and Coord_sys and Layout 
 *         class objects, respectively. A CAR_CU_Mesh::NCVF_d type vertex 
 *         array and a CAR_CU_Mesh::CCVF_i type cell_pair array are generated 
 *         directly. Each of these objects is required by the CAR_CU_Mesh 
 *         class constructor. The map_Nodes private member function is called 
 *         to eliminate superfluous node data that is present in the RTT_Format
 *         mesh file and the FC_Nodes private member function is called to add
 *         face-centered node data that is not provided in the file. The 
 *         seek_Adoption and the assign_Meshes private member functions are 
 *         also called to develop a hiearchy of generation-specific meshes 
 *         from the unstructured mesh data. This functionality was developed 
 *         to support proposed multigrid schemes based upon the CAR_CU_Mesh 
 *         class, but these features are not yet fully implemented herein.
 *
 *\sa The \ref amr_overview presents a summary of the capabilities and the
 *    intended usage of this mesh class.
 */     
class CAR_CU_Builder
{
private:
  // data from Parser needed to build mesh

  // coordinate system string
/*!
 * \brief The problem coordinate system (e.g., xyz).
 */
    string coord_system;
  // vertex (node) coordinate array
/*!
 * \brief The node coordinate values.
 */
    CAR_CU_Mesh::NCVF_d vertex;
  // cell vertexes
/*!
 * \brief The connections between cell faces.
 */
    CAR_CU_Mesh::CCVF_i cell_pair;

  // member functions for building CAR_CU_Mesh
  // build Layout helper functions
/*!
 * \brief Constructs a Layout class object.
 * \param coord An existing, initialized Coord_sys class object.
 * \param rttMesh Smart pointer to an existing, initialized RTT_Format
 *                class object.
 * \return Smart pointer to the new Layout class object.
 */
    SP<Layout> build_Layout(const Coord_sys & coord, 
			    const SP<RTT_Format> & rttMesh);

  // build Coord_sys helper functions
/*!
 * \brief Constructs a Coord_sys class object.
 * \return Smart pointer to the new Coord_sys class object.
 */
    SP<Coord_sys> build_Coord();

  // assign cells to generation levels in cell_gens and return a vector of 
  // sets of the cells present in each generation (gen_cells);
/*!
 * \brief Assigns cells to generation (i.e., refinement) levels.
 * \param rttMesh Smart pointer to an existing, initialized RTT_Format
 *                class object.
 * \param cell_gens Generation level of each cell (calculated).
 * \return Generation-specific cell sets.
 */
    vector<set<int> > assign_Generations(CAR_CU_Mesh::CCSF_i & cell_gens, 
					 const SP<RTT_Format> & rttMesh);

  // Establish a hiearchy for the mesh - modifies gen_cells directly and
  // returns a multimap of the new parents paired with their children.
/*!
 * \brief Generates new "parent cells" by grouping "child" (i.e., refined) 
 *        cells together. The process is repeated from the most refined
 *        generation level to the cells with no refinement. A newly formed
 *        parent cell can also be a child cell for a generation level of 
 *        lesser refinement.
 * \param gen_cells Generation-specific cell sets (returned by private member
 *                   function assign_Generations and modified directly by this
 *                   member function as new parent cells are created and 
 *                   added).
 * \param rttMesh Smart pointer to an existing, initialized RTT_Format
 *                class object.
 * \return Mapping between the parent cells and their children.
 */
    multimap<int, set<int> > seek_Adoption(vector<set< int> > & gen_cells, 
					   const SP<RTT_Format> & rttMesh);

  // Group the cells into meshes. Return vector of sets containing the cells
  // in each mesh. Can be used to assign child cells or to "connect" the
  // parents, depending upon whether the parents multimap is null or assigned,
  // respectively).
/*!
 * \brief Assigns either child or child and parent cells to 
 *        generation-specific meshes depending upon whether or not the 
 *        seek_Adoption private member function is called prior.
 * \param gen_cells Generation-specific cell sets (returned by private member
 *                   function assign_Generations and possibly modified by 
 *                   seek_Adoption).
 * \param parents Mapping between parent cells and their children (either 
 *                initialized on return from the private member function 
 *                seek_Adoption or null).
 * \param rttMesh Smart pointer to an existing, initialized RTT_Format
 *                class object.
 * \return Mapping between the parent cells and their children.
 */
    vector<set< int> > assign_Meshes(const vector<set< int> > & gen_cells, 
				     const multimap<int, set<int> > & parents,
				     const SP<RTT_Format> & rttMesh);

  // Renormalize node numbering to be specific to a mesh
/*!
 * \brief Eliminate superfluous node data that is present in the RTT_Format 
 *        mesh file.
 * \param rttMesh Smart pointer to an existing, initialized RTT_Format
 *                class object.
 * \return Mapping between the old and new node numbers.
 */
    vector<int> map_Nodes(const SP<RTT_Format> & rttMesh);


  // generate nodes that are centered on the cell faces, add them to the 
  // cell_pairs vector, calculate their coordinates, and add these to the 
  // vertex vector.
/*!
 * \brief Generate nodes that are centered on the cell faces, adds them to the
 *        cell_pairs vector, calculates their coordinates, and adds the 
 *        coordinates to the vertex vector. Numbering of the face-centererd
 *        nodes begins after all of the cell corner nodes.
 * \param nnodes Total number of nodes (corner nodes only on input and both
 *               corner and face-centered nodes on output).
 * \param rttMesh Smart pointer to an existing, initialized RTT_Format
 *                class object.
 */
    void FC_Nodes(int nnodes, const SP<RTT_Format> & rttMesh);

  // Begin_Doc os_builder-int.tex
  // Begin_Verbatim 

public:
  // constructor
/*!
 * \brief Constructs a CAR_CU_Builder class object.
 * \param interface Smart pointer to an existing, initialized CAR_CU_Interface
 *        class object.
 */
    template<class IT>
    inline explicit CAR_CU_Builder(SP<IT> interface);

/*!
 * \brief Constructs a CAR_CU_Mesh class object.
 * \param rttMesh Smart pointer to an existing, initialized RTT_Format
 *                class object.
 * \return Smart pointer to the new CAR_CU_Mesh class object.
 */
  // build Mesh member functions
    SP<CAR_CU_Mesh> build_Mesh(const SP<RTT_Format> & rttMesh);

  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// inline functions for CAR_CU_Builder
//---------------------------------------------------------------------------//

template<class IT>
inline CAR_CU_Builder::CAR_CU_Builder(SP<IT> interface)
{
    Require (interface);

  // get data arrays from CAR_CU_Interface needed to build CAR_CU_Mesh
    coord_system = interface->get_coordinates();
}

} // end namespace rtt_amr

#endif                          // __amr_CAR_CU_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Builder.hh
//---------------------------------------------------------------------------//
