//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Builder.hh
// B.T. Adams (bta@lanl.gov)
// 18 May 99
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
#include "CAR_CU_Mesh.hh"
#include "ds++/SP.hh"
#include "RTT_Format.hh"
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
using dsxx::SP;

// RTT_Format input file class
using rtt_format::RTT_Format;

class CAR_CU_Builder
{
private:
  // data from Parser needed to build mesh

  // coordinate system string
    string coord_system;
  // vertex (node) coordinate array
    CAR_CU_Mesh::NCVF_d vertex;
  // cell vertexes
    CAR_CU_Mesh::CCVF_i cell_pair;

  // Data needed for localized surface sources
    vector<string> ss_pos;
    vector< vector<int> > defined_surcells;

  // member functions for building CAR_CU_Mesh
  // build Layout helper functions
    SP<Layout> build_Layout(const Coord_sys &, const SP<RTT_Format> & rttMesh);

  // build Coord_sys helper functions
    SP<Coord_sys> build_Coord();

  // assign cells to generation levels in cell_gens and return a vector of 
  // sets of the cells present in each generation (gen_cells);
    vector<set<int> > assign_Generations(CAR_CU_Mesh::CCSF_i & cell_gens, 
					 const SP<RTT_Format> & rttMesh);

  // Establish a hiearchy for the mesh - modifies gen_cells directly and
  // returns a multimap of the new parents paired with their children.
    multimap<int, set<int> > seek_Adoption(vector<set< int> > & gen_cells, 
					   const SP<RTT_Format> & rttMesh);

  // Group the cells into meshes. Return vector of sets containing the cells
  // in each mesh. Can be used to assign child cells or to "connect" the
  // parents, depending upon whether the parents multimap is null or assigned,
  // respectively).
    vector<set< int> > assign_Meshes(const vector<set< int> > & gen_cells, 
				     const multimap<int, set<int> > & parents,
				     const SP<RTT_Format> & rttMesh);

  // Renormalize node numbering to be specific to a mesh
    vector<int> map_Nodes(const SP<RTT_Format> & rttMesh);


  // generate nodes that are centered on the cell faces, add them to the 
  // cell_pairs vector, calculate their coordinates, and add these to the 
  // vertex vector.
    void FC_Nodes(int nnodes, const SP<RTT_Format> & rttMesh);

  // Begin_Doc os_builder-int.tex
  // Begin_Verbatim 

public:
  // constructor
    template<class IT>
    inline explicit CAR_CU_Builder(SP<IT>);

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
//                              end of amr_mesh/CAR_CU_Builder.hh
//---------------------------------------------------------------------------//
