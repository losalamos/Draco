//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Builder.hh
// B.T. Adams (bta@lanl.gov)
// 18 May 99
//---------------------------------------------------------------------------//
// @> CAR_CU_Builder class header file
//---------------------------------------------------------------------------//

#ifndef __mc_CAR_CU_Builder_hh__
#define __mc_CAR_CU_Builder_hh__

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

namespace rtt_mc 
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

  // return the face name for given face number.
    string get_face_name(int face) const;

  // Begin_Doc os_builder-int.tex
  // Begin_Verbatim 

public:
  // constructor
    template<class IT>
    inline explicit CAR_CU_Builder(SP<IT>);

  // build Mesh member functions
    SP<CAR_CU_Mesh> build_Mesh(const SP<RTT_Format> & rttMesh);
    
  // return the number of grouped surface source cells
    int get_ss_size() { return ss_pos.size(); }
  // return the number of grouped surface source cells in a given set
    int get_ss_size(int surface) 
    { return defined_surcells[surface - 1].size(); }
  // return the position (lox, hix, etc.) of a set of grouped surface source 
  // cells
     string get_ss_pos(int surface) { return ss_pos[surface - 1]; }
  // return the positions (lox, hix, etc.) of the all of the grouped surface 
  // source cells
    vector<string> get_ss_pos() { return ss_pos; }
  // return the defined surface source cells in a given set
    vector<int> get_defined_surcells(int surface) 
    {
        vector<int> source_set(defined_surcells[surface - 1].size());
	for (int cell = 0; cell < defined_surcells[surface - 1].size(); cell++)
	    source_set[cell] = defined_surcells[surface - 1][cell];

        return source_set; 
    }
  // return all of the defined surface cells read from the RTT Format file
    vector<vector<int> > get_defined_surcells() { return defined_surcells; }

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
    ss_pos = interface->get_ss_pos();
    defined_surcells.resize(ss_pos.size());
}

} // end namespace rtt_mc

#endif                          // __mc_CAR_CU_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of mc/CAR_CU_Builder.hh
//---------------------------------------------------------------------------//
