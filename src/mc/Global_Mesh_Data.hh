//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Global_Mesh_Data.hh
 * \author Thomas M. Evans and Todd Urbatsch
 * \date   Thu Dec  4 17:36:07 2003
 * \brief  Global_Mesh_Data class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Global_Mesh_Data_hh
#define rtt_mc_Global_Mesh_Data_hh

#include <vector>
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "Topology.hh"

namespace rtt_mc
{

// Forward declarations.
class OS_Mesh;
class RZWedge_Mesh;
class Sphyramid_Mesh;

//===========================================================================//
/*!
 * \class Global_Mesh_Data
 * \brief Class that stores global data on meshes regardless of
 * decomposition.
 *
 * This class provides a mechanism for storing global mesh data regardless of
 * the parallel decomposition.  For example, it is sometimes necessary to
 * know the global-extents of a mesh.  This class can be handed an
 * arbitrarily decomposed mesh and the global extents can be calculated.
 *
 * The global mesh information that is currently available is:
 * - global spatial extents
 * .
 *
 * This class may do parallel communication in the constructor.
 *
 * \sa Global_Mesh_Data.t.hh, Global_Mesh_Data.cc for detailed descriptions.
 */
/*! 
 * \example mc/test/tstGlobal_Mesh_Data.cc 
 * 
 * description of example
 */
// revision history:
// -----------------
// 0) (Thu Dec  4 17:36:07 2003) Thomas M. Evans: original
// 
//===========================================================================//

template<class MT>
class Global_Mesh_Data 
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<Topology> SP_Topology;
    typedef std::vector<double>    sf_double;

  private:
    // >>> DATA 

    // Global mesh extents.
    sf_double spatial_extents;

    // Topology.
    SP_Topology topology;

  private:
    // >>> IMPLEMENTATION

    // Calculate global mesh data.
    void calc_global_mesh_data(const MT &);

  public:
    // Constructor.
    Global_Mesh_Data(SP_Topology, const MT &);

    // >>> ACCESSORS

    //! Get the global spatial_extents.
    const sf_double& get_spatial_extents() const { return spatial_extents; }
};

//---------------------------------------------------------------------------//
// MEMBER SPECIALIZATIONS ON MESH TYPE
//---------------------------------------------------------------------------//

// OS_Mesh specialization.
template<>
void Global_Mesh_Data<OS_Mesh>::calc_global_mesh_data(const OS_Mesh &mesh);

// RZWedge_Mesh specialization.
template<>
void Global_Mesh_Data<RZWedge_Mesh>::calc_global_mesh_data(
    const RZWedge_Mesh &mesh);

// Sphyramid_Mesh specialization.
template<>
void Global_Mesh_Data<Sphyramid_Mesh>::calc_global_mesh_data(
    const Sphyramid_Mesh &mesh);

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
// TEMPLATE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "Global_Mesh_Data.i.hh"

#endif // rtt_mc_Global_Mesh_Data_hh

//---------------------------------------------------------------------------//
//              end of mc/Global_Mesh_Data.hh
//---------------------------------------------------------------------------//
