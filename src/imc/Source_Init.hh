//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Source_Init.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 11 09:50:00 2000
 * \brief  Source_Init header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Source_Init_hh__
#define __imc_Source_Init_hh__

#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <string>
#include <vector>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Source_Init
 *
 * \brief Do preliminary source preparations and initializations that require
 * a global mesh.
 * 
 * This class is a prelude to parallel partitioning.  It is designed to
 * perform the following tasks: \arg calculate source particle numbers in
 * each cell so that optimal, parallel-spatial partitioning can be done..
 *
 * Preliminary source data estimated in this class can be used by the
 * Topology_Builder to calculate the parallel-problem topology.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class IT, class MT>
class Source_Init 
{
  public:
    typedef rtt_dsxx::SP<IT>                  SP_Interface;
    typedef rtt_dsxx::SP<MT>                  SP_Mesh;
    typedef std::vector<int>                  sf_int;
    typedef std::vector<std::vector<int> >    vf_int;
    typedef std::vector<std::string>          sf_string;
    typedef std::vector<double>               sf_double;
    typedef std::vector<std::vector<double> > vf_double;
    typedef std::string                       std_string;

  private:
    // Interface object.
    SP_Interface interface;
    
  public:
    //! Constructor.
    Source_Init(SP_Interface);

    // Calculate source numbers as a prelude to parallel, spatial
    // partitioning. 
    void calc_source_numbers() { Insist(0, "Not yet implemented!"); }
};

} // end namespace rtt_imc

#endif                          // __imc_Source_Init_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source_Init.hh
//---------------------------------------------------------------------------//
