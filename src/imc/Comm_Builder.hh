//----------------------------------*-C++-*----------------------------------//
// Comm_Builder.hh
// Thomas M. Evans
// Tue Jun  1 17:01:44 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Comm_Builder header file
//---------------------------------------------------------------------------//

#ifndef __imc_Comm_Builder_hh__
#define __imc_Comm_Builder_hh__
 
//===========================================================================//
// class Comm_Builder - 
//
// Purpose : Build a Communicator depending on the data given it, because
//           communicator is templated on Particle Type (PT) so must the
//           builder
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "Communicator.hh"
#include "ds++/SP.hh"
#include <vector>

namespace rtt_imc
{

template<class PT>
class Comm_Builder 
{
  private:
    // usefull typedefs in Comm_Builder
    typedef dsxx::SP<Communicator<PT> > SP_Comm;
    typedef std::vector<int> intvec;
    typedef std::vector<vector<int> > dintvec;

    // given a global cell and a processor, calculate a local cell
    // NOT YET IMPLEMENTED, WE WILL PROBABLY NEED DATA FOR THIS
    int local_cell() { return 0; }

  public:
    // default constructor
    Comm_Builder() {}

    // build_Communicator for full DD
    SP_Comm build_Communicator(const intvec &, const intvec &);

    // build Communicator for general DD/replication
    SP_Comm build_Communicator();
};

} // end namespace rtt_imc

#endif                          // __imc_Comm_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Comm_Builder.hh
//---------------------------------------------------------------------------//
