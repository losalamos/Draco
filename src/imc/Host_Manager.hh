//----------------------------------*-C++-*----------------------------------//
// Host_Manager.hh
// Thomas M. Evans
// Sun Aug  2 12:43:05 1998
//---------------------------------------------------------------------------//
// @> Host_Manager class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Host_Manager_hh__
#define __imc_Host_Manager_hh__

//===========================================================================//
// class Host_Manager - 
//
// Purpose : run IMC through an appropriate host code
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "Names.hh"
#include "Opacity_Builder.hh"
#include "Parallel_Source_Init.hh"
#include "Particle_Buffer.hh"
#include "Particle.hh"
#include "Source.hh"
#include "Tally.hh"
#include "Mat_State.hh"
#include "Opacity.hh"
#include "Communicator.hh"
#include "Global.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"

IMCSPACE

// draco components
using RNG::Rnd_Control;
using dsxx::SP;

template<class MT, class BT, class IT, class PT = Particle<MT> >
class Host_Manager
{
private:
  // objects needed for IMC_Transport
    SP<MT> mesh;
    SP<Opacity<MT> > opacity;
    SP<Mat_State<MT> > mat_state;
    SP<Rnd_Control> rnd_con;
    SP<Particle_Buffer<PT> > buffer;
    SP<Source<MT> > source;
    SP<Tally<MT> > tally;
    SP<Communicator<PT> > communicator;
   
  // initialization objects
    SP<Parallel_Source_Init<MT> > source_init;

  // census objects
    SP<typename Particle_Buffer<PT>::Census> new_census_bank;

  // problem variables
    double delta_t;
    int cycle;
    int dump_f;

  // some service functions for conserving memory
    template<class T> inline void kill(SP<T> &spref);

public:
  // default constructor
    Host_Manager(int);

  // run one cycle
    void execute_IMC(typename IT::Arguments &);

  // initialize this problem cycle
    void initialize(const typename IT::Arguments &);

  // run the problem for one time_cycle
    void step_IMC();

  // do any collections or cleanup necessary
    void regroup(typename IT::Arguments &);

  // print diagnostics
    void cycle_dump() const;
};

//---------------------------------------------------------------------------//
// inline functions
//---------------------------------------------------------------------------//
// release a smart pointer, also, if this is the last SP to a particular
// object, the whole thing will be deleted

template<class MT, class BT, class IT, class PT>
template<class T>
inline void Host_Manager<MT,BT,IT,PT>::kill(SP<T> &spref)
{
  // assigning this SP to a null SP
    spref = SP<T>();
    Ensure (!spref);
}

CSPACE

#endif                          // __imc_Host_Manager_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Host_Manager.hh
//---------------------------------------------------------------------------//
