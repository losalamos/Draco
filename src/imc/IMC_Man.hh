//----------------------------------*-C++-*----------------------------------//
// IMC_Man.hh
// Thomas M. Evans
// Wed Jun  3 10:36:11 1998
//---------------------------------------------------------------------------//
// @> IMC_Man class header file.
//---------------------------------------------------------------------------//

#ifndef __imc_IMC_Man_hh__
#define __imc_IMC_Man_hh__

//===========================================================================//
// class IMC_Man - 
//
// Purpose : Manager for running IMCTEST package.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imc/Names.hh"
#include "imc/Opacity_Builder.hh"
#include "imc/Source_Init.hh"
#include "imc/Parallel_Builder.hh"
#include "imc/Particle_Buffer.hh"
#include "imc/Particle.hh"
#include "imc/Source.hh"
#include "imc/Tally.hh"
#include "imc/Mat_State.hh"
#include "imc/Opacity.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"

IMCSPACE

// draco components
using RNG::Rnd_Control;

// template manager on: MT=mesh type; BT=mesh builder type; IT=interface
// type; PT=particle type
template<class MT, class BT, class IT, class PT = Particle<MT> >
class IMC_Man 
{
private:
  // objects used by all processors
    SP<MT> mesh;
    SP<Opacity<MT> > opacity;
    SP<Mat_State<MT> > mat_state;
    SP<Rnd_Control> rnd_con;
    SP<Parallel_Builder<MT> > parallel_builder;
    SP<Particle_Buffer<PT> > buffer;

  // objects used only by the host
    SP<Source_Init<MT> > source_init;
    SP<typename Mat_State<MT>::Shell> global_state;

  // problem variables
    double delta_t;
    int cycle;
    int max_cycle;

  // verbosity switch
    bool verbose;

  // some service functions for conserving memory
    template<class T> inline void kill(SP<T> &spref);
    
public:
  // default constructor
    IMC_Man(bool = false);

  // run everything over the requisite number of cycles
    void execute_IMC(char *);

  // initialize the problem on the host processors
    void host_init(char *);

  // initialize the problem on the IMC processors
    void IMC_init();

  // run the problem for one time-cycle
    void run_IMC();

  // do collect stuff at the end of the timestep
    void regroup();

  // do the output and end
    void output();
};

//---------------------------------------------------------------------------//
// inline functions
//---------------------------------------------------------------------------//
// release a smart pointer, also, if this is the last SP to a particular
// object, the whole thing will be deleted

template<class MT, class BT, class IT, class PT>
template<class T>
inline void IMC_Man<MT,BT,IT,PT>::kill(SP<T> &spref)
{
  // assigning this SP to a null SP
    spref = SP<T>();
    Ensure (!spref);
}

CSPACE

#endif                          // __imc_IMC_Man_hh__

//---------------------------------------------------------------------------//
//                              end of imc/IMC_Man.hh
//---------------------------------------------------------------------------//
