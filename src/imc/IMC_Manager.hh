//----------------------------------*-C++-*----------------------------------//
// IMC_Manager.hh
// Thomas M. Evans
// Wed Jun  3 10:36:11 1998
//---------------------------------------------------------------------------//
// @> IMC_Manager class header file.
//---------------------------------------------------------------------------//

#ifndef __imc_IMC_Manager_hh__
#define __imc_IMC_Manager_hh__

//===========================================================================//
// class IMC_Manager - 
//
// Purpose : Manager for running IMCTEST package as a standalone module.
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
#include "imc/Global_Tally.hh"
#include "imc/Communicator.hh"
#include "imc/Global.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <string>

IMCSPACE

// draco components
using RNG::Rnd_Control;
using dsxx::SP;

// stl components
using std::string;

// template manager on: MT=mesh type; BT=mesh builder type; IT=interface
// type; PT=particle type
template<class MT, class BT, class IT, class PT = Particle<MT> >
class IMC_Manager 
{
private:
  // objects used by all processors
    SP<MT> mesh;
    SP<Opacity<MT> > opacity;
    SP<Mat_State<MT> > mat_state;
    SP<Rnd_Control> rnd_con;
    SP<Parallel_Builder<MT> > parallel_builder;
    SP<Particle_Buffer<PT> > buffer;
    SP<Source<MT> > source;
    SP<Tally<MT> > tally;
    SP<Communicator<PT> > communicator;

  // objects used only by the host
    SP<Source_Init<MT> > source_init;
    SP<Global_Tally<MT> > global_state;

  // census objects
    SP<typename Particle_Buffer<PT>::Comm_Buffer> new_census_buffer;
    SP<typename Particle_Buffer<PT>::Census> new_census_bank;

  // problem variables
    double delta_t;
    int cycle;
    int max_cycle;
    int print_f;
    int dump_f;
    int rnstream;
    string parallel_scheme;

  // verbosity switch
    bool verbose;

  // convert the parallel_scheme to an int for easy sending
    inline int get_scheme(string) const;
    inline string get_scheme(int) const;

  // some service functions for conserving memory
    template<class T> inline void kill(SP<T> &spref);
    
public:
  // default constructor
    IMC_Manager(bool = false);

  // run everything over the requisite number of cycles
    void execute_IMC(char *);

  // initialize the problem on the host processors
    void host_init(char *);

  // initialize the problem on the IMC processors
    void IMC_init();

  // run the problem for one time-cycle
    void step_IMC_rep();
    void step_IMC_dd();

  // do collect stuff at the end of the timestep
    void regroup();

  // do the output and end
    void output();

  // print diagnostics
    void cycle_dump() const;
    void verbose_dump() const;	
};

//---------------------------------------------------------------------------//
// inline functions
//---------------------------------------------------------------------------//
// release a smart pointer, also, if this is the last SP to a particular
// object, the whole thing will be deleted

template<class MT, class BT, class IT, class PT>
template<class T>
inline void IMC_Manager<MT,BT,IT,PT>::kill(SP<T> &spref)
{
  // assigning this SP to a null SP
    spref = SP<T>();
    Ensure (!spref);
}

//---------------------------------------------------------------------------//
// convert the parallel_scheme string into an integer

template<class MT, class BT, class IT, class PT>
inline int IMC_Manager<MT,BT,IT,PT>::get_scheme(string ps) const
{
  // replication = 1 : DD = 2 : DD/replication = 3
    int value;
    if (ps == "replication")
	value = 1;
    else if (ps == "DD")
	value = 2;
    else if (ps == "DD/replication")
	value = 3;
    else 
	Check (0);
    return value;
}

//---------------------------------------------------------------------------//
// convert an int back into an integer

template<class MT, class BT, class IT, class PT>
inline string IMC_Manager<MT,BT,IT,PT>::get_scheme(int ps) const
{
    string value;
    if (ps == 1)
	value = "replication";
    else if (ps == 2)
	value = "DD";
    else if (ps == 3)
	value == "DD/replication";
    else 
	Check (0);
    return value;
}

CSPACE

#endif                          // __imc_IMC_Manager_hh__

//---------------------------------------------------------------------------//
//                              end of imc/IMC_Manager.hh
//---------------------------------------------------------------------------//
