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
// 1)  9-16-98 : added local cell -> global cell update for census particles
//               in step_IMC_dd
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
#include "c4/global.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <string>
#include <numeric>
#include <iostream>

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
public:
  // some usefull typedefs
    typedef typename Particle_Buffer<PT>::Comm_Buffer Comm_Buffer;
    typedef typename Particle_Buffer<PT>::Census      Census;
    typedef typename Particle_Buffer<PT>::Bank        Bank;
    typedef typename PT::Diagnostic                   PT_Diagnostic;

  // dd communication controller
    class DD_Comm
    {
    private:
      // accumulated status reports
	vector<int> status;
      // local status reports
	int *node_status;
      // local status indicator
	int local_status;
      // number of sends in this timestep
	int sends;
	
      // illegal services
	DD_Comm(const DD_Comm &);
	const DD_Comm& operator=(const DD_Comm &);

    public:
      // constructor
	inline DD_Comm();
	inline ~DD_Comm();
	
      // async communication accounting
	inline void update_send(const vector<int> &);
	inline void update_send(const int);
	inline void update_recv(const int);

      // communication terminating indicator
	inline bool terminate();

      // accessor functions
	int get_sends() const { return sends; }
    };
	
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
    SP<Comm_Buffer> new_census_buffer;
    SP<Census> new_census_bank;

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

  // DD transport functions
    void dd_loop(Bank &, SP<PT_Diagnostic>, int &, DD_Comm &); 
    void dd_Particle_transport(SP<PT>, SP<PT_Diagnostic>, Bank &, int &,
			       DD_Comm &);

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

//---------------------------------------------------------------------------//
// DD_Comm nested class inline functions
//---------------------------------------------------------------------------//
// constructor

template<class MT, class BT, class IT, class PT> inline 
IMC_Manager<MT,BT,IT,PT>::DD_Comm::DD_Comm()
    : status(C4::nodes()), local_status(0), sends(0)
{
  // make local node status object that we will use for updates between
  // synchronizations of the timestep
    node_status = new int[C4::nodes()];

  // initialize node status arrays to zero
    for (int i = 0; i < C4::nodes(); i++)
    {
	status[i]      = 0;
	node_status[i] = 0;
    }
    
    Ensure (!std::accumulate(status.begin(), status.end(), 0));
}

//---------------------------------------------------------------------------//
// destructor

template<class MT, class BT, class IT, class PT>
inline IMC_Manager<MT,BT,IT,PT>::DD_Comm::~DD_Comm()
{
  // reclaim memory
    delete [] node_status;
}

//---------------------------------------------------------------------------//
// update the node_status for single particle Communications

template<class MT, class BT, class IT, class PT>
inline void IMC_Manager<MT,BT,IT,PT>::DD_Comm::update_send(const int node)
{
    Require (node < C4::nodes());
    
  // update the node_status array with the node we sent a particle buffer to
    if (node >= 0)
    {
	node_status[node]++;
	local_status = 1;
	sends++;
    }
}

//---------------------------------------------------------------------------//
// update the node_status after a flush of Communications	

template<class MT, class BT, class IT, class PT> inline 
void IMC_Manager<MT,BT,IT,PT>::DD_Comm::update_send(const vector<int> &nodes) 
{
  // update the node_status array with the nodes we sent particles to
    int global_node;
    for (int i = 0; i < nodes.size(); i++)
    {
	global_node = nodes[i];
	Check (global_node < C4::nodes());
	if (global_node >= 0)
	{
	    node_status[global_node]++;
	    local_status = 1;
	    sends++;
	}
    }
}

//---------------------------------------------------------------------------//
// update the node_status if we have received something

template<class MT, class BT, class IT, class PT> inline 
void IMC_Manager<MT,BT,IT,PT>::DD_Comm::update_recv(const int num_arrived)
{
    Require (num_arrived < C4::nodes());

  // subtract the number of times this node has received data
    if (num_arrived > 0)
    {
	node_status[C4::node()] -= num_arrived;
	local_status = 1;
    }
}

//---------------------------------------------------------------------------//
// determine if we need to terminate or keep going

template<class MT, class BT, class IT, class PT> inline 
bool IMC_Manager<MT,BT,IT,PT>::DD_Comm::terminate()
{
  // add up the local status data on each node
    C4::gsum(node_status, C4::nodes());
    C4::gsum(local_status);

  // by default we will assume we are done and then find out otherwise
    bool all_done = true;
    for (int i = 0; i < nodes(); i++)
    {
	status[i] += node_status[i];
	if (status[i] != 0)
	    all_done = false;
	node_status[i] = 0;
    }

  // determine if we need to stop
    if (all_done && local_status == 0)
	return true;
    else 
    {
      // we are not done so return false and reset the local_status
	local_status = 0;
    }

  // we aren't done so return false
    return false;
}
	
CSPACE

#endif                          // __imc_IMC_Manager_hh__

//---------------------------------------------------------------------------//
//                              end of imc/IMC_Manager.hh
//---------------------------------------------------------------------------//
