//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/DD_Transporter.hh
 * \author Todd J. Urbatsch and Thomas M. Evans
 * \date   Wed Apr 19 15:37:14 2000
 * \brief  DD_Transporter header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_DD_Transporter_hh__
#define __imc_DD_Transporter_hh__

#include "Transporter.hh"
#include "mc/Topology.hh"
#include "c4/global.hh"
#include <vector>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class DD_Transporter
 *
 * \brief Domain Decomposition topology version of the IMC Transporter.

 * The DD_Transporter derived class is a full DD, or general DD, IMC
 * Transporter.  General DD refers to topologies in which some cells are
 * replicated across processors while others are not.  Full DD refers to
 * topologies in which each cell exists on only one processor.

 */
// revision history:
// -----------------
// 0) original
// 1) 07-24-00 : Incoming buffers are now checked after every single
//               transported particle instead of after every buffer_size
//               particles are run.  This is more robust, but maybe less
//               efficient .
// 
//===========================================================================//

template<class MT, class PT>
class DD_Transporter : public Transporter<MT,PT>
{
  public:
    // Useful typdefs.
    typedef rtt_dsxx::SP<MT>                      SP_Mesh;
    typedef rtt_dsxx::SP<Opacity<MT> >            SP_Opacity;
    typedef rtt_dsxx::SP<Mat_State<MT> >          SP_Mat_State;
    typedef rtt_dsxx::SP<Source<MT,PT> >          SP_Source;
    typedef rtt_dsxx::SP<Tally<MT> >              SP_Tally;
    typedef rtt_dsxx::SP<Communicator<PT> >       SP_Communicator;
    typedef typename Particle_Buffer<PT>::Census  Census;
    typedef typename Particle_Buffer<PT>::Bank    Bank;
    typedef rtt_dsxx::SP<Census>                  SP_Census;
    typedef std::string                           std_string;
    typedef rtt_dsxx::SP<rtt_mc::Topology>        SP_Topology;
    typedef rtt_dsxx::SP<Particle_Buffer<PT> >    SP_Buffer;
    typedef std::vector<C4::C4_Req>               sf_C4_Req;
    typedef std::vector<int>                      sf_int;
    typedef rtt_dsxx::SP<typename PT::Diagnostic> SP_PT_Diagnostic;

  private:
    // Mesh Type object.
    SP_Mesh mesh;

    // Opacity object.
    SP_Opacity opacity;
    
    // Mat_State object.
    SP_Mat_State mat_state;

    // Source object.
    SP_Source source;

    // Tally object.
    SP_Tally tally;

    // Communicator (not used in full replication topology).
    SP_Communicator communicator;

    // Topology (better be full replication)
    SP_Topology topology;

    // Particle Buffer.
    SP_Buffer buffer;

    // Particle counters.
    int finished;
    int num_run;
    int nsrc_run;
    int num_done;
    int num_to_run;

    // Cycle data.
    double delta_t;
    int    cycle;
    int    print_f;

    // C4 Requestors (the numbers refer to message tags).
    sf_C4_Req  rcv500_ndone;
    C4::C4_Req rcv501_fin;

    // Master processors vectors for num_done.
    sf_int recv_num_done;

    // IMC nodes's flag for finished status.
    int recv_finished;

  private:
    // >>> IMPLEMENTATION OF DD TRANSPORT

    // Asynchronous communication functions.
    void trans_src_async(SP_PT_Diagnostic, Bank &, SP_Census);
    void trans_domain_async(SP_PT_Diagnostic, Bank &, SP_Census);
    void update();
    void post_step_arecvs();
    void complete_step_arecvs();

  public:
    // Constructor.
    DD_Transporter(SP_Topology, SP_Buffer);

    // >>> PUBLIC INTERFACE
    
    // Set up objects for this transport step.
    void set(SP_Mesh, SP_Mat_State, SP_Opacity, SP_Source, SP_Tally,
	     SP_Communicator);

    // Unset objects (assign back to null pointers).
    void unset();

    // Run a transport step.
    SP_Census transport(double, int, int, int, bool);

    //! Determine the transport topology.
    std_string type() const { return topology->get_parallel_scheme(); }

    // Determine if we are ready for transport.
    bool ready() const;
};

} // end namespace rtt_imc

#endif                          // __imc_DD_Transporter_hh__

//---------------------------------------------------------------------------//
//                              end of imc/DD_Transporter.hh
//---------------------------------------------------------------------------//
