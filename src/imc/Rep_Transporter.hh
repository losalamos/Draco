//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Rep_Transporter.hh
 * \author Thomas M. Evans
 * \date   Thu Apr 13 11:41:37 2000
 * \brief  Rep_Transporter header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Rep_Transporter_hh__
#define __imc_Rep_Transporter_hh__

#include "Transporter.hh"
#include "mc/Topology.hh"

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Rep_Transporter
 *
 * \brief Full replication topology version of the IMC Transporter.
 *
 * The Rep_Transporter derived class is a full topology-replicated IMC
 * transporter.  It expects the full mesh and data to be copied on each
 * processor.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, class PT> 
class Rep_Transporter : public Transporter<MT,PT>
{
  public:
    // Useful typdefs.
    typedef rtt_dsxx::SP<MT>                            SP_Mesh;
    typedef rtt_dsxx::SP<Opacity<MT> >                  SP_Opacity;
    typedef rtt_dsxx::SP<Mat_State<MT> >                SP_Mat_State;
    typedef rtt_dsxx::SP<Source<MT,PT> >                SP_Source;
    typedef rtt_dsxx::SP<Tally<MT> >                    SP_Tally;
    typedef rtt_dsxx::SP<rtt_mc::Communicator<PT> >     SP_Communicator;
    typedef typename rtt_mc::Particle_Stack<PT>::Census Census;
    typedef rtt_dsxx::SP<Census>                        SP_Census;
    typedef std::string                                 std_string;
    typedef rtt_dsxx::SP<rtt_mc::Topology>              SP_Topology;

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

  public:
    // Constructor.
    Rep_Transporter(SP_Topology);

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

#endif                          // __imc_Rep_Transporter_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Rep_Transporter.hh
//---------------------------------------------------------------------------//
