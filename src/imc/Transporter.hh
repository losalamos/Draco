//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Transporter.hh
 * \author Thomas M. Evans
 * \date   Thu Apr 13 11:01:53 2000
 * \brief  Transporter class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Transporter_hh
#define rtt_imc_Transporter_hh

#include "ds++/SP.hh"
#include <string>

namespace rtt_mc
{

// Forward declarations.
template<class PT> class Communicator;
template<class PT> class Particle_Containers;

}

namespace rtt_imc
{

// Forward declarations.
class Extrinsic_Surface_Tracker;
template<class MT> class Mat_State;
template<class MT> class Tally;
template<class MT> class Random_Walk;
template<class MT, class FT> class Opacity;
template<class MT, class FT, class PT> class Source;
 
//===========================================================================//
/*!
 * \class Transporter
 *
 * \brief Perform IMC transport -> base class specification.
 *
 * The transporter class requires six fundamental objects: MT, Source,
 * Opacity, Mat_State, Tally, Random_Walk, and Communicator to perform IMC
 * transport.  The fundamental operation of the transporter is to transport
 * all IMC particles from a source in one timestep.
 */
/*!
 * \example imc/test/tstTransporter.cc
 *
 * Transporter unit test.
 */
// revision history:
// -----------------
// 0) original
// 1) 19-MAY-03 : added Random_Walk class to Transporter
// 2) 19-AUG-03 : added surface tracking to Transporter
// 
//===========================================================================//

template<class MT, class FT, class PT>
class Transporter 
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>                                 SP_Mesh;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >                    SP_Opacity;
    typedef rtt_dsxx::SP<Mat_State<MT> >                     SP_Mat_State;
    typedef rtt_dsxx::SP<Source<MT,FT,PT> >                  SP_Source;
    typedef rtt_dsxx::SP<Tally<MT> >                         SP_Tally;
    typedef rtt_dsxx::SP<Random_Walk<MT> >                   SP_Random_Walk;
    typedef rtt_dsxx::SP<Extrinsic_Surface_Tracker>          SP_Tracker;
    typedef rtt_dsxx::SP<rtt_mc::Communicator<PT> >          SP_Communicator;
    typedef typename rtt_mc::Particle_Containers<PT>::Census Census;
    typedef rtt_dsxx::SP<Census>                             SP_Census;
    typedef std::string                                      std_string;

  public:
    // Constructor.
    Transporter() {/*...*/}

    // Destructor.
    virtual ~Transporter() {/*...*/}

    // >>> PUBLIC INTERFACE

    //! Set objects.
    virtual void set(SP_Mesh, SP_Mat_State, SP_Opacity, SP_Source, SP_Tally,
		     SP_Random_Walk, SP_Tracker, SP_Communicator) = 0;

    //! Unset objects.
    virtual void unset() = 0;

    //! Do transport for one time cycle.
    virtual SP_Census transport(double, int, int, int, bool) = 0;

    //! Determine type of transporter.
    virtual std_string type() const = 0;

    //! Query to see if Transporter has all objects assigned.
    virtual bool ready() const = 0;

    //! Get the global, total, number of particles run in the transporter.
    virtual int get_num_run() const = 0;
};

} // end namespace rtt_imc

#endif                          // rtt_imc_Transporter_hh

//---------------------------------------------------------------------------//
//                              end of imc/Transporter.hh
//---------------------------------------------------------------------------//
