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

#ifndef __imc_Transporter_hh__
#define __imc_Transporter_hh__

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
template<class MT> class Mat_State;
template<class MT> class Tally;
template<class MT, class FT> class Opacity;
template<class MT, class FT, class PT> class Source;
 
//===========================================================================//
/*!
 * \class Transporter
 *
 * \brief Perform IMC transport -> base class specification.
 *
 * The transporter class requires five fundamental objects: MT, Source,
 * Opacity, Mat_State, Tally, and Communicator to perform IMC transport.  The
 * fundamental operation of the transporter is to transport all IMC particles
 * from a source in one timestep.
 */
// revision history:
// -----------------
// 0) original
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
		     SP_Communicator) = 0;

    //! Unset objects.
    virtual void unset() = 0;

    //! Do transport for one time cycle.
    virtual SP_Census transport(double, int, int, int, bool) = 0;

    //! Determine type of transporter.
    virtual std_string type() const = 0;

    //! Query to see if Transporter has all objects assigned.
    virtual bool ready() const = 0;
};

} // end namespace rtt_imc

#endif                          // __imc_Transporter_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Transporter.hh
//---------------------------------------------------------------------------//
