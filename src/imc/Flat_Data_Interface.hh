//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Flat_Data_Interface.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 14 16:33:35 2001
 * \brief  Interface definition for Flat_Mat_State_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Flat_Data_Interface_hh__
#define __imc_Flat_Data_Interface_hh__

#include "Flat_Data_Container.hh"
#include "ds++/SP.hh"

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Flat_Data_Interface
 *
 * \brief Additional interface functions required by the
 * Flat_Mat_State_Builder.
 *
 * This interface class contains additional definitions (beyond those
 * required by rtt_imc::Interface) required by the Flat_Mat_State_Builder.
 * To make a class an interface for the Flat_Mat_State_Builder it must
 * inherit from both rtt_imc::Interface and rtt_imc::Flat_Data_Interface.
 *
 * All fields are returned as std::vector<double> types.  Fields are
 * cell-centered. 
 *
 * There are no name collisions between this interface definition,
 * CDI_Data_Interface, and rtt_imc::Interface.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Flat_Data_Interface 
{ 
  private:
    // Useful typedefs.
    typedef rtt_dsxx::SP<Flat_Data_Container> SP_Flat_Data_Container;
    
  public:
    //! Constructor.
    Flat_Data_Interface() { /* no data to construct */ }

    //! Virtual constructor to make life happy down the inhertiance chain.
    virtual ~Flat_Data_Interface() { }

    // >>> PUBLIC INTERFACE REQUIRED BY FLAT_MAT_STATE_BUILDER

    //! Get a Smart-Pointer to the Flat_Data_Container.
    virtual SP_Flat_Data_Container get_flat_data_container() const = 0;
};

} // end namespace rtt_imc

#endif                          // __imc_Flat_Data_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Flat_Data_Interface.hh
//---------------------------------------------------------------------------//
