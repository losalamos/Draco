//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/NoPrecon.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 09:48:40 2000
 * \brief  No (identity) preconditioner
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_NoPrecon_hh__
#define __ConjGrad_NoPrecon_hh__

namespace rtt_ConjGrad
{
 
//===========================================================================//
/*!
 * \class NoPrecon
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Field>
class NoPrecon 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    //Default: NoPrecon();
    //Default: NoPrecon(const NoPrecon &rhs);
    //Default: ~NoPrecon();

    // MANIPULATORS
    
    //Default: NoPrecon& operator=(const NoPrecon &rhs);

    // ACCESSORS

    Field operator()(const Field b) const
    {
	return b;
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGrad

#endif                          // __ConjGrad_NoPrecon_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/NoPrecon.hh
//---------------------------------------------------------------------------//
