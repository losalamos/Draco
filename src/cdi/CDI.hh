//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:06 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_CDI_hh__
#define __cdi_CDI_hh__

#include <iostream>
#include <vector>

#include "ds++/SP.hh"

#include "Opacity.hh"

namespace rtt_cdi
{
 
using std::string;
using std::vector;
using rtt_dsxx::SP;
 
// DATA
    
//===========================================================================//
/*!
 * \class CDI
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class CDI 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    SP<Opacity> spOpacity;
    
  public:

    // CREATORS
    
    // constructor acts as creator?
    CDI( SP<Opacity> _spOpacity );
    
    //defaulted CDI(const CDI &rhs);
    ~CDI()
    {
	//delete pOpacity;
	cout << "Destroying CDI Object." << endl << endl;
    };

    // MANIPULATORS
    
    //defaulted CDI& operator=(const CDI &rhs);

    // ACCESSORS

    double getGrayOpacity( const double temp, const double density );

    string getOpacityDataFilename() { return spOpacity->getDataFilename(); };

    vector<int> getMatIDs();

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi

#endif                          // __cdi_CDI_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/CDI.hh
//---------------------------------------------------------------------------//
