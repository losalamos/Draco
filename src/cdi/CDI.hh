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
    
 enum OpType 
 {
     Gandolf, /*!< Obtain opacity information from an IPCRESS file
		using the Gandolf reader. */
     EOSPAC,  /*!< Obtain opacity information from a ???? file using 
		the EOSPAC reader. */
     Analytic /*!< Use analytic opacities. */
 };

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

    string opacityDataFilename;

    OpType opacityType; // enumerated value.

    SP<Opacity> spOpacity;
    
  public:

    // CREATORS
    
    CDI( OpType _opacity_type, string _opacity_data_filename );
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

    string getOpacityDataFilename() { return opacityDataFilename; };

    vector<int> getMatIDs();

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi

#endif                          // __cdi_CDI_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/CDI.hh
//---------------------------------------------------------------------------//
