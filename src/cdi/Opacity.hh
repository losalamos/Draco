//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/Opacity.hh
 * \author Kelly Thompson
 * \date   Thu Jun 23 13:55:06 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_Opacity_hh__
#define __cdi_Opacity_hh__

#include <iostream>
#include <vector>

namespace rtt_cdi
{
 
using std::string;
using std::cout;
using std::endl;
using std::vector;

//===========================================================================//
/*!
 * \class Opacity
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

// ---------------------------------------- //
//       Opacity Base Class (abstract)      //
// ---------------------------------------- //

class Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

  public:

    // CREATORS
    
    Opacity() 
    {
	// empty
	cout << "In Opacity::Opacity() constructor." << endl;
    }
    //defaulted Opacity(const CDI &rhs);
    virtual ~Opacity()
    {
	// empty
	cout << "Destroying Opacity Object." << endl << endl;
    }

    // MANIPULATORS
    
    //defaulted CDI& operator=(const CDI &rhs);

    // ACCESSORS

    virtual string getDataFilename() = 0;

    virtual vector<int> getMatIDs() = 0;

    virtual double getGray( const double temp, const double density ) = 0;       

    //  protected:     

    // DATA

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi

#endif                          // __cdi_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/Opacity.hh
//---------------------------------------------------------------------------//
