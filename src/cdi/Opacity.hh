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

#ifndef __cdi_Opacityhh__
#define __cdi_Opacity_hh__

#include <iostream>

namespace rtt_cdi
{
 
using std::string;
using std::cout;
using std::endl;

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
	cout << "In Opacity::Opacity() constructor." << endl;
    };
    //defaulted Opacity(const CDI &rhs);
    //defaulted ~Opacity();

    // MANIPULATORS
    
    //defaulted CDI& operator=(const CDI &rhs);

    // ACCESSORS

    virtual double getGray( const double temp, const double density ) = 0;       

    string getDataFilename();

    //  protected:     // DATA
    //    string dataFilename;

  private:
    
    // IMPLEMENTATION
};

// ---------------------------------------- //
//             Gandolf Opacity              //
// ---------------------------------------- //

 class GandolfOpacity : public Opacity
 {
     // NESTED CLASSES AND TYPEDEFS
     
     // DATA

     const string dataFilename;

   public:
     
     // CREATORS
     
     GandolfOpacity( string _data_filename );
     //defaulted GandolfOpacity(const CDI &rhs);
     //defaulted ~GandolfOpacity();
     
     // MANIPULATORS
     
     //defaulted CDI& operator=(const CDI &rhs);
     
     // ACCESSORS
     
     double getGray( const double temp, const double density );
     string getDataFilename() { return dataFilename; };
     
   private:
     
     // IMPLEMENTATION 
 };
 
} // end namespace rtt_cdi

#endif                          // __cdi_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/Opacity.hh
//---------------------------------------------------------------------------//
