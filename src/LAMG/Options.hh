//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/Options.hh
 * \author Randy M. Roberts
 * \date   Tue Jan 25 17:01:32 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMG_Options_hh__
#define __LAMG_Options_hh__

namespace rtt_LAMG
{
 
//===========================================================================//
/*!
 * \class Options
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Options 
{

    // NESTED CLASSES AND TYPEDEFS

  public:
    
    enum Outlevel { SILENT=0, LEVEL1=1, LEVEL2=2, LEVEL3=3};

    // DATA

  private:
    
    int itsmax_m;
    double tol_m;
    Outlevel levout_m;
    
  public:

    // CREATORS
    
    Options(int itsmax_in=1000, double tol_in=1.0e-4,
	    Outlevel levout_in=SILENT)
	: itsmax_m(itsmax_in), tol_m(tol_in), levout_m(levout_in)
    {
	// empty
    }
    
    //DEFAULT: Options(const Options &rhs);
    //DEFAULT: ~Options();

    // MANIPULATORS

    Options &itsmax(int itsmax_in)
    {
	itsmax_m = itsmax_in;
	return *this;
    }
    
    Options &tol(double tol_in)
    {
	tol_m = tol_in;
	return *this;
    }
    
    Options &levout(Outlevel levout_in)
    {
	levout_m = levout_in;
	return *this;
    }
    
    //DEFAULT: Options& operator=(const Options &rhs);

    // ACCESSORS

    int itsmax() const { return itsmax_m; }
    double tol() const { return tol_m; }
    Outlevel levout() const { return levout_m; }
    
  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_LAMG

#endif                          // __LAMG_Options_hh__

//---------------------------------------------------------------------------//
//                              end of LAMG/Options.hh
//---------------------------------------------------------------------------//
