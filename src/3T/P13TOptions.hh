//----------------------------------*-C++-*----------------------------------//
// P13TOptions.hh
// Randy M. Roberts
// Thu Mar 12 15:16:01 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_P13TOptions_hh__
#define __3T_P13TOptions_hh__

//===========================================================================//
// class P13TOptions - 
//     This class contains the options with which to run the
//     P13T class.
// 
//===========================================================================//

class P13TOptions {

    // DATA

  private:

    double p1TauMultiplier;
    bool   isCoupledMaterial;
    
  public:

    // CREATORS

    P13TOptions(double p1TauMultiplier_ = 1.0,
		bool isCoupledMaterial_ = true)
	: p1TauMultiplier(p1TauMultiplier_),
	  isCoupledMaterial(isCoupledMaterial_)
    {
	// empty
    }

    // MANIPULATORS

    // ACCESSORS

    double getP1TauMultiplier() const { return p1TauMultiplier; }

    bool getIsCoupledMaterial() const { return isCoupledMaterial; }
};

#endif                          // __3T_P13TOptions_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/P13TOptions.hh
//---------------------------------------------------------------------------//
