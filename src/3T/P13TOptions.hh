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

    bool electronConduction;
    bool ionConduction;

  public:

    // CREATORS

    P13TOptions(bool electronConduction_, bool ionConduction_);
    P13TOptions(const P13TOptions &rhs);
    ~P13TOptions();

    // MANIPULATORS

    P13TOptions& operator=(const P13TOptions &rhs);

    // ACCESSORS

    bool wantElectronConduction() const { return electronConduction; }
    bool wantIonConduction() const { return electronConduction; }
};

#endif                          // __3T_P13TOptions_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/P13TOptions.hh
//---------------------------------------------------------------------------//
