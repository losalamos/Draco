//----------------------------------*-C++-*----------------------------------//
// pcg_DB.hh
// Dave Nystrom
// Fri May  9 09:45:08 1997
//---------------------------------------------------------------------------//
// @> pcg descriptor block.
//---------------------------------------------------------------------------//

#ifndef __linalg_pcg_DB_hh__
#define __linalg_pcg_DB_hh__

#include "util/persist.hh"

#include "nml/Group.hh"
#include "nml/Items.hh"

#include "ds++/String.hh"

//===========================================================================//
// class pcg_DB - pcg descriptor block

// Holds namelist parameters for setting up and controlling pcg.
//===========================================================================//

class pcg_DB : public persist {
  public:
// Some way to name this thing.

    String name;

// Pick up the namelist variables.

#include "linalg/.nml_pcg.hh"

// Pick up the auxiliary variables.

// Methods.

  public:
    pcg_DB( String _name ) : name(_name) {}

    String Name() const { return name; }

    void setup_namelist( NML_Group& g );

    void save   ( ofstream& ofs );
    void restore( ifstream& ifs );
};

#endif                          // __linalg_pcg_DB_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/pcg_DB.hh
//---------------------------------------------------------------------------//
