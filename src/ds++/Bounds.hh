//----------------------------------*-C++-*----------------------------------//
// Bounds.hh
// Geoffrey Furnish
// Mon Jan 27 12:11:50 1997
//---------------------------------------------------------------------------//
// @> Class to represent index ranges for DS++ containers.
//---------------------------------------------------------------------------//

#ifndef __ds_Bounds_hh__
#define __ds_Bounds_hh__

#include "config.hh"

NAMESPACE_DS_BEG

//===========================================================================//
// class Bounds - Index range specifier class

// This class holds the inclusive bounds specifiers used to construct DS++
// containers when something other than zero based indexing is desired.  For
// example, Bounds(-1,3) represents an index range running from -1 through 3, 
// or a total of 5 valid indexes.  
//===========================================================================//

class Bounds {
    int a, b;
  public:
    explicit Bounds( int a_, int b_ ) : a(a_), b(b_) {}
    int min() const { return a; }
    int max() const { return b; }
    int len() const { return b - a + 1; }
};

NAMESPACE_DS_END

#endif                          // __ds_Bounds_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/Bounds.hh
//---------------------------------------------------------------------------//
