//----------------------------------*-C++-*----------------------------------//
// ADKey.hh
// Geoffrey Furnish
// Fri Dec 23 09:13:24 1994
//---------------------------------------------------------------------------//
// @> Amorphous data description key
//---------------------------------------------------------------------------//

#ifndef __util_ADKey_hh__
#define __util_ADKey_hh__

//===========================================================================//
// class ADKey - Amorphous data descriptor key

// This is the key used in conjunction with the Amorphous Data File, ADFile.
// Essentially the key is a large character array which can be filled in via
// the usual methods (sprintf, etc).  The keys are then used by the ADFile
// class for tagging data.
//===========================================================================//

class ADKey {
  public:
    enum { BUFSIZE = 100 };

    char s[ BUFSIZE ];

// Utility mehtods.

    int operator==( const ADKey& key ) const;
};

#endif                          // __util_ADKey_hh__

//---------------------------------------------------------------------------//
//                              end of util/ADKey.hh
//---------------------------------------------------------------------------//
