//----------------------------------*-C++-*----------------------------------//
// ADFile.hh
// Geoffrey Furnish
// Fri Dec 23 07:53:46 1994
//---------------------------------------------------------------------------//
// @> An "amorphous data file" for GTS.
//---------------------------------------------------------------------------//

#ifndef __util_ADFile_hh__
#define __util_ADFile_hh__

#include "ds++/String.hh"
#include "ds++/Mat.hh"

#include "util/Indx_File.hh"
#include "util/ADKey.hh"

//===========================================================================//
// class ADFile - Amorphous Data File

// This class enables you to read and write an indexed file of amorphous data.
// In otherwords, you use this when you expect to put data of all different
// types in the file, assigning each a key.  The index key, "ADKey", is
// essentially a flat, fairly large, character array.  This allows you to
// create keys which contain descriptive information about each data record
// written to the file.  This info can be useful in helping to identify the
// data and determine how best to read it back later.  An "Amorphous Data
// eXamination" (ADX) utility is provided seperately.
//===========================================================================//

class ADFile : public Indx_File<ADKey> {

// Nobody has any business trying to copy this thing that I can see.

    ADFile( const ADFile& );
    ADFile& operator=( const ADFile& );

    String fname;
    int nrecs, mode;

  public:
    ADFile( String _fname, int _mode, int _nrecs );

    Mat1<ADKey> Get_keys();
};

#endif                          // __util_ADFile_hh__

//---------------------------------------------------------------------------//
//                              end of util/ADFile.hh
//---------------------------------------------------------------------------//
