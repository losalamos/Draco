//----------------------------------*-C++-*----------------------------------//
// FileSet.hh
// Geoffrey Furnish
// Sat Aug 26 15:47:54 1995
//---------------------------------------------------------------------------//
// @> Class to help manage sets of files.
//---------------------------------------------------------------------------//

#ifndef __util_FileSet_hh__
#define __util_FileSet_hh__

#include <dirent.h>
#include <ctype.h>

#include "ds++/DynArray.hh"

#include <string>

//===========================================================================//
// class FileSet - Help manage sets of files

// This class encapsulates the scandir function in libc, providing an
// interface to that function which makes it easier to manage sets of files.
// In particular, we want to have ordered sequences of files which share a
// common prefix stem, and a common extension, and which are distinguished by
// a numeric field in the name.  For example, 
//        phi00.dat phi01.dat phi02.dat
// Functions are provided to find the existing files which match the pattern,
// find the first or last in the sequence, find the next file which should be
// in the sequence, and format a name with a given sequence number.
//===========================================================================//

class FileSet {

    FileSet( const FileSet& );
    FileSet& operator=( const FileSet& );

    std::string stem;
    std::string ext;
    int nwid;			// width of numeric field.

    dsxx::DynArray<std::string> names;
    int nfiles;

    int lstseq;

  public:
    FileSet( const char *stem_, const char *ext_, int nwid_ =2 );

    void scan();
    bool select( const struct dirent *pdir );

    int find_last_sequence_number();
    int      next_sequence_number();
    std::string   next_sequence_name();
    std::string format_sequence_name( int seq );

    int matching_files() const { return nfiles; }
    std::string name( int i ) const { return names[i]; }
};

#endif                          // __util_FileSet_hh__

//---------------------------------------------------------------------------//
//                              end of util/FileSet.hh
//---------------------------------------------------------------------------//
