//----------------------------------*-C++-*----------------------------------//
// FileSet.cc
// Geoffrey Furnish
// Sat Aug 26 15:47:54 1995
//---------------------------------------------------------------------------//
// @> Class to help manage sets of files.
//---------------------------------------------------------------------------//

#include "util/FileSet.hh"

#include "ds++/Assert.hh"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>

//---------------------------------------------------------------------------//
// Set up the file name components, and scan the directory for matches.
//---------------------------------------------------------------------------//

FileSet::FileSet( const char *stem_, const char *ext_, int nwid_ /*=2*/ )
    : stem(stem_), ext(ext_), nwid(nwid_)
{
    scan();
}

//---------------------------------------------------------------------------//
// Wade through the directory, looking for matches.
//---------------------------------------------------------------------------//

void FileSet::scan()
{
    std::vector<char *> vnames;

    DIR *pdir = opendir( "." );

    struct dirent *d;
    while( (d = readdir(pdir)) != NULL )
    {
	if (select(d))
	    vnames.push_back( d->d_name );
	free( d );
    }

    nfiles = vnames.size();

// Time to sort elements in vnames.;
    std::sort( vnames.begin(), vnames.end() );

// Now load the names DynArray.

    names.low(0);
    names.high(0);

    for( int i=0; i < nfiles; i++ ) {
	names[i] = vnames[i];
	free( vnames[i] );
    }

    find_last_sequence_number();
}

//---------------------------------------------------------------------------//
// Check to see if a given directory entry matches the pattern for this
// fileset. 
//---------------------------------------------------------------------------//

bool FileSet::select( const struct dirent *pdir )
{
    int i, j;
    const char *s = pdir->d_name;

    for( i=0; i < stem.length(); i++ )
	if (s[i] != stem[i]) return false;

    for( j=0; j < nwid; j++ )
	if (!isdigit(s[i+j])) return false;

    if (!strcmp( s+i+j, &ext[0] ))
	return true;

    return false;
}

//---------------------------------------------------------------------------//
// This method finds the sequence number for the last existing file in the
// sequence.  
//---------------------------------------------------------------------------//

int FileSet::find_last_sequence_number()
{
    lstseq = -1;

    if (nfiles) {
	char *s = &names[nfiles-1][0];
	char buf[40];
	s += stem.length();

	int i=0;
	while( isdigit(s[i]) ) {
	    buf[i] = s[i];
	    i++;
	}

	buf[i] = '\0';

	lstseq = atoi( buf );
    }

    return lstseq;
}

//---------------------------------------------------------------------------//
// Return the sequence number for the next file which should be created in
// this sequence, presuming unit sequence increment.  This does /not/ create
// the file, only tell you what it should be.
//---------------------------------------------------------------------------//

int FileSet::next_sequence_number()
{
    return lstseq+1;
}

//---------------------------------------------------------------------------//
// Returns the formatted name for the next sequence number.  Again, this does
// not create it, but rather only says what it should be.
//---------------------------------------------------------------------------//

std::string FileSet::next_sequence_name()
{
    return format_sequence_name( next_sequence_number() );
}

//---------------------------------------------------------------------------//
// Construct a correctly formatted filename in this sequence, given the
// sequence number.
//---------------------------------------------------------------------------//

std::string FileSet::format_sequence_name( int seq )
{
    char buf[20];

    sprintf( buf, "%0*d", nwid, seq );

    std::string r = stem;
    r += buf;
    r += ext;

    return r;
}

//---------------------------------------------------------------------------//
//                              end of FileSet.cc
//---------------------------------------------------------------------------//
