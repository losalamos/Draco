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

static FileSet *active;

int FileSet_select_vector( const struct dirent *d )
{
    return active->select( d );
}

//---------------------------------------------------------------------------//
// Set up the file name components, and scan the directory for matches.
//---------------------------------------------------------------------------//

FileSet::FileSet( String _stem, String _ext, int _nwid /*=2*/ )
    : stem(_stem), ext(_ext), nwid(_nwid)
{
    scan();
}

//---------------------------------------------------------------------------//
// Wade through the directory, looking for matches.
//---------------------------------------------------------------------------//

void FileSet::scan()
{
    int i;
    struct dirent **namelist;

    active = this;

#if 0

#define SELECT_HACK 
#define COMPAR_HACK 

#if defined(_CRAYMPP) || defined(__sgi) || defined(__DECCXX) || \
    (defined(__alpha) && defined(__osf__)) // DEC OSF1 on AlphaServer
#undef SELECT_HACK
#define SELECT_HACK (int (*)( struct dirent * ))
#endif

#if defined(__PGI) || defined(_POWER)
#undef SELECT_HACK
#undef COMPAR_HACK
#define SELECT_HACK (int (*)())
#define COMPAR_HACK (int (*)())
#endif

    nfiles = scandir( ".", &namelist,
		      SELECT_HACK FileSet_select_vector,
		      COMPAR_HACK alphasort );

    if (nfiles < 0) {
	throw( "scandir error!  Maybe invalid directory?" );
    }

    names.low(0);
    names.high(0);

    for( i=0; i < nfiles; i++ ) {
	names[i] = namelist[i]->d_name;
	free( namelist[i] );
    }
    free( namelist );

    find_last_sequence_number();

#endif
}

//---------------------------------------------------------------------------//
// Check to see if a given directory entry matches the pattern for this
// fileset. 
//---------------------------------------------------------------------------//

int FileSet::select( const struct dirent *pdir )
{
    int i, j;
     const char *s = pdir->d_name;

     for( i=0; i < stem.len(); i++ )
	 if (s[i] != stem[i]) return 0;

     for( j=0; j < nwid; j++ )
	 if (!isdigit(s[i+j])) return 0;

     if (!strcmp( s+i+j, &ext[0] ))
	 return 1;

     return 0;
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
	s += stem.len();

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

String FileSet::next_sequence_name()
{
    return format_sequence_name( next_sequence_number() );
}

//---------------------------------------------------------------------------//
// Construct a correctly formatted filename in this sequence, given the
// sequence number.
//---------------------------------------------------------------------------//

String FileSet::format_sequence_name( int seq )
{
    char buf[20];

    sprintf( buf, "%0*d", nwid, seq );

    String r = stem;
    r += buf;
    r += ext;

    return r;
}

//---------------------------------------------------------------------------//
//                              end of FileSet.cc
//---------------------------------------------------------------------------//
