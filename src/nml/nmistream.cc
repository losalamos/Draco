//---------------------------------*-C++-*---------------------------------//
// nmistream.cc
// Geoffrey Furnish
// 3 June 1992
//-------------------------------------------------------------------------//
// @> Implementation of abstracted input stream type for namelist
// @> library. Loosely based on version for old nml, but rewritten
// @> extensively to use libds++ features.
//-------------------------------------------------------------------------//

#define ATT_io
#include <iostream.h>
#include <fstream.h>

#include <ctype.h>

#include "nml/nmstream.hh"

#define MAX_ITEM_LEN 40
#define MAX_LINE_LEN 200

//-------------------------------------------------------------------------//
// This function reads in a namelist block from an input stream.
// The precise operation of the input stream on any given system is
// handled by the other member functions below.
//-------------------------------------------------------------------------//

nmistream& nmistream::operator>>( NML_Block& b )
{
    char c;

// The size is no longer a fixed requirement, just a heuristic to get
// us started.

    String vname( MAX_ITEM_LEN ), value( MAX_ITEM_LEN );

    char buf[ MAX_LINE_LEN ];

// We begin right after the block name identifier.

    while(get(c)) {
	if (isspace(c)) continue;

	if (c == '#') {		// "Comment to end of line" introducer.
	    getline( buf, 200 );
	    continue;
        }

	vname = c;		// Collect the value name.
	while(get(c)) {
	    if (c != '=' && !isspace(c))
	      vname += c;
	    else
	      break;
	}

	// Check for the terminator string(s).  Currently a terminator
	// _is_ required.  It might be nice to drop this restriction
	// eventually, but I will have to think about that quite a bit
	// more before I can implement something good.  So, for the
	// time being...

	if ( vname == "$end" || vname == "$" )
	  break;

	while(isspace(c)) get(c); // Skip any whitespace after
				    // a token identifier.
	
	if (c != '=') {		// Attempt disaster recover.  Who
				// knows if this will work?
	    cout << "Unexpected character :" << c
		 << ":.  Skipping to =.\n";
	    while( c != '=' ) get(c);
	}

// Find the item by this name.

	list<NML_Item *> itmlist( b.Itemlist() );	
	int found = 0;

	for( list<NML_Item *>::iterator ili = itmlist.begin();
	     ili != itmlist.end(); ili++ ) {
	    NML_Item *pi = *ili;

	    if ( pi->Name() == vname ) {
		found=1;
		pi->Read_Value( *this ); // Now read it in.
		break;
	    }
	}

	if (!found) {
	    cerr << "** ERROR ** No item " << vname << " found in block "
		 << b.Name() << '\n';

	// At this point we have found a bogus variable name and an =.  So we
	// really ought to try to parse past the supplied value in a best
	// effort to avoid botching the whole input deck.  If we just
	// continue from here, the supplied value will be another bogus input
	// char, and it will just keep going south from there...

	    cerr << "Attempting to skip forward to next valid token.\n";

	    do {
		get(c);
	    } while( isspace(c) );

	    if (c == '\"') {
	    // Absorb up to the next " char.
		do {
		    get(c);
		} while( c != '\"' );
	    }
	    else if (c == '{') {
	    // Absorb up to the closing curly brace.
		do {
		    get(c);
		} while( c != '}' );
	    }
	    else {
	    // Doesn't appear to be any form of quoted entity, so just suck
	    // up characters until the next whitespace.

		do {
		    get(c);
		} while( !isspace(c) );
	    }
	}
    }
    return *this;
}

//-------------------------------------------------------------------------//
// Here we read in a namelist group.
//-------------------------------------------------------------------------//

nmistream& nmistream::operator>>( NML_Group& g )
{
    char c;
    char buf[   MAX_LINE_LEN ];

    String bname( MAX_ITEM_LEN );

    while( get(c) ) {		// Find a namelist block.
	if (isspace(c))
	  continue;

	if (c == '#') {		// Comment to end of line.
	    getline( buf, 200 );
	    continue;
        }

	if (c == '$') {		// Start of a block name?
	    bname = "";
	    // Get the word...
	    while(get(c)) {
		if (isspace(c)) break;
		bname += c;
	    }

	    list<NML_Block *> blklist = g.Blocklist();

	    int found = 0;

	    for( list<NML_Block *>::iterator bli = blklist.begin();
		 !found && bli != blklist.end(); bli++ ) {
		NML_Block *pb = *bli;
		if (bname == pb->Name() ) {
		    found = 1;
		    (*this) >> *pb;
		}
	    }

	    if (!found) {
		cout << "No block by the name of :" << bname <<
		  ": has been defined.\n";
		cout << "Probably gonna get totally fouled"
		  " up from here on...\n";
		cout << "Will attempt recovery by skipping to end of next\n";
		cout << "   token which begins with a '$'.\n";
		cout << "   Cross your fingers.\n";

		while( c != '$' ) get(c);
		while( !isspace(c) ) get(c);
		// Done all we can do.	Hope for the best...;
	    }
	    continue;
	}
	cout << "Character :" << c << ": not expected. \n";
    }
    return *this;
}

//-------------------------------------------------------------------------//
// Now comes the dicey stuff, for different platforms.
//-------------------------------------------------------------------------//

#ifdef ATT_io

nmistream::nmistream( const char *name, char *mode )
{
    inf.open( name, ios::in );

    if (!inf)
	cerr << "Failed to open input stream " << name << ".\n";
}

nmistream::~nmistream(void)
{
    inf.close();
}

void nmistream::close(void)
{
    inf.close();
}

int nmistream::get( char &c )
{
    if (inf.get(c))
      return 1;
    else 
      return 0;
}

int nmistream::getline( char *buf, int len )
{
    if (inf.getline( buf, len ))
      return 1;
    else
      return 0;
}

void nmistream::putback( char c )
{
    inf.putback( c );
}

#endif				// ATT_io

//-------------------------------------------------------------------------//
//                         end of nmistream.cc
//-------------------------------------------------------------------------//
