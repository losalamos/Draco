//---------------------------------*-C++-*---------------------------------//
// nmostream.cc
// Geoffrey Furnish
// 3 June 1992
//-------------------------------------------------------------------------//
// @> Implementation of abstracted output stream type for namelist
// @> library.  Loosely based on version for old nml, but rewritten
// @> extensively to use libds++ features.
//-------------------------------------------------------------------------//

#define ATT_io
#include <iostream.h>
#include <fstream.h>

#include "nml/nmstream.hh"
#include <list>
using std::list;

using dsxx::String;

//-------------------------------------------------------------------------//
// Dump out a block.
//-------------------------------------------------------------------------//

nmostream& nmostream::operator<<( const NML_Block& b )
{
    (*this) << "\n$" << b.Name() << "\n\t";

    int pos = 2 + b.Name().len();
    int j = (pos-1) % 8;
    if (j) pos += 8-j;
// This is all wrong.  Must recompute the value.


// output all the items.;

    list<NML_Item *> itmlist = b.Itemlist();
    for( list<NML_Item *>::iterator ili = itmlist.begin();
	 ili != itmlist.end(); ili++ ) {
	NML_Item *pi = *ili;

	String out = pi->File_Rep();
	int len = out.len();

	if (out[len-2] == '\n' && out[len-1] == '\t') {
	    (*this) << out;
	    pos = 9;
	    continue;
	}
	if (pos + out.len() > 77 ) {
	    (*this) << "\n\t";
	    pos = 9;
	}
	(*this) << out << " ";
	pos += out.len() + 1;
    }

//    if (pos > 65) (*this) << '\n';

    (*this) << "\n$end\n\n";

    return *this;
}

//-------------------------------------------------------------------------//
// Write an entire namelist group out to disk file.
//-------------------------------------------------------------------------//

nmostream& nmostream::operator<<( const NML_Group& g )
{
    (*this) << "#\n";
    (*this) << "# " << g.Name() << '\n';
    (*this) << "#\n";

    list<NML_Block *> blklst = g.Blocklist();
    for( list<NML_Block *>::iterator bli = blklst.begin();
	 bli != blklst.end(); bli++ ) {
	NML_Block *pb = *bli;
	(*this) << *pb;
    }

    (*this) << "#\n";
    (*this) << "# End of group " << g.Name() << '\n';
    (*this) << "#\n";
    
    return *this;
}

//-------------------------------------------------------------------------//
// Now comes the dicey stuff, the output for each crazy architecture.
//-------------------------------------------------------------------------//

#ifdef ATT_io
     
nmostream::nmostream( const char *name, char *mode )
{
    of.open( name, ios::out );

    if (!of) 
      cout << "Failed to open " << name << " for output.\n";
}

nmostream::~nmostream()
{
    close();
}

void nmostream::close()
{
    of.close();
}

nmostream& nmostream::operator<<( const char& c )
{
    of << c;
    return *this;
}

nmostream& nmostream::operator<<( const char *s )
{
    of << s;
    return *this;
}

nmostream& nmostream::operator<<( const String& s )
{
    const char *p = s;
    return (*this) << p;
}

#endif				// ATT_io

//-------------------------------------------------------------------------//
//                         end of nmostream.cc
//-------------------------------------------------------------------------//
