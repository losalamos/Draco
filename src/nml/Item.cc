//----------------------------------*-C++-*----------------------------------//
// Item.cc
// Geoffrey Furnish
// 15 January 1994
//---------------------------------------------------------------------------//
// @> Implements the NML_Item base class.
//---------------------------------------------------------------------------//

#include "Item.hh"
using std::list;

#include "nmstream.hh"

#include <ctype.h>

using rtt_dsxx::String;

void NML_Item::set_default()
{
    value = default_value;
    update_value();
}

void NML_Item::set_new_default( const String& ndef, int update_now )
{
    default_value = ndef;

    if (update_now)
	set_default();
}

//---------------------------------------------------------------------------//
// This method is responsible for producing the representation that an item
// has on disk.  Used for writing.  If the simplistic representation provided
// by this base class is inadequate, a derived class may overload.
//---------------------------------------------------------------------------//

String NML_Item::File_Rep()
{
    String r = name + "=" + value;

    return r;
}

//---------------------------------------------------------------------------//
// Read the new value from the input stream.  This provides a basic facility
// which is probably appropriate for most items.  However, individual item
// types may wish to use a special/peculiar specification format.  Such cases
// can be handled by overloading this method in the derived class.
//---------------------------------------------------------------------------//

void NML_Item::Read_Value( nmistream& nis )
{
    String v;
    char c;
    
    // On entry, we're on the '='.  Move past it.

    while( nis.get(c) ) {
	if (isspace(c)) continue;
	nis.putback(c);
	break;
    }

    // Now we're at the beginning of the value string.

    while( nis.get(c) ) {
	if (isspace(c)) {
	    nis.putback(c);
	    break;
	} else
	    v += c;
    }

    set_value( v );
    update_value();
    do_callbacks();
}

//---------------------------------------------------------------------------//
// Some methods for dealing with callbacks.
//---------------------------------------------------------------------------//

void NML_Item::clear_callbacks()
{
    cblist.erase( cblist.begin(), cblist.end() );
}

int NML_Item::add_callback( NML_Callback& cb )
{
    cblist.push_back( cb );

    return 1;			// What could go wrong ???
}

void NML_Item::delete_callback( NML_Callback& cb )
{
    for( list< NML_Callback >::iterator cbli = cblist.begin();
	 cbli != cblist.end(); cbli++ )
	if ( cb == *cbli )
	    cblist.erase( cbli );
}

void NML_Item::do_callbacks()
{
    for( list< NML_Callback >::iterator cbli = cblist.begin();
	 cbli != cblist.end(); cbli++ )
	(*cbli).invoke( value );
}

//---------------------------------------------------------------------------//
//                              end of Item.cc
//---------------------------------------------------------------------------//
