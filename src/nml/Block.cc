//----------------------------------*-C++-*----------------------------------//
// Block.cc
// Geoffrey Furnish
// 15 January 1994
//---------------------------------------------------------------------------//
// @> A class for managing a set of Items.
//---------------------------------------------------------------------------//

#include "Block.hh"
#include "Items.hh"

#include "ds++/Assert.hh"

#include <iostream>
using namespace std;

using rtt_dsxx::String;

NML_Block::~NML_Block()
{
    for( list<NML_Item *>::iterator i = itmlist.begin();
         i != itmlist.end(); i++ )
        delete *i;
}

int NML_Block::add( NML_Item *pitm )
{
//     cout << "NML_Block(" << name << ")::add, Adding "
// 	 << pitm->Name() << endl << flush;

    if ( pitm->Name() == "PageBreak" ||
	 pitm->Name() == "BoxStart" ||
	 pitm->Name() == "BoxEnd" ) {

	itmlist.push_back( pitm );
	return 1;
    }

// Disallow duplicate names.

    for( list<NML_Item *>::iterator ili = itmlist.begin();
	 ili != itmlist.end(); ili++ ) {

	NML_Item *pi = *ili;

	if ( pi->Name() == pitm->Name() ) {
	    cerr << "Attempt to add duplicate name to namelist block.\n";
	    cerr << "Request ignored.\n";
	    cerr << "pi->Name() :" << pi->Name() << ": pitm->Name() :"
		 << pitm->Name() << ":.\n";
	    return 0;
	}
    }

    itmlist.push_back( pitm );
    return 1;
}

int NML_Block::modify_defaults( va_list *pap )
{
    for(;;) {
	String itmnam = va_arg( *pap, char * );
	
	if (itmnam == NMI_MODIFY_END)
	  return 1;		// assume nothing went wrong for now.
	
	for( list<NML_Item *>::iterator ili = itmlist.begin();
	     ili != itmlist.end(); ili++ ) {

	    NML_Item *pi = *ili;

	    if ( itmnam == pi->Name() )
		pi->modify_default( pap );
	}
    }
}

int NML_Block::add_callback( const String& itmnam, NML_Callback& cb )
{
    int item_found = 0, retcode = 0;

    for( list<NML_Item *>::iterator ili = itmlist.begin();
	 ili != itmlist.end() && !item_found; ili++ ) {

	NML_Item *pi = *ili;

	if ( pi->Name() == itmnam ) {
	    item_found = 1;
	    retcode = pi->add_callback( cb );
	}
    }

    return retcode;
}

//---------------------------------------------------------------------------//
// Find the item with the specified name, or else throw an exception.
//---------------------------------------------------------------------------//

NML_Item *NML_Block::Get_item( String itmnam )
{
    for( list<NML_Item *>::iterator ili = itmlist.begin();
	 ili != itmlist.end(); ili++ ) {

	NML_Item *pi = *ili;

	if (pi->Name() == itmnam)
	    return pi;
    }

    throw( "Item not found on list." );
}

//---------------------------------------------------------------------------//
// Revert all items back to their default value.
//---------------------------------------------------------------------------//

void NML_Block::set_defaults()
{
    for( list<NML_Item *>::iterator ili = itmlist.begin();
	 ili != itmlist.end(); ili++ )

	(*ili)->set_default();
}

//---------------------------------------------------------------------------//
//                              end of Block.cc
//---------------------------------------------------------------------------//
