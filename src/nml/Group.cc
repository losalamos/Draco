//----------------------------------*-C++-*----------------------------------//
// Group.cc
// Geoffrey Furnish
// 15 January 1994
//---------------------------------------------------------------------------//
// @> A class for managing entire namelist groups.
//---------------------------------------------------------------------------//

#include <map>
using std::map;

#include <list>
using std::list;

#include "nml/Group.hh"
#include "nml/Items.hh"
#include "nml/nmstream.hh"

#include <stdarg.h>

#include "ds++/String.hh"
#include "ds++/Assert.hh"

using dsxx::String;

// Maintain a list of pointers to all available namelist groups, so
// that the itcl/tk interface can figure out where to find its data.

map<String, NML_Group *> NML_grplist;

//---------------------------------------------------------------------------//
// Construct a new Group with given name.
//---------------------------------------------------------------------------//

NML_Group::NML_Group( const String& n )
    : name(n),
      user_item_creator( NULL ),
      user_item_modifier( NULL )
{
    if (NML_grplist[ name ] != 0)
	throw "NML_Group: A group by this name already exists.";

// Blunder onward, who knows what can go wrong next.

    NML_grplist[ name ] = this;
}

//---------------------------------------------------------------------------//
// Clean up.
//---------------------------------------------------------------------------//

NML_Group::~NML_Group()
{
    NML_grplist.erase( name );

    for( list<NML_Block *>::iterator i = blklist.begin();
         i != blklist.end(); i++ )
        delete *i;
}

// From the SunOS stdarg manpage:

// The argument list (or its remainder) can be passed to another
// function using a pointer to a variable of type va_list- in which
// case a call to va_arg in the subroutine advances the argument-list
// pointer with respect to the caller as well.

// THE PLAN:

// Okay, the old namelist had horendously complex addblock() and
// modify_block_defaults(), since they processed the entire argument
// list for construction of items of all known types, etc ad nauseum.
// This was ugly, and error prone.

// The new namelist will make this much cleaner by merely handing off
// the argument list to the appropriate constructor, and let it handle
// parsing of the args needed to init itself.  This shoudl make for
// /much/ simpler code all the way around.

//---------------------------------------------------------------------------//
// addblock

// Add a new namelist block to the group.  The variable length
// argument list specifies the type and initialization list for each
// item to be added to the block.  The available item types are listed
// in "nml/Items.hh", but there is also a facility for user-definable
// types (not finished yet).  The form of the initialization data is
// type dependent; see "nml/Items.cc".
//---------------------------------------------------------------------------//

int NML_Group::addblock( const char *name, ... )
{
// See if there is an existing block by this name.

    NML_Block *pb = NULL;

// If not, start a new one.

    pb = blkmap[ name ];
    if (!pb) {
	pb = new NML_Block( name );
	blklist.push_back( pb );
	blkmap[ name ] = pb;
    }

// Now process the args, creating items in turn.

    va_list ap;
    va_start(ap,name);

    int done=0;
    NML_Item *pitm;

    while( !done ) {
	int itmtype = va_arg(ap,int);
	pitm = NULL;

	switch( itmtype ) {

	case NMI_INT:
	    pitm = new nmi_int( &ap );
	    break;

	case NMI_FLOAT:
	    pitm = new nmi_float( &ap );
	    break;

	case NMI_DOUBLE:
	    pitm = new nmi_double( &ap );
	    break;

	case NMI_INT_SET:
	    pitm = new nmi_int_set( &ap );
	    break;

	case NMI_INT_ARRAY:
	    pitm = new nmi_int_array( &ap );
	    break;

	case NMI_INT_DYNARRAY:
	    pitm = new nmi_DynArray_int( &ap );
	    break;

	case NMI_BOOL:
	    pitm = new nmi_bool( &ap );
	    break;

	case NMI_FLOG:
	    pitm = new nmi_flog( &ap );
	    break;

	case NMI_STRING:
	    pitm = new nmi_String( &ap );
	    break;

// Add more standard types here.

	case NMI_PAGE_BREAK:
	    pitm = new nmi_PageBreak( &ap );
	    break;

	case NMI_BOX_START:
	    pitm = new nmi_BoxStart( &ap );
	    break;

	case NMI_BOX_END:
	    pitm = new nmi_BoxEnd( &ap );
	    break;

	case NMI_END:
	    done = 1;
	    break;

// We allow the user to create his own Items if he chooses, by
// specifying a handler to parse the unrecognized parameters.
// Alternatively, signal an error, probably severe.

	default:
	    if (user_item_creator)
		pitm = (*user_item_creator) ( itmtype, &ap );
	    else {
		cerr << "Unknown item type :" << itmtype << ":, discarding\n";
		cerr << " token.  This error is probably unrecoverable,\n";
		cerr << "but being the hopelessly naive and optimistic\n";
		cerr << "chaps that we are, we'll just keep blasting\n";
		cerr << "along...\n\n";
	    }
	    break;
	}

	if (pitm) pb->add( pitm );
    }

    va_end(ap);

    return 1;			// Should return something to indicate
				// if we had trouble or not.
}

//---------------------------------------------------------------------------//
// modify_block_defaults

//  Change the default value for specified items in the specified
//  block.  
//---------------------------------------------------------------------------//

void NML_Group::modify_block_defaults( const char *name, ... )
{
// See if there is an existing block by this name.

    map< String, NML_Block *>::iterator i = blkmap.find( name );
    NML_Block *pb = (*i).second;
    if ( !pb ) throw "No such block.";

// Now start processing of the va_list, but leave it to the block to
// do the actual parsing.

    va_list ap;
    va_start(ap,name);

    pb->modify_defaults( &ap );

    va_end(ap);
}

//---------------------------------------------------------------------------//
// delete_block

// This method is responsible for discarding a namelist block which is
// no longer needed with this group.  Since the List<T> class does not
// provide a simple way to delete items on a List, this is going to be
// a fair amount of work.
//---------------------------------------------------------------------------//

int NML_Group::delete_block( const String& name )
{
    list< NML_Block * > nblist;

    int block_deleted = 0;

    for( list<NML_Block *>::iterator bli = blklist.begin();
	 bli != blklist.end(); bli++ ) {
	NML_Block *pb = *bli;

	if ( pb->Name() != name )
	    nblist.push_back( pb );
	else {
	    delete pb;
	    block_deleted++;
	}
    }

// Okay, we've constructed a new list containing all the old blocks
// except for the one to be discarded (and we've destroyed that one),
// so now we should be able to just assign the new list to the
// permanent one.  I just hope this does't recursively delete the
// members of the old list, thus wrecking the new one.  Have to check
// that sometime.

    blklist = nblist;

    blkmap.erase( name );

    return block_deleted == 1;
}

//---------------------------------------------------------------------------//
// Read a namelist group from a file.
//---------------------------------------------------------------------------//

int NML_Group::readgroup(  const String& file )
{
    nmistream s( file );

    s >> *this;

    return 1;
}

//---------------------------------------------------------------------------//
// Write this group out to a file.
//---------------------------------------------------------------------------//

int NML_Group::writegroup(  const String& file )
{
    nmostream of( file, "w" );

    of << *this;

    return 1;
}

//---------------------------------------------------------------------------//
// Reset all items in all blocks back to their default values.
//---------------------------------------------------------------------------//

void NML_Group::set_defaults()
{
    for( list<NML_Block *>::iterator bli = blklist.begin();
	 bli != blklist.end(); bli++ )
	(*bli)->set_defaults();
}

//---------------------------------------------------------------------------//
// add_callback

// Add a callback to a specific item of a specific block in this group.
//---------------------------------------------------------------------------//

int  NML_Group::add_callback( const String& blknam, const String& itmnam,
			      NML_Callback& cb )
{
    return blkmap[ blknam ]->add_callback( itmnam, cb );
}

//---------------------------------------------------------------------------//
// Get_item

// Find and return a namelist item in a specified block.  Note:
// getting a bit more brave--throw an exception if we have trouble
// locating the data.  Note 2:  This function would be about 3 lines
// if we were using an exception enabled Map class to store the
// info...  Have to look into that.
//---------------------------------------------------------------------------//

NML_Item *NML_Group::Get_item( String blknam, String itmnam )
{
    return blkmap[ blknam ]->Get_item( itmnam );
}

void NML_Group::Set_widget( String blknam, String itmnam,
			    String widget )
{
    NML_Item *pi = Get_item( blknam, itmnam );

    pi->Set_widget( widget );
}

//---------------------------------------------------------------------------//
//                              end of Group.cc
//---------------------------------------------------------------------------//
