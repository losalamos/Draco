//----------------------------------*-C++-*----------------------------------//
// Group.hh
// Geoffrey Furnish
// 15 January 1994
//---------------------------------------------------------------------------//
// @> A class for managing entire namelist groups.
//---------------------------------------------------------------------------//

#ifndef __nml_Group_hh__
#define __nml_Group_hh__

#include <list>
using std::list;

#include <map>
using std::map;

#include "nml/Block.hh"
#include "ds++/String.hh"

#include <stdarg.h>

typedef NML_Item * (*NML_Item_handler) ( int type, va_list *pap );

class NML_Group {
    String name;
    list<NML_Block *> blklist;
    map<String, NML_Block *> blkmap;

    NML_Item_handler user_item_creator;
    NML_Item_handler user_item_modifier;

  public:
    NML_Group( const String& n );
    ~NML_Group();

    int addblock( const char *name, ... );
    void modify_block_defaults( const char *name, ... );

    int delete_block( const String& name );

    int readgroup( const String& file );
    int writegroup( const String& file );

    void set_defaults();

    void write_Tk_info( const String& file );
    void popup_Tk( int argc, char *argv[] ); // {}

    int  add_callback( const String& blknam, const String& itmnam,
		       NML_Callback& cb );
    void add_callback( char *s, void (*f)(void *d), void *_data ); // {}
    void perform_callback( int argc, char *argv[] );

    String Name() const { return name; }
    const list<NML_Block *>& Blocklist() const { return blklist; }

    NML_Item *Get_item( String blknam, String itmnam );
    void Set_widget( String blknam, String itmnam, String widget );
};

#endif				// __nml_Group_hh__

//---------------------------------------------------------------------------//
//                              end of nml/Group.hh
//---------------------------------------------------------------------------//
