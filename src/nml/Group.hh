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
// using std::list;

#include <map>
// using std::map;

#include "Block.hh"
#include "ds++/String.hh"

#include <stdarg.h>

typedef NML_Item * (*NML_Item_handler) ( int type, va_list *pap );

class NML_Group {
    rtt_dsxx::String name;
    std::list<NML_Block *> blklist;
    std::map<rtt_dsxx::String, NML_Block *> blkmap;

    NML_Item_handler user_item_creator;
    NML_Item_handler user_item_modifier;

  public:
    NML_Group( const rtt_dsxx::String& n );
    ~NML_Group();

    int addblock( const char *name, ... );
    void modify_block_defaults( const char *name, ... );

    int delete_block( const rtt_dsxx::String& name );

    int readgroup( const rtt_dsxx::String& file );
    int writegroup( const rtt_dsxx::String& file );

    void set_defaults();

    void write_Tk_info( const rtt_dsxx::String& file );
    void popup_Tk( int argc, char *argv[] ); // {}

    int  add_callback( const rtt_dsxx::String& blknam, const rtt_dsxx::String& itmnam,
		       NML_Callback& cb );
    void add_callback( char *s, void (*f)(void *d), void *_data ); // {}
    void perform_callback( int argc, char *argv[] );

    rtt_dsxx::String Name() const { return name; }
    const std::list<NML_Block *>& Blocklist() const { return blklist; }

    NML_Item *Get_item( rtt_dsxx::String blknam, rtt_dsxx::String itmnam );
    void Set_widget( rtt_dsxx::String blknam, rtt_dsxx::String itmnam, rtt_dsxx::String widget );
};

#endif				// __nml_Group_hh__

//---------------------------------------------------------------------------//
//                              end of Group.hh
//---------------------------------------------------------------------------//
