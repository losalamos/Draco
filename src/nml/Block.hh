//----------------------------------*-C++-*----------------------------------//
// Block.hh
// Geoffrey Furnish
// 15 January 1994
//---------------------------------------------------------------------------//
// @> A class for managing a set of Items.
//---------------------------------------------------------------------------//

#ifndef __nml_Block_hh__
#define __nml_Block_hh__

#include <list>
// using std::list;

#include "Item.hh"
#include "ds++/String.hh"

#include <stdarg.h>

class NML_Block {
    rtt_dsxx::String name;
    std::list< NML_Item * > itmlist;

  public:
    NML_Block( const rtt_dsxx::String& n ) : name(n) {}
    ~NML_Block();

    int add( NML_Item *pitm );
    int modify_defaults( va_list *pap );

    int add_callback( const rtt_dsxx::String& itmnam, NML_Callback& cb );

    const rtt_dsxx::String& Name() const { return name; }
    const std::list<NML_Item *>& Itemlist() const { return itmlist; }

    NML_Item *Get_item( rtt_dsxx::String itmnam );

    void set_defaults();
};

#endif				// __nml_Block_hh__

//---------------------------------------------------------------------------//
//                              end of Block.hh
//---------------------------------------------------------------------------//
