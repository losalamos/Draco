//----------------------------------*-C++-*----------------------------------//
// Item.hh
// Geoffrey Furnish
// 15 January 1994
//---------------------------------------------------------------------------//
// @> Declares the NML_Item base class.
//---------------------------------------------------------------------------//

#ifndef __nml_Item_hh__
#define __nml_Item_hh__

#include <list>
// using std::list;

#include "ds++/String.hh"

#include "nml/Callback.hh"

#include <stdarg.h>

class nmistream;

class NML_Item {

  protected:
    dsxx::String name;
    dsxx::String value;		// (The current value).
    dsxx::String default_value;
    std::list< NML_Callback > cblist;

    dsxx::String widget;
    int active;

  public:
    int defered_callback;
    dsxx::String callback_binding;

  public:
    NML_Item( const dsxx::String& _name,
	      const dsxx::String& _defval )
	: name(_name), default_value(_defval),
	  defered_callback(0)
    {}
    virtual ~NML_Item() {}

// Since most derived types will be determining their name and
// default_value through va_arg processing, it isn't really reasonable
// to expect them to know these during the initialilzation phase.  So
// we provide a default ctor, and permit "lazy init's".

    NML_Item() : widget("Default") {}

    void set_default();
    void set_new_default( const dsxx::String& ndef, int update_now =1 );
    virtual void modify_default( va_list *pap ) =0;

    void set_value( const dsxx::String& newval ) { value = newval; }
    virtual void update_value() =0;

    virtual dsxx::String File_Rep();
    virtual void Read_Value( nmistream& nis );

// Some methods for dealing with callbacks.

    void clear_callbacks();
    int  add_callback( NML_Callback& cb );
    void delete_callback( NML_Callback& cb );
    void do_callbacks();

    const dsxx::String& Name() const { return name; }
    const dsxx::String& Value() const { return value; }
    const dsxx::String& Default() const { return default_value; }

// Some methods for helping with Tk interaction.

    virtual dsxx::String Type() const =0;
    virtual dsxx::String Choices() const  { return ""; }
    virtual dsxx::String Widget() const { return widget; }

    void Set_widget( dsxx::String w ) { widget = w; }

    // Might need a way to pass in configuration data too.  Could add
    // a config datum to this class, and provide a way to export it to
    // the Tcl side.  Or, could try to figure out how to use the
    // option database, etc.  Let's wait for [incr Tk] before doing
    // anything rash on this front.

// Initial support for item supression.

    int Active() const { return active; }
    void Activate() { active = 1; }
    void Deactivate() { active = 0; }
};

#endif				// __nml_Item_hh__

//---------------------------------------------------------------------------//
//                              end of nml/Item.hh
//---------------------------------------------------------------------------//
