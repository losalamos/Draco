//----------------------------------*-C++-*----------------------------------//
// Items.cc
// Geoffrey Furnish
// 25 January 1994
//---------------------------------------------------------------------------//
// @> Implements the basic Item types.
//---------------------------------------------------------------------------//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "Items.hh"
#include "nmstream.hh"

#include <map>
using std::map;

using rtt_dsxx::DynArray;
using rtt_dsxx::String;

//===========================================================================//
// class nmi_int

// A basic namelist class for handling integers.
//===========================================================================//

nmi_int::nmi_int( va_list *pap )
{
    name = va_arg( *pap, char * );
    int dv = va_arg( *pap, int );
    pv = va_arg( *pap, int * );

    char buf[40];
    sprintf( buf, "%d", dv );
    default_value = buf;

    widget = "Entry";

    set_default();
}

void nmi_int::modify_default( va_list *pap )
{
    char buf[40];
    sprintf( buf, "%d", va_arg( *pap, int ) );

    set_new_default( buf );    
}

void nmi_int::update_value()
{
    *pv = atoi( value );	// String class needs cast to 
				// (const char *)
}

//===========================================================================//
// class nmi_float

// A basic namelist class for handling floats.
//===========================================================================//

nmi_float::nmi_float( va_list *pap )
{
    name = va_arg( *pap, char * );
    char buf[40];
    sprintf( buf, "%f", va_arg( *pap, double ) ); // float->double
						  // promo in va_list.
    pv = va_arg( *pap, float * );
    default_value = buf;

    widget = "Entry";

    set_default();
}

void nmi_float::modify_default( va_list *pap )
{
    char buf[40];
    sprintf( buf, "%f", va_arg( *pap, double ) );

    set_new_default( buf );
}

void nmi_float::update_value()
{
    *pv = atof( value );
}

//===========================================================================//
// class nmi_double

// A basic namelist class for handling doubles.
//===========================================================================//

nmi_double::nmi_double( va_list *pap )
{
    name = va_arg( *pap, char * );
    char buf[40];
    sprintf( buf, "%lf", va_arg( *pap, double ) );
    pv = va_arg( *pap, double * );
    default_value = buf;

    widget = "Entry";

    set_default();
}

void nmi_double::modify_default( va_list *pap )
{
    char buf[40];
    sprintf( buf, "%lf", va_arg( *pap, double ) );

    set_new_default( buf );
}

void nmi_double::update_value()
{
    *pv = atof( value );	// Loss of precision?
}

//===========================================================================//
// class nmi_int_set

// A basic namelist class for handling sets of named integers.
//
// Currently this class works by having the user specify the /numeric/
// value of the default to the ctor, as well as to the modify default
// method.  This seems to be most convenient for the user, who is
// presumably using enums or some such.  But it requires inverting the
// map, which is expensive.  Since this class is the prototype for
// nmi_bool and nmi_flog, I may want to simplify this by changing the
// type of the default specification to String.  Have to think about
// it some more.
//===========================================================================//

// This ctor is used by derived classes who know what the option set
// should look like upon entry.  Probably they won't know the desired
// default value or storage location, however, so taking care of those
// things  is left up to them.

nmi_int_set::nmi_int_set( int undef, int nvals, ... )
{
    va_list ap;
    va_start(ap, nvals);

    for( int i=0; i < nvals; i++ ) {
	String vname = va_arg( ap, char * );
	int    value = va_arg( ap, int );
	m[vname] = value;
    }

// pv not set!  Must be done by derived class.
}

// This is the publicly visible ctor.  Argument list must be of the
// form:
//     char *            name of item
//     char *            default string
//     int               "undefined" value: returned on invalid input
//     int               number of members
//     int *             storage location
//     <char *,int> ...  pairs of name/value

nmi_int_set::nmi_int_set( va_list *pap )
{
    name          = va_arg( *pap, char * );
    default_value = va_arg( *pap, char * );
    int undef     = va_arg( *pap, int );
    int nvals     = va_arg( *pap, int );
    pv            = va_arg( *pap, int * );

    map<String,int> mm;

    for( int i=0; i < nvals; i++ ) {
	String vn = va_arg( *pap, char * );
	int val   = va_arg( *pap, int );
	mm[vn] = val;
    }

    m = mm;

    widget = "NML_IntSet";
    widget = "IntSet";

    set_default();
}

void nmi_int_set::modify_default( va_list *pap )
{
    set_new_default( va_arg( *pap, char * ) );
}

void nmi_int_set::update_value()
{
    *pv = m[value];
}

String nmi_int_set::Choices() const
{
    String choices;

    for( map<String,int>::const_iterator p = m.begin();
	 p != m.end(); p++ )
	choices += (*p).first + ' ';

    return choices;
}

//===========================================================================//
// class nmi_int_array

// A basic namelist class for handling arrays of integers.
//===========================================================================//

nmi_int_array::nmi_int_array( va_list *pap )
{
    name = va_arg( *pap, char * );
    default_value = va_arg( *pap, char * );
    pv = va_arg( *pap, int * );
    max_vals = va_arg( *pap, int );
    cur_vals = va_arg( *pap, int * );

    widget = "Entry";

    set_default();
}

void nmi_int_array::modify_default( va_list *pap )
{
    default_value = va_arg( *pap, char * );
}

#define MAX_ITEM_LEN 40

void nmi_int_array::update_value()
{
    int i=0, j, k=0;
    char str[ MAX_ITEM_LEN ];

    while( value[i] && k < max_vals ) {
    // Although we know any values we read in have no leading
    // whitespace, this is not necessarily true of the default
    // value.  So have to check. for leading whitespace.

	while( isspace(value[i]) ) i++;
	if (!value[i]) break;

    // Currently there is not supposed to be any space between
    // numbers, since a space terminates the token parser, but
    // hopefully I'll get around to relaxing that eventually.
    // So, just go ahead and make this understand formats of
    // the form #,#,# or # # #.

	j=0;
	while( value[i] && value[i] != ',' && !isspace(value[i]) ) 
	  str[j++] = value[i++]; // Copy the nest # into str.

	str[j] = '\0';

	pv[k++] = atoi( str );
	if (value[i] == ',')
	  i++;			// Skip over the , if it was the seperator.
    }
    *cur_vals = k;
}

//===========================================================================//
// class nmi_bool

// A basic namelist class for boolean values.
//===========================================================================//

nmi_bool::nmi_bool( va_list *pap )
    : nmi_int_set( 0, 16,
		   "off", 0,
		   "n", 0,
		   "no", 0,
		   "f", 0,
		   ".f", 0,
		   "false", 0,
		   ".false",0,
		   "0", 0,
		   "on", 1,
		   "y", 1,
		   "yes", 1,
		   "t", 1,
		   ".t", 1,
		   "true", 1,
		   ".true", 1,
		   "1", 1 )
{
    name = va_arg( *pap, char * );
    default_value = va_arg( *pap, char * );
    pv = va_arg( *pap, int * );

    widget = "Bool";

    set_default();
}

void nmi_bool::modify_default( va_list *pap )
{
    set_new_default( va_arg( *pap, char * ) );
}

void nmi_bool::update_value()
{
    *pv = m[value];
}

//===========================================================================//
// class nmi_flog

// A basic namelist class for working with Fortran logicals.
//===========================================================================//

nmi_flog::nmi_flog( va_list *pap )
    : nmi_int_set( 0, 14,
		   "off", FLOG_FALSE,
		   "n", FLOG_FALSE,
		   "no", FLOG_FALSE,
		   "f", FLOG_FALSE,
		   ".f", FLOG_FALSE,
		   "false", FLOG_FALSE,
		   ".false",FLOG_FALSE,
		   "on", FLOG_TRUE,
		   "y", FLOG_TRUE,
		   "yes", FLOG_TRUE,
		   "t", FLOG_TRUE,
		   ".t", FLOG_TRUE,
		   "true", FLOG_TRUE,
		   ".true", FLOG_TRUE )
{
    name = va_arg( *pap, char * );

    default_value = va_arg( *pap, char * );
    pv = va_arg( *pap, int * );

    widget = "Bool";

    set_default();
}

void nmi_flog::modify_default( va_list *pap )
{
    set_new_default( va_arg( *pap, char * ) );
}

void nmi_flog::update_value()
{
    *pv = m[value];
}

//===========================================================================//
// class nmi_String

// A basic namelist class for handling String values.
//===========================================================================//

nmi_String::nmi_String( va_list *pap )
{
    name = va_arg( *pap, char * );

    default_value = va_arg( *pap, char * );
    ps = va_arg( *pap, String * );

    widget = "Entry";

    set_default();
}

void nmi_String::modify_default( va_list *pap )
{
    set_new_default( va_arg( *pap, char * ) );
}

void nmi_String::update_value()
{
    *ps = value;
}

//---------------------------------------------------------------------------//
// Since pretty much anything could be in the string, we output using the
// quote delimited method.  See Read_Value() below.
//---------------------------------------------------------------------------//

String nmi_String::File_Rep()
{
    String r = name + "=\"" + value + "\"";

    return r;
}

//---------------------------------------------------------------------------//
// Since a string value can be pathological, such as by being empty, or by
// having space characters, etc, we need special handling.  We treat two
// cases:  if the first character is a " sign, we slurp to the next one,
// otherwise, we slurp to the next whitespace.
//---------------------------------------------------------------------------//

void nmi_String::Read_Value( nmistream& nis )
{
    String v;
    char c;

// We start on the =.  Skip over any whitespace.

    while( nis.get(c) ) {
	if (isspace(c)) continue;
	break;
    }

    if (c == '\"') {		// Quote delimited.
	while( nis.get(c) ) {

	// Should handle embedded quotes and other special characters, but
	// just haven't gotten to that yet.  For now, terminate on next ".

	    if (c == '\"')
		break;
	    else
		v += c;
	}
    } else {			// Whitespace delimited.
	v = c;
	while( nis.get(c) ) {
	    if (isspace(c)) {
		nis.putback(c);
		break;
	    }
	    v += c;
	}
    }

    set_value( v );
    update_value();
    do_callbacks();
}

//===========================================================================//
// class nmi_DynArray_int - For i/o of DynArray<int>

// A class for robust namelist i/o of integer arrays.  We use a DS++
// DynArray<int> class to ensure that there is no dependence on the number of
// values specified, etc.
//===========================================================================//

nmi_DynArray_int::nmi_DynArray_int( va_list *pap )
{
    name = va_arg( *pap, char * );

    default_value = va_arg( *pap, char * );
    pdai = va_arg( *pap, DynArray<int> * );

    widget = "Entry";

    set_default();
}

void nmi_DynArray_int::modify_default( va_list *pap )
{
    set_new_default( va_arg( *pap, char * ) );
}

//---------------------------------------------------------------------------//
// A simple converter to extract the numbers from the value and assign them to
// ascending elements of the DynArray<int>.  This is very rough.  Valid input
// should give valid output, but invalid input is anyone's guess.
//---------------------------------------------------------------------------//

void nmi_DynArray_int::update_value()
{
    DynArray<int>& d = *pdai;

// Code to parse new value;

    DynArray<int> r(1);
    int i=0, k=0;

    while( value[i] ) {

    // Chop of leader, whitespace, etc, if any.

	while( isspace(value[i]) || value[i] == '{' || value[i] == '[' )
	    i++;
	if (!value[i] || value[i] == '{' || value[i] == ']' ) break;

    // Must be a value, so get it.

	String num;
	while( isdigit(value[i]) ) {
	    num += value[i];
	    i++;
	}
	r[k++] = atoi( &num[0] );

	if (value[i] == ',') i++;
    }

    d = r;
}

//---------------------------------------------------------------------------//
// Produce the canonical form for representing a DynArray<int> on disk.
//---------------------------------------------------------------------------//

String nmi_DynArray_int::File_Rep()
{
    String v = name + "={";
    DynArray<int>& d = *pdai;
    char buf[15];

    for( int i = d.low(); i <= d.high(); i++ ) {
	if (i > d.low()) v += ',';
	sprintf( buf, "%d", d[i] );
	v += buf;
    }

    v += '}';
    return v;
}

//---------------------------------------------------------------------------//
// A very simple parser for inputing a DynArray<int>.  Valid forms are:
//   var = #,#,#,#     whitespace delimited (leading space optional)
//   var = {#,#,#}     where {} can also be "" or [] and , can also be space.
//---------------------------------------------------------------------------//

void nmi_DynArray_int::Read_Value( nmistream& nis )
{
    String v;
    char c;

// We start on the =.  Skip over any whitespace.

    while( nis.get(c) ) {
	if (isspace(c)) continue;
	break;
    }

// Select the delimiter.

    int delimited = 0;
    char delim = ' ';
    if (c == '\"')
	delim = '\"', delimited=1, nis.get(c);

    else if (c == '{')
	delim = '}', delimited=1, nis.get(c);

    else if (c == '[')
	delim = ']', delimited=1, nis.get(c);

// Absorb to closing delimiter.

    while( c &&
	   (delimited && c != delim) || (!delimited && !isspace(c)) ) {
	v += c;
	nis.get(c);
    }

// If whitespace delimited, leave it in the stream.

    if (c == ' ') nis.putback(c);

    // Complete the read.

    set_value( v );
    update_value();
    do_callbacks();
}

//===========================================================================//
// class nmi_PageBreak

// Format control.
//===========================================================================//

nmi_PageBreak::nmi_PageBreak( va_list *pap )
{
    name = "PageBreak";

    char *descrip = va_arg( *pap, char * );

    value = default_value = descrip;

}

void nmi_PageBreak::modify_default( va_list *pap )
{
// Meaningless.
}

void nmi_PageBreak::update_value()
{
// Meaningless.
}

String nmi_PageBreak::File_Rep()
{
    String r = "\n\t# " + value + "\n\t";

    return r;
}

//===========================================================================//
// class nmi_BoxStart

// Format control.
//===========================================================================//

nmi_BoxStart::nmi_BoxStart( va_list *pap )
{
    name = "BoxStart";
    default_value = "";

// Do nothing for now.  Later we might try to add ...
}

void nmi_BoxStart::modify_default( va_list *pap )
{
// Meaningless.
}

void nmi_BoxStart::update_value()
{
// Meaningless.
}

String nmi_BoxStart::File_Rep()
{
    String r = "\n\t";

    return r;
}

//===========================================================================//
// class nmi_BoxEnd

// Format control.
//===========================================================================//

nmi_BoxEnd::nmi_BoxEnd( va_list *pap )
{
    name = "BoxEnd";
    default_value = "";

// Do nothing for now.  Later we might try to add ...
}

void nmi_BoxEnd::modify_default( va_list *pap )
{
// Meaningless.
}

void nmi_BoxEnd::update_value()
{
// Meaningless.
}

String nmi_BoxEnd::File_Rep()
{
    String r = "\n\t";

    return r;
}

//---------------------------------------------------------------------------//
//                              end of Items.cc
//---------------------------------------------------------------------------//
