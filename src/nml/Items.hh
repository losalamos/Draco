//----------------------------------*-C++-*----------------------------------//
// Items.hh
// Geoffrey Furnish
// 25 January 1994
//---------------------------------------------------------------------------//
// @> Declare a variety of useful basic Item types.
//---------------------------------------------------------------------------//

#ifndef __nml_Items_hh__
#define __nml_Items_hh__

#include <map>
// using std::map;

#include "Item.hh"
#include "ds++/DynArray.hh"

#include <stdarg.h>

#define rUNDEFINED -999999.
#define iUNDEFINED -999999

#define NMI_INT         1
#define NMI_FLOAT       2
#define NMI_DOUBLE      3
#define NMI_INT_SET     4
#define NMI_INT_ARRAY   5
#define NMI_FLOAT_ARRAY 6
#define NMI_BOOL        7
#define NMI_FLOG        8
#define NMI_STRING      9
#define NMI_INT_DYNARRAY 10

// Some special cases to help with formatting for groups with many items.

#define NMI_PAGE_BREAK 50
#define NMI_BOX_START  51
#define NMI_BOX_END    52

#define NMI_ADD_END   100
#define NMI_END       100

#define NMI_MODIFY_END "nmi_modify_block_defaults_end"
#define NMI_MOD_END    NMI_MODIFY_END

#define N_int(a,b) NMI_INT, #a, b, &a

class nmi_int : public NML_Item {

    int *pv;

  public:
    nmi_int( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String Type() const { return "int"; }
};

#define N_float(a,b) NMI_FLOAT, #a, b, &a

class nmi_float : public NML_Item {

    float *pv;

  public:
    nmi_float( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String Type() const { return "float"; }
};

#define N_double(a,b) NMI_DOUBLE, #a, b, &a

class nmi_double : public NML_Item {

    double *pv;

  public:
    nmi_double( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String Type() const { return "double"; }
};

class nmi_int_set : public NML_Item {

  protected:
    int *pv;
    std::map<rtt_dsxx::String, int> m;

    nmi_int_set( int undef, int nvals, ... );

  public:
    nmi_int_set( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String Type() const { return "int_set"; }
    virtual rtt_dsxx::String Choices() const;
};

class nmi_int_array : public NML_Item {

    int *pv;			// Pointer to storage array of length
    int max_vals;		// max_vals.
    int *cur_vals;		// # of vals user actually specified.

  public:
    nmi_int_array( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String Type() const { return "int_array"; }
};

#define N_bool(a,b) NMI_BOOL, #a, #b, &a

class nmi_bool : public nmi_int_set {

  public:
    nmi_bool( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String Type() const { return "bool"; }
};

// This is known to be right for Linux, HP/UX, AIX, SUNOS, and UNICOS
// on Cray 2 series machines.  We\'ll assume it works on all other
// machines unless we know specifically otherwise.

#ifndef FLOG_TRUE
#define FLOG_TRUE 1
#endif

#ifndef FLOG_FALSE
#define FLOG_FALSE 0
#endif

// Cray Y-MP series machines are "different" :-).

#ifdef CRAY1
#undef FLOG_TRUE
#undef FLOG_FALSE
#define FLOG_TRUE -1
#define FLOG_FALSE 0
#endif

#define N_flog(a,b) NMI_FLOG, #a, #b, &a

class nmi_flog : public nmi_int_set {

  public:
    nmi_flog( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String Type() const { return "flog"; }
};

#define N_String(a,b) NMI_STRING, #a, b, &a

class nmi_String : public NML_Item {

    rtt_dsxx::String *ps;

  public:
    nmi_String( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String File_Rep();
    void Read_Value( nmistream& nis );

    rtt_dsxx::String Type() const { return "String"; }
};

#define N_DynArray_int(a,b) NMI_INT_DYNARRAY, #a, b, &a

class nmi_DynArray_int : public NML_Item {

    rtt_dsxx::DynArray<int> *pdai;

  public:
    nmi_DynArray_int( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();

    rtt_dsxx::String File_Rep();
    void Read_Value( nmistream& nis );

    rtt_dsxx::String Type() const { return "DynArray<int>"; }
};

// The page formatting item types.

class nmi_PageBreak : public NML_Item {
  public:
    nmi_PageBreak( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();
    rtt_dsxx::String File_Rep();

    rtt_dsxx::String Type() const { return "PageBreak"; }
};

class nmi_BoxStart : public NML_Item {
  public:
    nmi_BoxStart( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();
    rtt_dsxx::String File_Rep();

    rtt_dsxx::String Type() const { return "BoxStart"; }
};

class nmi_BoxEnd : public NML_Item {
  public:
    nmi_BoxEnd( va_list *pap );
    void modify_default( va_list *pap );
    void update_value();
    rtt_dsxx::String File_Rep();

    rtt_dsxx::String Type() const { return "BoxEnd"; }
};

#endif				// __nml_Items_hh__

//---------------------------------------------------------------------------//
//                              end of Items.hh
//---------------------------------------------------------------------------//
