#include "3T/testP13T/testFullP13T.hh"
#include "3T/P13T.cc"

using XTM::testFullP13T;

typedef rtt_matprops::MarshakMaterialProps UMCMP;
typedef testFullP13T<UMCMP>::MT MT;
typedef testFullP13T<UMCMP>::DS DS;
typedef testFullP13T<UMCMP>::MatStateFC MSFFC;
typedef testFullP13T<UMCMP>::MatStateCC MSFCC;

template class P13T<MT,MSFCC,MSFFC,DS>;


typedef rtt_matprops::InterpedMaterialProps UMCMP2;
typedef testFullP13T<UMCMP2>::MT MT2;
typedef testFullP13T<UMCMP2>::DS DS2;
typedef testFullP13T<UMCMP2>::MatStateFC MSFFC2;
typedef testFullP13T<UMCMP2>::MatStateCC MSFCC2;

template class P13T<MT2,MSFCC2,MSFFC2,DS2>;
