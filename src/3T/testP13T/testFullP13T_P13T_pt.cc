#include "3T/testP13T/testFullP13T.hh"
#include "3T/P13T.cc"

using XTM::testFullP13T;

typedef testFullP13T::MT MT;
typedef testFullP13T::MP MP;
typedef testFullP13T::DS DS;

template class P13T<MT,MP,DS>;
