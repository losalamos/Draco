//----------------------------------*-C++-*----------------------------------//
// Sort.pt
// Dave Nystrom
// 13 January 1997
//---------------------------------------------------------------------------//

#include "Sort.cc"

template void HeapSort( int, int * );
template void index( int n, int *v, int *idx );

template void HeapSort( int, float * );
template void index( int n, float *v, int *idx );

template void HeapSort( int, double * );
template void index( int n, double *v, int *idx );
