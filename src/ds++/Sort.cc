//----------------------------------*-C++-*----------------------------------//
// Sort.cc
// Geoffrey Furnish
// Thu Sep 28 09:31:08 1995
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Sort.hh"

#include <iostream>

using namespace std;

// Implement the various sort functions.

// Heapsort, ala NR.  (Was SORT.C)

template<class T>
void HeapSort( int n, T *ra )
{
    int l,j,ir,i;
    T rra;

    l=(n >> 1)+1;
    ir=n;
    for (;;) {
	if (l > 1)
	    rra=ra[--l];
	else {
	    rra=ra[ir];
	    ra[ir]=ra[1];
	    if (--ir == 1) {
		ra[1]=rra;
		return;
	    }
	}
	i=l;
	j=l << 1;
	while (j <= ir) {
	    if (j < ir && ra[j] < ra[j+1]) ++j;
	    if (rra < ra[j]) {
		ra[i]=ra[j];
		j += (i=j);
	    }
	    else j=ir+1;
	}
	ra[i]=rra;
    }
}

// Produce an ordered index for an array of T's.  (From NR, INDEX.C)

template<class T>
void index( int n, T *v, int *idx )
// void indexx(n,arrin,indx)
// int n,indx[];
// float arrin[];
{
// The sucky NR code assumes arrays based at 1, sheesh.  
// Hack it like this (snort):

    T *arrin = v-1;
    int *indx = idx-1;

    int l, j, ir, indxt, i;
    T q;

    for( j=1; j<=n; j++ )
	indx[j] = j;

    l=(n >> 1) + 1;
    ir=n;
    for (;;) {
	if (l > 1)
	    q=arrin[(indxt=indx[--l])];
	else {
	    q=arrin[(indxt=indx[ir])];
	    indx[ir]=indx[1];
	    if (--ir == 1) {
		indx[1]=indxt;
	    // return; Oh gag, this is sucky, can't leave here, the indexes
	    //         are all computed to be 1 based.  Sheesh!
		break;
	    }
	}
	i=l;
	j=l << 1;
	while (j <= ir) {
	    if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
	    if (q < arrin[indx[j]]) {
		indx[i]=indx[j];
		j += (i=j);
	    }
	    else j=ir+1;
	}
	indx[i]=indxt;
    }

// Okay, done sorting the indexes, now fix them for conventional C style
// indexing.  grrrrrrrrrrrrrrr.

    for( j=1; j <= n; j++ )
	indx[j]--;
}

//---------------------------------------------------------------------------//
//                              end of Sort.cc
//---------------------------------------------------------------------------//
