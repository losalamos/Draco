//---------------------------------*-C++-*---------------------------------//
// random.cc
// Geoffrey Furnish
// 7 May 1992
//-------------------------------------------------------------------------//
// @> Random number generation functions from Numerical Recipes.
//-------------------------------------------------------------------------//

#include <math.h>
#ifndef m88k
#include <stdlib.h>
#endif

#include "util/random.hh"

float ran0(int *idum)
{
    static float y,maxran,v[98];
    static int iff=0;
    int j;
    unsigned i,k;

    if (*idum < 0 || iff == 0) {
	iff=1;
	i=2;
	do {
	    k=i;
	    i<<=1;
	} while (i);
	maxran=k;
	srand(*idum);
	*idum=1;
	for (j=1;j<=97;j++)
	    rand();
	for (j=1;j<=97;j++)
	    v[j]=rand();
	y=rand();
    }
    j=int( 1+97.0*y/maxran );
//	if (j > 97 || j < 1) nrerror("RAN0: This cannot happen.");
    y=v[j];
    v[j]=rand();
    return y/maxran;
}

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

float ran1(int *idum)
{
    static long ix1,ix2,ix3;
    static float r[98];
    float temp;
    static int iff=0;
    int j;

    if (*idum < 0 || iff == 0) {
	iff=1;
	ix1=(IC1-(*idum)) % M1;
	ix1=(IA1*ix1+IC1) % M1;
	ix2=ix1 % M2;
	ix1=(IA1*ix1+IC1) % M1;
	ix3=ix1 % M3;
	for (j=1;j<=97;j++) {
	    ix1=(IA1*ix1+IC1) % M1;
	    ix2=(IA2*ix2+IC2) % M2;
	    r[j]=(ix1+ix2*RM2)*RM1;
	}
	*idum=1;
    }
    ix1=(IA1*ix1+IC1) % M1;
    ix2=(IA2*ix2+IC2) % M2;
    ix3=(IA3*ix3+IC3) % M3;
    j=int( 1 + ((97*ix3)/M3) );
//	if (j > 97 || j < 1) nrerror("RAN1: This cannot happen.");
    temp=r[j];
    r[j]=(ix1+ix2*RM2)*RM1;
    return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

#define M 714025
#define IA 1366
#define IC 150889

float ran2(long *idum)
{
    static long iy,ir[98];
    static int iff=0;
    int j;

    if (*idum < 0 || iff == 0) {
	iff=1;
	if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
	for (j=1;j<=97;j++) {
	    *idum=(IA*(*idum)+IC) % M;
	    ir[j]=(*idum);
	}
	*idum=(IA*(*idum)+IC) % M;
	iy=(*idum);
    }
    j=int( 1 + 97.0*iy/M );
//	if (j > 97 || j < 1) nrerror("RAN2: This cannot happen.");
    iy=ir[j];
    *idum=(IA*(*idum)+IC) % M;
    ir[j]=(*idum);
    return (float) iy/M;
}

#undef M
#undef IA
#undef IC

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(int *idum)
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;

    if (*idum < 0 || iff == 0) {
	iff=1;
	mj=MSEED-(*idum < 0 ? -*idum : *idum);
	mj %= MBIG;
	ma[55]=mj;
	mk=1;
	for (i=1;i<=54;i++) {
	    ii=(21*i) % 55;
	    ma[ii]=mk;
	    mk=mj-mk;
	    if (mk < MZ) mk += MBIG;
	    mj=ma[ii];
	}
	for (k=1;k<=4;k++)
	    for (i=1;i<=55;i++) {
		ma[i] -= ma[1+(i+30) % 55];
		if (ma[i] < MZ) ma[i] += MBIG;
	    }
	inext=0;
	inextp=31;
	*idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;
    return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

float gasdev(int *idum)
{
    static int iset=0;
    static float gset;
    float fac,r,v1,v2;
    
    if  (iset == 0) {
	do {
	    v1=2.0*ran1(idum)-1.0;
	    v2=2.0*ran1(idum)-1.0;
	    r=v1*v1+v2*v2;
	} while (r >= 1.0);
	fac=sqrt(-2.0*log(r)/r);
	gset=v1*fac;
	iset=1;
	return v2*fac;
    } else {
	iset=0;
	return gset;
    }
}

//-------------------------------------------------------------------------//
//                              end of random.cc
//-------------------------------------------------------------------------//
