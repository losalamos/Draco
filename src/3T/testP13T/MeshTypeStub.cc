//----------------------------------*-C++-*----------------------------------//
// MeshTypeStub.cc
// Randy M. Roberts
// Fri Mar 20 13:40:29 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/MeshTypeStub.hh"

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

MeshTypeStub::scalar operator+(const MeshTypeStub::scalar &lop,
			       double rop)
{
    MeshTypeStub::scalar sf(lop);
    sf.val += rop;
    return sf;
}
MeshTypeStub::scalar operator+(double lop,
			       const MeshTypeStub::scalar &rop)
{
    MeshTypeStub::scalar sf(rop);
    sf.val += lop;
    return sf;
}
    
MeshTypeStub::scalar operator+(const MeshTypeStub::scalar &lop,
			       const MeshTypeStub::scalar &rop)
{
    MeshTypeStub::scalar sf(lop);
    sf.val += rop.val;
    return sf;
}

MeshTypeStub::scalar &operator+=(MeshTypeStub::scalar &lop,
				 const MeshTypeStub::scalar &rop)
{
    lop.val += rop.val;
    return lop;
}

MeshTypeStub::scalar operator-(const MeshTypeStub::scalar &lop,
			       double rop)
{
    MeshTypeStub::scalar sf(lop);
    sf.val -= rop;
    return sf;
}
MeshTypeStub::scalar operator-(double lop,
			       const MeshTypeStub::scalar &rop)
{
    MeshTypeStub::scalar sf(rop.spmesh);
    sf.val = lop - rop.val;
    return sf;
}
    
MeshTypeStub::scalar operator-(const MeshTypeStub::scalar &lop,
			       const MeshTypeStub::scalar &rop)
{
    MeshTypeStub::scalar sf(lop);
    sf.val -= rop.val;
    return sf;
}
    
MeshTypeStub::scalar operator*(const MeshTypeStub::scalar &lop,
			       const MeshTypeStub::scalar &rop)
{
    MeshTypeStub::scalar sf(lop);
    sf.val *= rop.val;
    return sf;
}

MeshTypeStub::scalar operator*(const MeshTypeStub::scalar &lop,
			       double rop)
{
    MeshTypeStub::scalar sf(lop);
    sf.val *= rop;
    return sf;
}

MeshTypeStub::scalar operator*(double lop,
			       const MeshTypeStub::scalar &rop)
{
    MeshTypeStub::scalar sf(rop);
    sf.val *= lop;
    return sf;
}

MeshTypeStub::scalar &operator*=(MeshTypeStub::scalar &lop,
				 double rop)
{
    lop.val *= rop;
    return lop;
}

MeshTypeStub::scalar operator/(const MeshTypeStub::scalar &lop,
			       const MeshTypeStub::scalar &rop)
{
    MeshTypeStub::scalar sf(lop);
    sf.val /= rop.val;
    return sf;
}

MeshTypeStub::scalar operator/(const MeshTypeStub::scalar &lop,
			       double rop)
{
    MeshTypeStub::scalar sf(lop);
    sf.val /= rop;
    return sf;
}

MeshTypeStub::scalar operator/(double lop,
			       const MeshTypeStub::scalar &rop)
{
    MeshTypeStub::scalar sf(rop.spmesh);
    sf.val = lop / rop.val;
    return sf;
}
    
MeshTypeStub::scalar &operator/=(MeshTypeStub::scalar &lop,
				 double rop)
{
    lop.val /= rop;
    return lop;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of MeshTypeStub.cc
//---------------------------------------------------------------------------//
