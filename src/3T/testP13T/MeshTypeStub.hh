//----------------------------------*-C++-*----------------------------------//
// MeshTypeStub.hh
// Randy M. Roberts
// Fri Mar 20 13:40:29 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_MeshTypeStub_hh__
#define __3T_testP13T_MeshTypeStub_hh__

#include "ds++/SP.hh"

#include <iostream>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class MeshTypeStub - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class MeshTypeStub
{

    // NESTED CLASSES AND TYPEDEFS

  public:
    
    struct scalar
    {
	double val;
	SP<MeshTypeStub> spmesh;
	
	scalar(const SP<MeshTypeStub> &spmesh_)
	    : spmesh(spmesh_), val(0.0)
	{
	    // empty
	}

	scalar &operator=(double rhs) { val = rhs; return *this;}

	int size() const { return 1; }
    };
    
    typedef scalar ccsf;
    typedef scalar bsbf;
    typedef scalar ncvf;
    typedef scalar fcdsf;

    // DATA
    
  public:

    // CREATORS
    
    // MANIPULATORS
    
    // ACCESSORS

    void print(std::ostream &os) const
    {
	os << "in MeshTypeStub::print(), this: "
	   << (void *)(this);
    }

  private:
    
    // IMPLEMENTATION
};

MeshTypeStub::scalar operator+(const MeshTypeStub::scalar &lop,
			       const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar operator+(const MeshTypeStub::scalar &lop, double rop);
MeshTypeStub::scalar operator+(double lop, const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar &operator+=(MeshTypeStub::scalar &lop,
				 const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar operator-(const MeshTypeStub::scalar &lop,
			       const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar operator-(const MeshTypeStub::scalar &lop, double rop);
MeshTypeStub::scalar operator-(double lop, const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar &operator-=(MeshTypeStub::scalar &lop, double rop);
MeshTypeStub::scalar operator*(const MeshTypeStub::scalar &lop,
			       const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar operator*(const MeshTypeStub::scalar &lop, double rop);
MeshTypeStub::scalar operator*(double lop, const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar &operator*=(MeshTypeStub::scalar &lop, double rop);
MeshTypeStub::scalar operator/(const MeshTypeStub::scalar &lop,
			       const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar operator/(const MeshTypeStub::scalar &lop, double rop);
MeshTypeStub::scalar operator/(double lop, const MeshTypeStub::scalar &rop);
MeshTypeStub::scalar &operator/=(MeshTypeStub::scalar &lop, double rop);

inline std::ostream &operator<<(std::ostream &os,
				const MeshTypeStub::scalar &rop)
{
    return os << rop.val;
}

END_NS_XTM  // namespace XTM

#endif                          // __3T_testP13T_MeshTypeStub_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/MeshTypeStub.hh
//---------------------------------------------------------------------------//
