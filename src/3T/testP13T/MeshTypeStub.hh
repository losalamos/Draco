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
#include "traits/ContainerTraits.hh"

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
	typedef const double *const_iterator;
	typedef double *iterator;
	typedef double value_type;

	double val;
	SP<MeshTypeStub> spmesh;
	
	scalar(const SP<MeshTypeStub> &spmesh_, double val_=0.0)
	    : spmesh(spmesh_), val(val_)
	{
	    // empty
	}

	scalar &operator=(double rhs) { val = rhs; return *this;}

	int size() const { return 1; }

	iterator begin()
	{
	    return &val;
	}
	const_iterator begin() const
	{
	    return &val;
	}
	iterator end()
	{
	    return &val+1;
	}
	const_iterator end() const
	{
	    return &val+1;
	}
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

    std::ostream &print(std::ostream &os) const
    {
	os << "(MeshTypeStub::this: " << (void *)(this) << ")";
	return os;
    }

  private:
    
    // IMPLEMENTATION
};

inline std::ostream &operator<<(std::ostream &os,
				const MeshTypeStub &rhs)
{
    return rhs.print(os);
}

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

template <>
class ContainerTraits<MeshTypeStub::scalar>
{
  public:
    typedef const double *const_iterator;
    typedef double *iterator;
    static inline iterator begin(MeshTypeStub::scalar &a)
    {
	return a.begin();
    }
    static inline const_iterator begin(const MeshTypeStub::scalar &a)
    {
	return a.begin();
    }
    static inline iterator end(MeshTypeStub::scalar &a)
    {
	return a.end();
    }
    static inline const_iterator end(const MeshTypeStub::scalar &a)
    {
	return a.end();
    }
    static inline bool conformal(MeshTypeStub::scalar a, MeshTypeStub::scalar b)
    {
	return true;
    }
};

END_NS_XTM  // namespace XTM

#endif                          // __3T_testP13T_MeshTypeStub_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/MeshTypeStub.hh
//---------------------------------------------------------------------------//
