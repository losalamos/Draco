//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/Solver.hh
 * \author Randy M. Roberts
 * \date   Tue Jan 25 14:05:44 2000
 * \brief  The interface to the LAMG F90 Solver.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMG_Solver_hh__
#define __LAMG_Solver_hh__

#include <LAMG/config.h>
#include "Options.hh"
#include "CompressedRowStorage.hh"
#include "ds++/SP.hh"
#include "ds++/Mat.hh"
#include "traits/ContainerTraits.hh"

#include <algorithm>

namespace rtt_LAMG
{

//===========================================================================//
/*!
 * \class Solver
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Solver 
{

  public:
    
    // NESTED CLASSES AND TYPEDEFS

    // We keep the representation in a separate class/contained object
    // so that we may pass LamgMatrixDcsrR's around be value
    // without having to create and destroy stuff on the F90 side.
    
    class Representation;

    typedef RTT_F90_INTEGER F90Pointer;

  public:

    typedef rtt_dsxx::Mat1<RTT_F90_INTEGER> IMat1;
    typedef rtt_dsxx::Mat1<RTT_F90_DOUBLE> DMat1;

    // DATA

  private:
    
    rtt_dsxx::SP<Representation> rep;
    
  public:

    // CREATORS
    
    inline Solver(const Options &opts);
    
    //DEFAULT: Solver(const Solver &rhs);
    //DEFAULT: ~Solver();

    // MANIPULATORS

    //DEFAULT: Solver& operator=(const Solver &rhs);

    void setOptions(const Options &opts);

    void solve(DMat1 &x, const CompressedRowStorage &crs, const DMat1 &b);

    template<class FT>
    void solve(FT &x, const CompressedRowStorage &crs, const FT &b)
    {
	typedef rtt_traits::ContainerTraits<FT> CT;
	Require(CT::conformal(x, b));

	DMat1 m1x(x.size());
	std::copy(x.begin(), x.end(), m1x.begin());

	DMat1 m1b(b.size());
	std::copy(b.begin(), b.end(), m1b.begin());

	solve(m1x, crs, m1b);

	std::copy(m1x.begin(), m1x.end(), x.begin());
    }
    
    // ACCESSORS

  private:
    
    // IMPLEMENTATION

    inline F90Pointer opaquePointer() const; 
};

class Solver::Representation
{
    // NESTED CLASSES AND TYPEDEFS

    // DATA

    F90Pointer ptr;

  public:

    // CREATORS

    Representation(const Options &options_in);
    ~Representation();
    
    // ACCESSORS

    F90Pointer opaquePointer() const { return ptr; }
};

Solver::Solver(const Options &options_in)
    : rep(new Representation(options_in))
{
    // empty
}

Solver::F90Pointer Solver::opaquePointer() const
{
    return rep->opaquePointer();
}

} // end namespace rtt_LAMG

#endif                          // __LAMG_Solver_hh__

//---------------------------------------------------------------------------//
//                              end of LAMG/Solver.hh
//---------------------------------------------------------------------------//
