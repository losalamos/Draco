//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/test/TestConjGradMatVec.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 10:25:10 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_test_TestConjGradMatVec_hh__
#define __ConjGrad_test_TestConjGradMatVec_hh__

#include "ds++/Mat.hh"
#include "ds++/Assert.hh"

namespace rtt_ConjGrad_test
{
 
//===========================================================================//
/*!
 * \class TestConjGradMatVec
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class T>
class TestConjGradMatVec 
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef T value_type;
    typedef rtt_dsxx::Mat1<T> Mat;

    // Forward Reference
    class JacobiPrecon;

    friend class JacobiPrecon;

  private:

    // DATA

    int nxs;
    int nys;
    
  public:

    // CREATORS
    
    TestConjGradMatVec(int nxs_in, int nys_in)
	: nxs(nxs_in), nys(nys_in)
    {
	// empty
    }
    
    //Defaulted: TestConjGradMatVec(const TestConjGradMatVec &rhs);
    //Defaulted: ~TestConjGradMatVec();

    // MANIPULATORS
    
    //Defaulted: TestConjGradMatVec& operator=(const TestConjGradMatVec &rhs);

    // ACCESSORS

    Mat operator()(const Mat &x) const
    {
	Assert(x.size() == nxs*nys);
	
	Mat b(x);
	
	for( int ix = 0; ix < nxs; ix++ )
	{
	    for( int iy = 0; iy < nys; iy++ )
	    {
		int indva   = ix     + nxs*(iy);
		int indvqr0 = ix     + nxs*(iy);
		int indvqr1 = ix - 1 + nxs*(iy);
		int indvqr2 = ix + 1 + nxs*(iy);
		int indvqr3 = ix     + nxs*(iy-1);
		int indvqr4 = ix     + nxs*(iy+1);

		b(indva) = 4.0*x(indvqr0);

		if( ix !=     0 ) b(indva) -= x(indvqr1);
		if( ix != nxs-1 ) b(indva) -= x(indvqr2);
		if( iy !=     0 ) b(indva) -= x(indvqr3);
		if( iy != nys-1 ) b(indva) -= x(indvqr4);
	    }
	}
	return b;
    }
	

  private:
    
    // IMPLEMENTATION
};

template<class T>
class TestConjGradMatVec<T>::JacobiPrecon
{
    typedef TestConjGradMatVec<T>::Mat Mat;
    
    int nxs;
    int nys;

  public:
    
    JacobiPrecon(const TestConjGradMatVec<T> &matvec)
	: nxs(matvec.nxs), nys(matvec.nys)
    {
	// empty
    }

    Mat operator()(const Mat &x) const
    {
	Require(x.size() == nxs*nys);

	Mat b(x);
	for (Mat::iterator bit = b.begin(); bit != b.end(); ++bit)
	    (*bit) /= 4.0;
	return b;
    }
    
};

} // end namespace rtt_ConjGrad_test

#endif                          // __ConjGrad_test_TestConjGradMatVec_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/test/TestConjGradMatVec.hh
//---------------------------------------------------------------------------//
