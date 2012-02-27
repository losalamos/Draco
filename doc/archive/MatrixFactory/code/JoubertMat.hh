/*-----------------------------------*-C++-*---------------------------------*/
/* JoubertMat.hh */
/* Randy M. Roberts */
/* Thu May 27 15:48:02 1999 */
/*---------------------------------------------------------------------------*/
/* @> */
/*---------------------------------------------------------------------------*/

#ifndef __MatrixFactory_JoubertMat_hh__
#define __MatrixFactory_JoubertMat_hh__

#include <vector>
#include <iostream>

#include "ds++/Assert.hh"

namespace JoubertNS
{

typedef double F90Real;
typedef int F90Int;

class JoubertMat
{
    F90Int nrow;
    F90Int ncol;
    F90Int nnz;
    std::vector<F90Int> ia;
    std::vector<F90Int> ja;
    std::vector<F90Real> a;

    friend std::ostream &operator<<(std::ostream &os, const JoubertMat &mat)
    {
	return mat.print(os);
    }
    
  public:

    JoubertMat(int nrow_, int ncol_, int nnz_,
	       const std::vector<F90Int> &ia_,
	       const std::vector<F90Int> &ja_,
	       const std::vector<F90Real> &a_)
	: nrow(nrow_), ncol(ncol_), nnz(nnz_),
	  ia(ia_), ja(ja_), a(a_)
    {
	Require(ia.size() == nrow+1);
	Require(ja.size() == nnz);
	Require(a.size() == nnz);
    }

    // ~JoubertMat()
    // {
    //    std::cerr << "In JoubertNS::~JoubertMat" << std::endl;
    // }
    
    std::ostream &print(std::ostream &os) const
    {
	os << " nrows: " << nrow << std::endl
	   << " ncols: " << ncol << std::endl
	   << " nnz: " << nnz << std::endl;

	os << "ia: ";
	std::copy(ia.begin(), ia.end(), std::ostream_iterator<int>(os, " "));
	os << std::endl;

	for (int ir=1; ir<=nrow; ir++)
	{
	    for (int iic=ia[ir-1]; iic<ia[ir]; iic++)
	    {
		int ic = ja[iic-1];
		double data = a[iic-1];

		os << "iic: " << iic << ", ";
		os << "A(" << ir << "," << ic << "): " << data << std::endl;
	    }
	}

	os << "Done" << std::endl;
	return os;
    }
};

} // namespace JoubertNS

#endif    /* __MatrixFactory_JoubertMat_hh__ */

/*---------------------------------------------------------------------------*/
/*    end of JoubertMat.hh */
/*---------------------------------------------------------------------------*/
