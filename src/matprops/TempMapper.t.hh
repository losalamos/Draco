//----------------------------------*-C++-*----------------------------------//
// TempMapper.t.cc
// Randy M. Roberts
// Thu Oct  8 09:44:09 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "TempMapper.hh"

#include "traits/ContainerTraits.hh"
#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_matprops
{

 template<class MT>
 TempMapper<MT>::TempMapper(double power_, double gamma_,
                            const SP<MT> &spMesh_,
			    const FieldConstructor &fCtor_)
    : power(power_), gamma(gamma_), spMesh(spMesh_), fCtor(fCtor_)
 {
     Require(gamma >= 0.0 && gamma <= 1.0);
     Require(power != 0.0);
 }

 template<class MT>
 template<class FCV, class CCV, class OpAssignV>
 void TempMapper<MT>::tempCC2FC(FCV &faceTempsByMat, const CCV &cellTempsByMat,
				const typename MT::ccsf &cellTempsAvg,
				const OpAssignV &opAssignV) const
 {
     // Copy the average CC temperatures to the cell faces.
     
     MT::fcdsf faceTempsAvg(fCtor);
     MT::gather(faceTempsAvg, cellTempsAvg, MT::OpAssign());

     // Calculate the multiplicative ratio value "alpha"

     MT::fcdsf alpha(fCtor);

     calculateAlpha(alpha, faceTempsAvg);

     // Copy the cell center material temps to the faces

     MT::gather(faceTempsByMat, cellTempsByMat, opAssignV);

     FCV::iterator ftmit = faceTempsByMat.begin();
     MT::fcdsf::const_iterator ait = alpha.begin();

     while (ftmit != faceTempsByMat.end())
     {
	 double alpha = *ait;

	 typedef rtt_traits::ContainerTraits<FCV::value_type> CT;

	 for (CT::iterator tempit = CT::begin(*ftmit);
	      tempit != CT::end(*ftmit); tempit++)
	 {
	     *tempit *= alpha;
	 }

	 ftmit++;
	 ait++;
     }
 }

 template<class MT>
 void TempMapper<MT>::calculateAlpha(typename MT::fcdsf &alpha,
				     const typename MT::fcdsf &faceTempsAvg)
     const
 {
     // Get the temperatures in the cell adjacent to each face

     MT::fcdsf adjacentFaceTemp(fCtor);
     MT::swap_faces(adjacentFaceTemp, faceTempsAvg);

     MT::fcdsf::const_iterator aFTit = adjacentFaceTemp.begin();
     MT::fcdsf::const_iterator fTit = faceTempsAvg.begin();
     MT::fcdsf::iterator ait = alpha.begin();

     while (aFTit != adjacentFaceTemp.end())
     {
	 double Tavg = *fTit;
	 double Thigh = std::max(*aFTit, *fTit);
	 double Tlow = std::min(*aFTit, *fTit);

	 double &alpha = *ait;

	 using std::pow;
	 alpha = gamma*pow(Thigh, power) + (1.0 - gamma)*pow(Tlow, power);
	 alpha = pow(alpha, 1.0/power);
	 alpha /= Tavg;

	 aFTit++;
	 fTit++;
	 ait++;
     }
 }
}

//---------------------------------------------------------------------------//
//                              end of TempMapper.t.cc
//---------------------------------------------------------------------------//
