//----------------------------------*-C++-*----------------------------------//
// MultiMatCellMatProps.t.cc
// John McGhee
// Mon Sep 14 08:38:38 1998
//---------------------------------------------------------------------------//
// @> A material property class for cells which contain more than one material.
//---------------------------------------------------------------------------//

#include "matprops/MultiMatCellMatProps.hh"
#include "ds++/Assert.hh"
#include <algorithm>

namespace rtt_matprops 
{

 // CREATORS

 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 MultiMatCellMatProps<UMCMP>::MaterialStateField<FT, FT1, FT2>::
 MaterialStateField( const MultiMatCellMatProps<UMCMP> &matprops_,
		     const FT1 &density_, const FT1 &electronTemp_,
		     const FT1 &ionTemp_, const FT1 &volumeFraction_,
		     const FT2 &matId_)
{
    Require(density_.size() > 0);
    // Check to be sure all input has the same cell count.
    Require(density_.size() == electronTemp_.size());
    Require(density_.size() == ionTemp_.size());
    Require(density_.size() == volumeFraction_.size());
    Require(density_.size() == matId_.size());
    // Loop over cells.

    FT1::const_iterator etit  = electronTemp_.begin();
    FT1::const_iterator itit  = ionTemp_.begin();
    FT1::const_iterator vfit  = volumeFraction_.begin();
    FT2::const_iterator midit = matId_.begin();
    for (FT1::const_iterator denit = density_.begin(); 
	 denit != density_.end(); denit++, etit++, 
	     itit++, vfit++, midit++)
    {
	Require( (*denit).size() > 0 );
	using std::vector;
	typedef FT2::value_type::value_type value_type2;

	vector<value_type1> density((*denit).begin(), (*denit).end());
	vector<value_type1> electronTemp((*etit).begin(), (*etit).end());
	vector<value_type1> ionTemp((*itit).begin(), (*itit).end());
	vector<value_type2> matId((*midit).begin(), (*midit).end());

	
	// Load a material state field for this cell.
	matState.push_back(matprops_.spumcmp->getMaterialState(
	    density, electronTemp, ionTemp, matId));
	// Check to be sure that the input for this cell has a 
	// material count consistent with the other input.
	Require( (*vfit).size() == (*denit).size() );
	// Loop over materials, checking the volume fraction for
	// positivity and accumulating a normalization factor.
	FT1::value_type::const_iterator vfmait = (*vfit).begin();
	value_type1 xsum = 0.;
	int nmat = (*vfit).size();
	for (int imat = 0; imat < nmat; imat++, vfmait++)
	{
	    Require(*vfmait >= 0.);
	    xsum += *vfmait;
	}
	Require(xsum > 0.);
	xsum = 1./xsum;
	// Loop over materials, normalizing the input material volumes.
	std::vector<value_type1> voltmp(nmat);
	vfmait = (*vfit).begin();
	for (int imat = 0; imat < nmat; imat++, vfmait++)
	{
	    voltmp[imat] = (*vfmait)*xsum;
	}
	// Load the volume fractions for this cell.
	volumeFraction.push_back(voltmp);
    }
}


 // ACCESSORS


 // Pre and Post processing utilites

 template<class UMCMP>
 template<class FT, class FT1, class FT2>
 void  MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT,FT1,FT2>::mapAvgTemp(FT1 &matTemps, 
					    const FT &avgTemp) const
{
    cpCell2Mats(matTemps, avgTemp);
}

 template<class UMCMP>
 template<class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT,FT1,FT2>::cpCell2Mats(FT1 &matValues,
					     const FT &cellValue) const
{
    int ncell = volumeFraction.size();
    Require(ncell == cellValue.size());
    Require(ncell == matValues.size());
    FT::const_iterator cit = cellValue.begin();
    FT1::iterator mvit = matValues.begin();
    for (int icell = 0; icell < ncell; icell++, cit++, mvit++)
    {
	int nmat = volumeFraction[icell].size();
	Require ( nmat == (*mvit).size());
	for (FT1::value_type::iterator mvitit = (*mvit).begin(); 
	     mvitit != (*mvit).end(); mvitit++)
	{
	    *mvitit = *cit;
	}
    }
}

 template<class UMCMP>
 template<class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT,FT1,FT2>::
 getElectronEnergyDepbyMat(FT1 &results, 
			   const FT &newElectronTemp,
			   const FT &cellVolume ) const
{
    int ncell = volumeFraction.size();
    Require(ncell == newElectronTemp.size());
    Require(ncell == results.size());
    Require(ncell == cellVolume.size());

    FT1 temp = results;
    FT1 cv   = results;
    getElectronTempByMat(temp);
    getElectronSpecificHeatByMat(cv);

    FT1::iterator tempit = temp.begin();
    FT1::iterator cvit   = cv.begin();
    FT1::iterator reit   = results.begin();
    FT::const_iterator volit   = cellVolume.begin();
    FT::const_iterator ntit = newElectronTemp.begin();
    for (int icell = 0; icell < ncell; icell++, ntit++, volit++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*reit).size() == nmat );
	FT1::value_type::iterator tempitit = (*tempit++).begin();
	FT1::value_type::iterator cvitit   = (*cvit++).begin();
	FT1::value_type::iterator reitit   = (*reit++).begin();
	for (int imat=0; imat < nmat; imat++)
	{
	    *reitit++ = volumeFraction[icell][imat] * (*cvitit++) * 
		(*volit) * ( (*ntit) - (*tempitit++) ) ;
	}
    }
}


 template<class UMCMP>
 template<class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT,FT1,FT2>::
 getIonEnergyDepbyMat(FT1 &results, 
		      const FT &newIonTemp, 
		      const FT &cellVolume) const
{
    int ncell = volumeFraction.size();
    Require(ncell == newIonTemp.size());
    Require(ncell == results.size());
    Require(ncell == cellVolume.size());

    FT1 temp = results;
    FT1 cv   = results;
    getIonTempByMat(temp);
    getIonSpecificHeatByMat(cv);

    FT1::iterator tempit = temp.begin();
    FT1::iterator cvit   = cv.begin();
    FT1::iterator reit   = results.begin();
    FT::const_iterator volit   = cellVolume.begin();
    FT::const_iterator ntit = newIonTemp.begin();
    for (int icell = 0; icell < ncell; icell++, ntit++, volit++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*reit).size() == nmat );
	FT1::value_type::iterator tempitit = (*tempit++).begin();
	FT1::value_type::iterator cvitit   = (*cvit++).begin();
	FT1::value_type::iterator reitit   = (*reit++).begin();
	for (int imat=0; imat < nmat; imat++)
	{
	    *reitit++ = volumeFraction[icell][imat] * (*cvitit++) * 
		(*volit) * ( (*ntit) - (*tempitit++) ) ;
	}
    }
}

 // Cell average properties
    
 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT,FT1,FT2>::getSigmaAbsorption(int group, FT &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	*resit = 0.;
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> sigtmp(nmat);
	matState[icell].getSigmaAbsorption(group, sigtmp);
	for (int imat = 0; imat < nmat; imat++)
	{
	    *resit += volumeFraction[icell][imat]*sigtmp[imat];
	}
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT,FT1,FT2>::getSigmaTotal(int group, FT &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	*resit = 0.;
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> sigtmp(nmat);
	matState[icell].getSigmaTotal(group, sigtmp);
	for (int imat = 0; imat < nmat; imat++)
	{
	    *resit += volumeFraction[icell][imat]*sigtmp[imat];
	}
    }
}

 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT,FT1,FT2>::getSigmaEmission(int group, FT &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	*resit = 0.;
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> sigtmp(nmat);
	matState[icell].getSigmaEmission(group, sigtmp);
	for (int imat = 0; imat < nmat; imat++)
	{
	    *resit += volumeFraction[icell][imat]*sigtmp[imat];
	}
    }
}



 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getElectronIonCoupling(FT &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	*resit = 0.;
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> proptmp(nmat);
	matState[icell].getElectronIonCoupling(proptmp);
	for (int imat = 0; imat < nmat; imat++)
	{
	    *resit += volumeFraction[icell][imat]*proptmp[imat];
	}
    }
}

 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getElectronConductionCoeff(FT &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	*resit = 0.;
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> proptmp(nmat);
	matState[icell].getElectronConductionCoeff(proptmp);
	for (int imat = 0; imat < nmat; imat++)
	{
	    *resit += volumeFraction[icell][imat]*proptmp[imat];
	}
    }
}

 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getIonConductionCoeff(FT &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	*resit = 0.;
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> proptmp(nmat);
	matState[icell].getIonConductionCoeff(proptmp);
	for (int imat = 0; imat < nmat; imat++)
	{
	    *resit += volumeFraction[icell][imat]*proptmp[imat];
	}
    }
}

 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getElectronSpecificHeat(FT &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	*resit = 0.;
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> proptmp(nmat);
	matState[icell].getElectronSpecificHeat(proptmp);
	for (int imat = 0; imat < nmat; imat++)
	{
	    *resit += volumeFraction[icell][imat]*proptmp[imat];
	}
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getIonSpecificHeat(FT &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	*resit = 0.;
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> proptmp(nmat);
	matState[icell].getIonSpecificHeat(proptmp);
	for (int imat = 0; imat < nmat; imat++)
	{
	    *resit += volumeFraction[icell][imat]*proptmp[imat];
	}
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getElectronTemperature(FT &results) const
{
    typedef typename FT1::value_type::value_type value_type1;
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> cvtmp(nmat);
	matState[icell].getElectronSpecificHeat(cvtmp);

	std::vector<value_type1> tetmp(nmat);
	matState[icell].getElectronTemperature(tetmp);
	
	value_type1 xsum = 0.;
	value_type1 xnrm = 0.;
	for (int imat = 0; imat < nmat; imat++)
	{
            value_type1 xx = cvtmp[imat]*volumeFraction[icell][imat];
	    xsum += xx*tetmp[imat];
	    xnrm += xx;
	}
	Require(xsum > 0.);
	Require(xnrm > 0.);
	*resit = xsum/xnrm;
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getIonTemperature(FT &results) const
{
    typedef typename FT1::value_type::value_type value_type1;
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	std::vector<value_type1> cvtmp(nmat);
	matState[icell].getIonSpecificHeat(cvtmp);

	std::vector<value_type1> titmp(nmat);
	matState[icell].getIonTemperature(titmp);
	
	value_type1 xsum = 0.;
	value_type1 xnrm = 0.;
	for (int imat = 0; imat < nmat; imat++)
	{
            value_type1 xx = cvtmp[imat]*volumeFraction[icell][imat];
	    xsum += xx*titmp[imat];
	    xnrm += xx;
	}
	Require(xsum > 0.);
	Require(xnrm > 0.);
	*resit = xsum/xnrm;
    }
}


 // Properties by material


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getElectronSpecificHeatByMat(FT1 &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT1::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*resit).size() == nmat );
	std::vector<value_type1> cvtmp(nmat);
	matState[icell].getElectronSpecificHeat(cvtmp);
	std::copy( cvtmp.begin(), cvtmp.end(), (*resit).begin() );
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getIonSpecificHeatByMat(FT1 &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT1::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*resit).size() == nmat );
	std::vector<value_type1> cvtmp(nmat);
	matState[icell].getIonSpecificHeat(cvtmp);
	std::copy( cvtmp.begin(), cvtmp.end(), (*resit).begin() );
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getElectronTempByMat(FT1 &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT1::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*resit).size() == nmat );
	std::vector<value_type1> tetmp(nmat);
	matState[icell].getElectronTemperature(tetmp);
	std::copy( tetmp.begin(), tetmp.end(), (*resit).begin() );
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getIonTempByMat(FT1 &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT1::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*resit).size() == nmat );
	std::vector<value_type1> titmp(nmat);
	matState[icell].getIonTemperature(titmp);
	std::copy( titmp.begin(), titmp.end(), (*resit).begin() );
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getDensity(FT1 &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT1::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*resit).size() == nmat );
	std::vector<value_type1> tmp(nmat);
	matState[icell].getDensity(tmp);
	std::copy( tmp.begin(), tmp.end(), (*resit).begin() );
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getVolumeFraction(FT1 &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT1::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*resit).size() == nmat );
	std::copy( volumeFraction[icell].begin(),
		   volumeFraction[icell].end(), (*resit).begin() );
    }
}


 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 void MultiMatCellMatProps<UMCMP>::
 MaterialStateField<FT, FT1, FT2>::getMatId(FT2 &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT2::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*resit).size() == nmat );
	std::vector<FT2::value_type::value_type> tmp(nmat);
	matState[icell].getMatId(tmp);
	std::copy( tmp.begin(), tmp.end(), (*resit).begin() );
    }
}


}

//---------------------------------------------------------------------------//
//                              end of MultiMatCellMatProps.t.cc
//---------------------------------------------------------------------------//
