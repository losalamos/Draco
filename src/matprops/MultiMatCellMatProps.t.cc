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

// CREATORS

template <class UMCMP>
template <class FT, class FT1>
template <class FT2>
MultiMatCellMatProps<UMCMP>::MaterialStateField<FT, FT1>::
MaterialStateField( const MultiMatCellMatProps<UMCMP> &matprops_,
		    const FT1 &density_, const FT1 &electronTemp_,
		    const FT1 &ionTemp_, const FT1 &volumeFraction_,
		    const FT2 &matId_)
{
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


// Cell average properties
    
template <class UMCMP>
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getSigmaAbsorption(int group, FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getSigmaTotal(int group, FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getSigmaEmission(int group, FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getElectronIonCoupling(FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getElectronConductionCoeff(FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getIonConductionCoeff(FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getElectronSpecificHeat(FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getIonSpecificHeat(FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getElectronTemperature(FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getIonTemperature(FT &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getElectronSpecificHeatByMat(FT1 &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getIonSpecificHeatByMat(FT1 &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getElectronTempByMat(FT1 &results) const
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
template <class FT, class FT1>
void MultiMatCellMatProps<UMCMP>::
MaterialStateField<FT, FT1>::getIonTempByMat(FT1 &results) const
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


//---------------------------------------------------------------------------//
//                              end of MultiMatCellMatProps.t.cc
//---------------------------------------------------------------------------//
