//----------------------------------*-C++-*----------------------------------//
// MultiMatCellMatProps.hh
// John McGhee
// Mon Sep 14 08:38:38 1998
//---------------------------------------------------------------------------//
// @> A material property class for cells which contain more than one material.
//---------------------------------------------------------------------------//

#ifndef __matprops_MultiMatCellMatProps_hh__
#define __matprops_MultiMatCellMatProps_hh__

#include <vector>
#include <algorithm>
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "units/Units.hh" 

//===========================================================================//
// class MultiMatCellMatProps - Material properties for multi-material cells.

// This class provides a means to deal with computational cells which 
// contain more than one material. This is most often seen in Eulerian codes.

// 
//===========================================================================//

template <class UMCMP>  // Uni-Material Cell Material Properties
class MultiMatCellMatProps {

// NESTED  CLASSES AND TYPEDEFS

  public:

    template<class FT, class FT1> class MaterialStateField;

// FRIENDS
    
  private:

    template <class U>
    template <class FT, class FT1>
    friend class MultiMatCellMatProps<U>::MaterialStateField;

// DATA

  private:
    
    dsxx::SP<UMCMP> spumcmp;

// CREATORS

  public:
    MultiMatCellMatProps(const dsxx::SP<UMCMP> &spumcmp_)
	:spumcmp(spumcmp_)
    {
    //empty
    }

    ~MultiMatCellMatProps()
    {
    //empty
    }

// ACCESSORS

  public:

    const XTM::Units &getUnits() const
    {
	return spumcmp->getUnits();
    }

    int getNumGroups(int materialId) const
    {
	return spumcmp->getNumGroups(materialId);
    }

    int getMaxScatteringPnOrder(int materialId) const
    {
	return spumcmp->getMaxScatteringPnOrder(materialId);
    }

//------------------------------------------------------------------------//
// getMaterialState:
//    Return a material state field from density, temperature,
//    and materia id fields.
//------------------------------------------------------------------------//

    template<class FT, class FT1, class FT2>
    MaterialStateField<FT, FT1> 
    getMaterialState(const FT1 &density_,
		     const FT1 &electronTemp_,
		     const FT1 &ionTemp_,
		     const FT1 &volumeFraction_,
		     const FT2 &matId_) const
    {
	return MaterialStateField<FT,FT1>(*this, density_, electronTemp_,
					  ionTemp_, volumeFraction_, matId_);
    }
};


//===========================================================================//
// class MaterialStateField
//===========================================================================//

template <class UMCMP>
template <class FT, class FT1>
class MultiMatCellMatProps<UMCMP>::MaterialStateField
{

// NESTED TYPEDEFS

    typedef typename FT::value_type value_type;
    typedef typename FT1::value_type::value_type value_type1;
    typedef std::vector<value_type1> VV1;
    typedef typename UMCMP:: template MaterialStateField<VV1> UmcMsf;
    

// FRIENDS
    
    friend class MultiMatCellMatProps;

// DATA
    
  private:

//  std::vector<UMCMP::MaterialStateField<VV1> >  matState;
    std::vector<UmcMsf>  matState;
    
    std::vector<std::vector<value_type1> > volumeFraction;
    
// matState.size() = volumeFraction.size() = ncells
// matState[i].size() = volumeFraction[i].size() = nmat in cell i

// CREATORS

  private:

    template<class FT2>
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
	// Load a material state field for this cell.
	    matState.push_back(matprops_.spumcmp->getMaterialState(
		*denit,*etit,*itit,*midit));
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
	
    
  public:
    
    ~MaterialStateField()
    {
    // empty
    }

// ACCESSORS

  public:

// Average properties
    
    inline void getSigmaAbsorption(int group, FT &results) const;
    inline void getSigmaTotal(int group, FT &results) const;
    inline void getSigmaEmission(int group, FT &results) const;
    inline void getElectronIonCoupling(FT &results) const;
    inline void getElectonConductionCoeff(FT &results) const;
    inline void getIonConductionCoeff(FT &results) const;
    inline void getElectronSpecificHeat(FT &results) const;
    inline void getIonSpecificHeat(FT &results) const;
    inline void getElectronTemp(FT &results) const;
    inline void getIonTemp(FT &results) const;
    
// Properties by material

    inline void getElectronSpecificHeatByMat(FT1 &results) const;
    inline void getIonSpecificHeatByMat(FT1 &results) const;
    inline void getElectronTempByMat(FT1 &results) const;
    inline void getIonTempByMat(FT1 &results) const;

};

// INLINE IMPLEMENTATIONS
    
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
MaterialStateField<FT, FT1>::getElectronSpecificHeatByMat(FT1 &results) const
{
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
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
MaterialStateField<FT, FT1>::getElectronTemp(FT &results) const
{
    typedef typename FT1::value_type::value_type value_type1;
    Require(matState.size() == results.size());
    int icell = 0;
    for (FT::iterator resit = results.begin(); 
	 resit != results.end(); resit++, icell++)
    {
	int nmat = volumeFraction[icell].size();
	Require( (*resit).size() == nmat );
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


#endif                          // __matprops_MultiMatCellMatProps_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/MultiMatCellMatProps.hh
//---------------------------------------------------------------------------//
