//----------------------------------*-C++-*----------------------------------//
// MaterialPropertiesStub.hh
// Randy M. Roberts
// Mon Mar 23 08:41:02 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_MaterialPropertiesStub_hh__
#define __3T_testP13T_MaterialPropertiesStub_hh__

#include "ds++/SP.hh"
#include "3T/Units.hh"
#include <iostream>
#include <vector>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class MaterialPropertiesStub - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class MaterialPropertiesStub
{

    // NESTED CLASSES AND TYPEDEFS

  public:
    
    struct MatVals
    {
	double sigmaTotal;
	double sigmaAbsorption;
	double sigmaEmission;
	double electronIonCoupling;
	double electronConductionCoeff;
	double ionConductionCoeff;
	double electronSpecificHeat;
	double ionSpecificHeat;

	MatVals(double sT=0.0, double sA=0.0, double sE=0.0,
		double eic=0.0, double ecc=0.0, double icc=0.0,
		double esh=0.0, double ish=0.0)
	    : sigmaTotal(sT), sigmaAbsorption(sA), sigmaEmission(sE),
	      electronIonCoupling(eic), electronConductionCoeff(ecc),
	      ionConductionCoeff(icc), electronSpecificHeat(esh),
	      ionSpecificHeat(ish)
	{
	    // empty
	}
	      
    };

    typedef int MaterialStateField;

    // DATA

    std::vector<MatVals> matVals;

  private:

    Units units;
    SP<MT> spMesh;

    double TElectron;
    double TIon;
    
  public:

    // CREATORS
    
    MaterialPropertiesStub(const Units &units_, const SP<MT> &spMesh_,
			   double TElectron_, double TIon_,
			   const std::vector<MatVals> &matVals_);

    // MANIPULATORS
    
    // ACCESSORS

    std::ostream &print(std::ostream &os) const
    {
	os << "(MaterialPropertiesStub::this: "
	   << (void *)(this)
	   << " spMesh: " << *spMesh << ")";
	return os;
    }
    
    const SP<MT> getMesh() const { return spMesh; }
    Units getUnits() const { return units; }

    //------------------------------------------------------------------------//
    // getSigmaTotal:
    //     It is the material properties' responsibility to do
    //     any averaging of temperatures, etc. to achieve the correct
    //     resulting sigmaTotal.
    //------------------------------------------------------------------------//

    template<class FT, class GT>
    void getSigmaTotal(const MaterialStateField &matState,
		       GT group, FT &sigmaTotal) const;

    template<class FT, class GT>
    void getSigmaAbsorption(const MaterialStateField &matState,
			    GT group, FT &sigmaAbsorption) const;

    template<class FT, class GT>
    void getSigmaEmission(const MaterialStateField &matState,
			  GT group, FT &sigmaEmission) const;

    template<class FT>
    void getElectronIonCoupling(const MaterialStateField &matState,
				FT &gamma) const;

    template<class FT>
    void getElectronTemperature(const MaterialStateField &matState,
				FT &TElectron) const;

    template<class FT>
    void getIonTemperature(const MaterialStateField &matState,
			   FT &TIon) const;

    template<class FT>
    void getElectronConductionCoeff(const MaterialStateField &matState,
				    FT &kappaElectron) const;

    template<class FT>
    void getIonConductionCoeff(const MaterialStateField &matState,
			       FT &kappaIon) const;

    template<class FT>
    void getElectronSpecificHeat(const MaterialStateField &matState,
				 FT &CvElectron) const;

    template<class FT>
    void getIonSpecificHeat(const MaterialStateField &matState,
			    FT &CvIon) const;

  private:
    
    // IMPLEMENTATION
};

template<class MT>
inline std::ostream &operator<<(std::ostream &os,
				const MaterialPropertiesStub<MT> &rhs)
{
    return rhs.print(os);
}

END_NS_XTM  // namespace XTM

#endif                          // __3T_testP13T_MaterialPropertiesStub_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/MaterialPropertiesStub.hh
//---------------------------------------------------------------------------//
