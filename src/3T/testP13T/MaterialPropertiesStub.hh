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
    typedef int MaterialStateField;

    // DATA

  private:

    Units units;
    SP<MT> spMesh;
    
  public:

    // CREATORS
    
    MaterialPropertiesStub(const Units &units_, const SP<MT> &spMesh_);

    // MANIPULATORS
    
    // ACCESSORS

    void print(std::ostream &os) const
    {
	os << "in MaterialPropertiesStub::print(), this: "
	   << (void *)(this)
	   << " spMesh: ";
	spMesh->print(os);
    }
    
    const SP<MT> getMesh() const { return spMesh; }
    Units getUnits() const { return units; }

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

END_NS_XTM  // namespace XTM

#endif                          // __3T_testP13T_MaterialPropertiesStub_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/MaterialPropertiesStub.hh
//---------------------------------------------------------------------------//
