//----------------------------------*-C++-*----------------------------------//
// InterpedMaterialProps.hh
// Randy M. Roberts
// Wed Apr 15 08:44:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_InterpedMaterialProps_hh__
#define __matprops_InterpedMaterialProps_hh__

#include "matprops/BilinearInterpGrid.hh"
#include "matprops/BilinearInterpTable.hh"
#include "ds++/SP.hh"
#include "3T/Units.hh"
#include <vector>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class InterpedMaterialProps - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class InterpedMaterialProps
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    template<class FT>
    class MaterialStateField;

  private:

    struct GroupedTable
    {
	double              energy;
	BilinearInterpTable table;
	bool operator<(const GroupedTable &rhs) { return energy < rhs.energy; }
	bool operator<(double rhs) { return energy < rhs; }
    };
	
    // DATA

  private:
    
    Units units;

    BilinearInterpGrid               grid;
    std::vector<GroupedTable>        sigmaTotal;
    std::vector<GroupedTable>        sigmaAbsorption;
    std::vector<GroupedTable>        sigmaEmission;
    BilinearInterpTable              electronIonCoupling;
    BilinearInterpTable              electronConductionCoeff;
    BilinearInterpTable              ionConductionCoeff;
    BilinearInterpTable              electronSpecificHeat;
    BilinearInterpTable              ionSpecificHeat;

    
  public:

    // CREATORS
    
    InterpedMaterialProps(MaterialPropsFactory &factory);

    // MANIPULATORS
    
    // *** none ***

    // ACCESSORS

    Units getUnits() const { return units; }

    template<class FT>
    MaterialStateField<FT> getMaterialState(const FT &density_,
					    const FT &electronTemp_,
					    const FT &ionTemp_) const;    

    //------------------------------------------------------------------------//
    // getSigmaTotal:
    //------------------------------------------------------------------------//

    template<class FT>
    void getSigmaTotal(const MaterialStateField<FT> &matState,
		       int group, FT &results) const;

    template<class FT>
    void getSigmaTotal(const MaterialStateField<FT> &matState,
		       double group, FT &results) const;

    template<class FT>
    void getSigmaAbsorption(const MaterialStateField<FT> &matState,
			    int group, FT &results) const;

    template<class FT>
    void getSigmaAbsorption(const MaterialStateField<FT> &matState,
			    double group, FT &results) const;

    template<class FT>
    void getSigmaEmission(const MaterialStateField<FT> &matState,
			  int group, FT &results) const;

    template<class FT>
    void getSigmaEmission(const MaterialStateField<FT> &matState,
			  double group, FT &results) const;

    template<class FT>
    void getElectronIonCoupling(const MaterialStateField<FT> &matState,
				FT &results) const;

    template<class FT>
    void getElectronTemperature(const MaterialStateField<FT> &matState,
				FT &results) const;

    template<class FT>
    void getIonTemperature(const MaterialStateField<FT> &matState,
			   FT &results) const;

    template<class FT>
    void getElectronConductionCoeff(const MaterialStateField<FT> &matState,
				    FT &results) const;

    template<class FT>
    void getIonConductionCoeff(const MaterialStateField<FT> &matState,
			       FT &results) const;

    template<class FT>
    void getElectronSpecificHeat(const MaterialStateField<FT> &matState,
				 FT &results) const;

    template<class FT>
    void getIonSpecificHeat(const MaterialStateField<FT> &matState,
			    FT &results) const;
    
  private:
    
    // IMPLEMENTATION
};

//===========================================================================//
// class MaterialStateField
//===========================================================================//

template<class FT>
class InterpedMaterialProps::MaterialStateField
{
    friend class InterpedMaterialProps;
    
  private:
    
    const InterpedMaterialProps &matprops;

    std::vector<BilinearInterpTable::Memento> memento;
    FT                                        density;
    FT                                        electronTemp;
    FT                                        ionTemp;

  private:

    MaterialStateField(const InterpedMaterialProps &matprops_,
		       const FT &density_, const FT &electronTemp_,
		       const FT &ionTemp_);
    
};

END_NS_XTM  // namespace XTM

#endif                          // __matprops_InterpedMaterialProps_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/InterpedMaterialProps.hh
//---------------------------------------------------------------------------//
