//----------------------------------*-C++-*----------------------------------//
// InterpedMaterialProps.hh
// Randy M. Roberts
// Wed Apr 15 08:44:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_InterpedMaterialProps_hh__
#define __matprops_InterpedMaterialProps_hh__

#include "BilinearInterpGrid.hh"
#include "BilinearInterpTable.hh"
#include "ds++/SP.hh"
#include "units/Units.hh"

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

namespace rtt_matprops {

// Forward Reference

class MaterialPropsReader;


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

    // Forward declaration
    
    template<class FT, class FT2> class MaterialStateField;

    // FRIENDS

    template<class FT, class FT2>
    friend class MaterialStateField;

  private:

    // Forward declaration
    
    struct GroupedTable;

    // Forward declaration
    
    struct MaterialTables;

    // Definine the type that contains a dictionary of material tables
    // keyed via the integer material id number.
    
    typedef std::map<int, MaterialTables> MatTabMap;


    // DATA

  private:
    
    rtt_units::Units units;

    // This is the dictionary of material tables
    // keyed via the integer material id number.
    
    MatTabMap materials;
    

    // CREATORS
    
  public:

    InterpedMaterialProps(const std::vector<int> &materialIds,
			  MaterialPropsReader &reader);


    // MANIPULATORS
    
    // *** none ***


    // ACCESSORS

  public:

    const rtt_units::Units &getUnits() const { return units; }

    //------------------------------------------------------------------------//
    // getMaterialState:
    //    Return a material state field from density, temperature,
    //    and materia id fields.
    //------------------------------------------------------------------------//

    template<class FT, class FT2>
    MaterialStateField<FT,FT2> getMaterialState(const FT &density_,
						const FT &electronTemp_,
						const FT &ionTemp_,
						const FT2 &matId_) const;

    inline const std::string &getMaterialName(int materialId) const;

    inline int getNumGroups(int materialId) const;

    inline const std::vector<double> &getDensityGrid(int materialId) const;
    
    inline int getNumDensities(int materialId) const;

    inline const std::vector<double> &getTemperatureGrid(int materialId) const;
    
    inline int getNumTemperatures(int materialId) const;

    inline int getMaxScatteringPnOrder(int materialId) const;

    
    // IMPLEMENTATION

  private:
    
    //------------------------------------------------------------------------//
    // hasMaterialTables:
    //   Return whether the material tables specified by the material id exists.
    //------------------------------------------------------------------------//

    bool hasMaterialTables(int matId) const;
	
    //------------------------------------------------------------------------//
    // getMaterialTables:
    //   Return the material tables specified by the material id.
    //------------------------------------------------------------------------//

    const MaterialTables &getMaterialTables(int matId) const;
	
    //------------------------------------------------------------------------//
    // getMaterialTables:
    //   Return the material tables specified by the material id.
    //------------------------------------------------------------------------//

    MaterialTables &getMaterialTables(int matId);

    //------------------------------------------------------------------------//
    // interpolate:
    //   Given a material state, a group number, and a pointer to a 
    //   GroupedTable member, return a unary operation on the
    //   interpolated results from the table indicated by the method pointer
    //   and group.
    //------------------------------------------------------------------------//

    template<class FT, class FT2, class UnaryOperation>
    void interpolate(const MaterialStateField<FT,FT2> &matState, int group,
		     const GroupedTable MaterialTables::* pTable,
		     UnaryOperation op, FT &results) const;

    //------------------------------------------------------------------------//
    // interpolate:
    //   Given a material state, and a pointer to a pointer to a 
    //   BilinearInterpTable member, return a unary operation on the
    //   interpolated results from the table indicated by the method pointer.
    //------------------------------------------------------------------------//

    template<class FT, class FT2, class UnaryOperation>
    void interpolate(const MaterialStateField<FT,FT2> &matState,
		     const BilinearInterpTable MaterialTables::* pTable,
		     UnaryOperation op, FT &results) const;
};

//========================================================================//
// class InterpedMaterialProps::GroupedTable
//    This class is defined for material properties which vary
//    by group or continuous energy.
//    For energy, one must interpolate between groups (not yet implemented).
//========================================================================//

struct InterpedMaterialProps::GroupedTable
{
    // DATA

    int ngroups;
	
    // all of these vectors are ngroups long.
	
    std::vector<double>              energyUpperbounds;
    std::vector<double>              energyLowerbounds;
    std::vector<BilinearInterpTable> tables;

    // CREATORS

    GroupedTable() : ngroups(0) {};
    
    GroupedTable(const std::vector<double> &energyUpperbounds_,
		 const std::vector<double> &energyLowerbounds_,
		 const std::vector<BilinearInterpTable> tables_)
	: energyUpperbounds(energyUpperbounds_),
	  energyLowerbounds(energyLowerbounds_), tables(tables_)
    {
	ngroups = energyUpperbounds.size();
	Require(energyLowerbounds.size() == ngroups);
	Require(tables.size() == ngroups);
    }
    
    // ACCESSORS
	
    int numGroups() const { return ngroups; }

    const BilinearInterpTable &getTable(int groupNo) const
    {
	// Groups are one-based.
	return tables[groupNo-1];
    }

    // MANIPULATORS
	
    BilinearInterpTable &getTable(int groupNo)
    {
	// Groups are one-based.
	return tables[groupNo-1];
    }
};

//===========================================================================//
// class InterpedMaterialProps::MaterialTables
//===========================================================================//

struct InterpedMaterialProps::MaterialTables
{
    // DATA

    std::string              materialName;
    int                      numGroups;
    int                      numDensities;
    int                      numTemperatures;
    int                      maxScatteringPnOrder;
    
    rtt_dsxx::SP<BilinearInterpGrid>   spGrid;

    GroupedTable             sigmaTotal;
    GroupedTable             sigmaAbsorption;
    GroupedTable             sigmaScattering;
    GroupedTable             sigmaEmission;
    BilinearInterpTable      electronIonCoupling;
    BilinearInterpTable      electronConductionCoeff;
    BilinearInterpTable      ionConductionCoeff;
    BilinearInterpTable      electronSpecificHeat;
    BilinearInterpTable      ionSpecificHeat;

    // CREATORS

    MaterialTables() {};
    
    MaterialTables(const std::string &materialName_,
		   const rtt_dsxx::SP<BilinearInterpGrid> &spGrid_,
		   const GroupedTable &sigmaTotal_,
		   const GroupedTable &sigmaAbsorption_,
		   const GroupedTable &sigmaScattering_,
		   const GroupedTable &sigmaEmission_,
		   const BilinearInterpTable &electronIonCoupling_,
		   const BilinearInterpTable &electronConductionCoeff_,
		   const BilinearInterpTable &ionConductionCoeff_,
		   const BilinearInterpTable &electronSpecificHeat_,
		   const BilinearInterpTable &ionSpecificHeat_)
	: materialName(materialName_),
	  spGrid(spGrid_),
	  sigmaTotal(sigmaTotal_),
          sigmaAbsorption(sigmaAbsorption_),
          sigmaScattering(sigmaScattering_),
	  sigmaEmission(sigmaEmission_),
          electronIonCoupling(electronIonCoupling_),
          electronConductionCoeff(electronConductionCoeff_),
          ionConductionCoeff(ionConductionCoeff_),
          electronSpecificHeat(electronSpecificHeat_),
          ionSpecificHeat(ionSpecificHeat_)
    {
	numGroups = sigmaTotal.numGroups();
	numTemperatures = spGrid->size(1);
	numDensities = spGrid->size(2);
	maxScatteringPnOrder = 0;
    }

    // ACCESSORS
	
    const std::string &getMaterialName() const
    {
	return materialName;
    }
    
    const BilinearInterpGrid &getGrid() const
    {
	return *spGrid;
    }
};

//===========================================================================//
// class InterpedMaterialProps::MaterialStateField
//===========================================================================//

template<class FT, class FT2>
class InterpedMaterialProps::MaterialStateField
{
    // FRIENDS
    
    friend class InterpedMaterialProps;

    // Nested Classes and Typedefs
    
  private:

    typedef BilinearInterpTable::Memento Memento;

    typedef typename FT::value_type value_type;

    // A small unary operation class
    
    class MultByDensity
    {
	const MaterialStateField &field;
	int index;

      public:

	MultByDensity(const MaterialStateField &field_)
	    : field(field_), index(0)
	{ /* empty */ }
	
	value_type operator()(const value_type &rhs)
	{
	    return field.getDensity(index++) * rhs;
	}
    };

    friend class MultByDensity;
    
    // DATA

  private:

    int theSize;
    
    const InterpedMaterialProps *pMatprops;

    std::vector<Memento>         memento;
    
    std::vector<value_type>      density;
    std::vector<value_type>      electronTemp;
    std::vector<value_type>      ionTemp;
    std::vector<int>             matId;
	

    // CREATORS

  private:

    MaterialStateField(const InterpedMaterialProps &matprops_,
		       const FT &density_, const FT &electronTemp_,
		       const FT &ionTemp_, const FT2 &matId_);

    // ACCESSORS
	
  public:

    int size() const { return theSize; }

    const InterpedMaterialProps &getProps() const { return *pMatprops; }
    
    const rtt_units::Units &getUnits() const { return getProps().getUnits(); }

    void getElectronTemperature(FT &results) const;

    void getIonTemperature(FT &results) const;

    void getDensity(FT &results) const;

    void getMatId(FT2 &results) const;

    void getSigmaTotal(int group, FT &results) const
    {
	getProps().interpolate(*this, group, &MaterialTables::sigmaTotal,
			       MultByDensity(*this), results);
    }

    void getSigmaTotal(double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    void getSigmaAbsorption(int group, FT &results) const
    {
	getProps().interpolate(*this, group,
			       &MaterialTables::sigmaAbsorption,
			       MultByDensity(*this), results);
    }

    void getSigmaAbsorption(double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    void getSigmaScattering(int group, FT &results) const
    {
	getProps().interpolate(*this, group,
			       &MaterialTables::sigmaScattering,
			       MultByDensity(*this), results);
    }

    void getSigmaScattering(double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    void getSigmaEmission(int group, FT &results) const
    {
	getProps().interpolate(*this, group, &MaterialTables::sigmaEmission,
			       MultByDensity(*this), results);
    }

    void getSigmaEmission(double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    void getElectronIonCoupling(FT &results) const
    {
	getProps().interpolate(*this, &MaterialTables::electronIonCoupling,
			       MultByDensity(*this), results);
    }

    void getElectronConductionCoeff(FT &results) const
    {
	getProps().interpolate(*this,
			       &MaterialTables::electronConductionCoeff,
			       MultByDensity(*this), results);
    }

    void getIonConductionCoeff(FT &results) const
    {
	getProps().interpolate(*this, &MaterialTables::ionConductionCoeff,
			       MultByDensity(*this), results);
    }

    void getElectronSpecificHeat(FT &results) const
    {
	getProps().interpolate(*this, &MaterialTables::electronSpecificHeat,
			       MultByDensity(*this), results);
    }

    void getIonSpecificHeat(FT &results) const
    {
	getProps().interpolate(*this, &MaterialTables::ionSpecificHeat,
			       MultByDensity(*this), results);
    }

    std::ostream &print(std::ostream &os)
    {
	for (int i=0; i<memento.size(); i++)
	    os << memento[i];
	return os;
    }
    
    // IMPLEMENTATION

  private:
    
    const Memento &getMemento(int i) const { return memento[i]; }
    
    const value_type &getDensity(int i) const { return density[i]; }
    const value_type &getElectronTemp(int i) const { return electronTemp[i]; }
    const value_type &getIonTemp(int i) const { return ionTemp[i]; }
    const int &getMatId(int i) const { return matId[i]; }

};

// Inline definitions.

inline
const std::string &InterpedMaterialProps::getMaterialName(int materialId) const
{
    return getMaterialTables(materialId).materialName;
}

inline
int InterpedMaterialProps::getNumGroups(int materialId) const
{
    return getMaterialTables(materialId).numGroups;
}

inline
int InterpedMaterialProps::getNumDensities(int materialId) const
{
    return getMaterialTables(materialId).numDensities;
}

inline
int InterpedMaterialProps::getNumTemperatures(int materialId) const
{
    return getMaterialTables(materialId).numTemperatures;
}

inline
int InterpedMaterialProps::getMaxScatteringPnOrder(int materialId) const
{
    return getMaterialTables(materialId).maxScatteringPnOrder;
}
    
inline
const std::vector<double>
&InterpedMaterialProps::getDensityGrid(int materialId) const
{
    return getMaterialTables(materialId).spGrid->getGrid(2);
}
    
inline
const std::vector<double>
&InterpedMaterialProps::getTemperatureGrid(int materialId) const
{
    return getMaterialTables(materialId).spGrid->getGrid(1);
}

} // end of rtt_matprops namespace

#endif                          // __matprops_InterpedMaterialProps_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/InterpedMaterialProps.hh
//---------------------------------------------------------------------------//
