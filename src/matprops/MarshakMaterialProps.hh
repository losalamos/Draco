//----------------------------------*-C++-*----------------------------------//
// MarshakMaterialProps.hh
// Randy M. Roberts
// Wed Sep  9 10:09:41 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_MarshakMaterialProps_hh__
#define __matprops_MarshakMaterialProps_hh__

#include "units/PhysicalConstants.hh"
#include "units/Units.hh"
#include <vector>
#include <cmath>
#include <limits>

namespace XTM
{
 
 //===========================================================================//
 // class MarshakMaterialProps - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 class MarshakMaterialProps
 {

     // NESTED CLASSES AND TYPEDEFS

   public:

     // Forward declaration
    
     template<class FT> class MaterialStateField;

     // FRIENDS

     template<class FT>
     friend class MaterialStateField;

     // DATA

   private:
    
     Units units;

     double kappa0;
     double abar;
     int kappaPower;
     double gamma;
     
   public:

     // CREATORS
    
     MarshakMaterialProps(const Units &units_,
			  double kappa0_=10.0, double abar_=1.0,
			  int kappaPower_=3, double gamma_ = 5.0/3.0)
	 : units(units_), kappa0(kappa0_), abar(abar_),
	   kappaPower(kappaPower_), gamma(gamma_)
     {
	 // empty
     }
     
     // MarshakMaterialProps(const MarshakMaterialProps &rhs);
     // ~MarshakMaterialProps();

     // MANIPULATORS
    
     // MarshakMaterialProps& operator=(const MarshakMaterialProps &rhs);

     // ACCESSORS

  public:

    const Units &getUnits() const { return units; }

    //------------------------------------------------------------------------//
    // getMaterialState:
    //    Return a material state field from density, temperature,
    //    and materia id fields.
    //------------------------------------------------------------------------//

    template<class FT, class FT2>
    MaterialStateField<FT> getMaterialState(const FT &density_,
					    const FT &electronTemp_,
					    const FT &ionTemp_,
					    const FT2 &matId_) const
     {
	 return MaterialStateField<FT>(*this, density_, electronTemp_,
				       ionTemp_, matId_);
     }
     
   private:
    
     // IMPLEMENTATION
 };

//===========================================================================//
// class MarshakMaterialProps::MaterialStateField
//===========================================================================//

template<class FT>
class MarshakMaterialProps::MaterialStateField
{
    // FRIENDS
    
    friend class MarshakMaterialProps;

    // Nested Classes and Typedefs
    
  private:

    typedef typename FT::value_type value_type;

    // DATA

  private:

    int theSize;
    
    const MarshakMaterialProps  *pMatprops;

    std::vector<value_type>      density;
    std::vector<value_type>      electronTemp;
    std::vector<value_type>      ionTemp;
    std::vector<int>             matId;
	

    // CREATORS

  private:

    template<class FT2>
    MaterialStateField(const MarshakMaterialProps &matprops_,
		       const FT &density_, const FT &electronTemp_,
		       const FT &ionTemp_, const FT2 &matId_)
	: pMatprops(&matprops_), theSize(density_.size()),
	  density(density_.begin(), density_.end()),
	  electronTemp(electronTemp_.begin(), electronTemp_.end()),
	  ionTemp(ionTemp_.begin(), ionTemp_.end()),
	  matId(matId_.begin(), matId_.end())
    {
	Require(density_.size() == electronTemp_.size());
	Require(density_.size() == ionTemp_.size());
	Require(density_.size() == matId_.size());
    }

    // ACCESSORS
	
  public:

    int size() const { return theSize; }

    const MarshakMaterialProps &getProps() const { return *pMatprops; }
    
    const Units &getUnits() const { return getProps().getUnits(); }

    inline void getElectronTemperature(FT &results) const;
    inline void getIonTemperature(FT &results) const;
    inline void getDensity(FT &results) const;

    inline void getSigmaAbsorption(int group, FT &results) const;

    void getSigmaAbsorption(double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    void getSigmaTotal(int group, FT &results) const
    {
	getSigmaAbsorption(group, results);
    }

    void getSigmaTotal(double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    void getSigmaEmission(int group, FT &results) const
    {
	getSigmaAbsorption(group, results);
    }

    void getSigmaEmission(double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    void getElectronIonCoupling(FT &results) const
    {
	const double gammaSI = 1.0e+20; // seconds^-1
	results = getUnits().ConvertTime(gammaSI);
    }

    void getElectronConductionCoeff(FT &results) const
    {
	results = 0;
    }

    void getIonConductionCoeff(FT &results) const
    {
	results = 0;
    }

    inline void getElectronSpecificHeat(FT &results) const;

    void getIonSpecificHeat(FT &results) const
    {
	getElectronSpecificHeat(results);
    }

    std::ostream &print(std::ostream &os)
    {
	for (int i=0; i<memento.size(); i++)
	    os << memento[i];
	return os;
    }
    
    // IMPLEMENTATION

  private:
    
    const value_type &getDensity(int i) const { return density[i]; }
    const value_type &getElectronTemp(int i) const { return electronTemp[i]; }
    const value_type &getIonTemp(int i) const { return ionTemp[i]; }
    const int &getMatId(int i) const { return matId[i]; }

};

#define MPMSF MarshakMaterialProps::MaterialStateField

template<class FT>
void MPMSF<FT>::getElectronTemperature(FT &results) const
{
    Require(size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < size(); i++)
	*resit++ = getElectronTemp(i);
}

template<class FT>
void MPMSF<FT>::getIonTemperature(FT &results) const
{
    Require(size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < size(); i++)
	*resit++ = getIonTemp(i);
}

template<class FT>
void MPMSF<FT>::getDensity(FT &results) const
{
    Require(size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < size(); i++)
	*resit++ = getDensity(i);
}

template<class FT>
void MPMSF<FT>::getSigmaAbsorption(int group, FT &results) const
{
    Require(size() == results.size());
	
    double kappa0 = pMatprops->kappa0;
    int kappaPower = pMatprops->kappaPower;

    using std::pow;
    
    FT::iterator resit = results.begin();

    for (int i=0; i < size(); i++)
	*resit++ = getDensity(i)*kappa0/pow(getElectronTemp(i), kappaPower);
}
    
template<class FT>
void MPMSF<FT>::getElectronSpecificHeat(FT &results) const
{
    Require(size() == results.size());
    
    const double abar = pMatprops->abar;
    const double protonMass =
	getUnits().InvertMass(PhysicalConstants::protonMassSI);
    
    const double ionMass = abar * protonMass;
    
    const double gamma = pMatprops->gamma;
    const double kSI = PhysicalConstants::boltzmannSI;

    const double k = getUnits().InvertEnergy(
	getUnits().ConvertTemperature(kSI));

    double Cv = k/(gamma-1.0);

    // Convert from a 3T Cv to a 2T Cv

    Cv /= 2.0;
    
    FT::iterator resit = results.begin();
    
    for (int i=0; i < size(); i++)
	*resit++ = Cv * getDensity(i) / ionMass;

}

#undef MPMSF

} // namespace XTM

#endif                          // __matprops_MarshakMaterialProps_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/MarshakMaterialProps.hh
//---------------------------------------------------------------------------//
