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
#include "ds++/SP.hh"
#include "units/Units.hh" 

namespace rtt_matprops {

 //===========================================================================//
 // class MultiMatCellMatProps - Material properties for multi-material cells.

 // This class provides a means to deal with computational cells which 
 // contain more than one material. This is most often seen in Eulerian codes.
 // The method used is a simplistic one modeled on a
 // "Temperature-Equilibration" principle. This class is templated on a
 // Uni-Material-Cell Mat-Props (UMCMP) class. Each cell has a 
 // UMCMP::materialStateField that is nmats long, where nmats is the number
 // of materials in the cell.
 //
 // Usage is simple. Just construct a smart pointer to your favorite
 // uni-material-cell material properties object (spumcmp), and then 
 // use that to constuct the multi-material-cell material properties
 // object (mmcmp). Next construct a MMCMP::MaterialStateField (msf) using
 // density, temperature, volume fraction, and material id's. You can
 // then get cell average or material specific data such as opacities, 
 // specific heats, etc., etc. from the methods of the msf.
 // 
 //===========================================================================//

 template <class UMCMP>  // Uni-Material Cell Material Properties
 class MultiMatCellMatProps {

     // NESTED  CLASSES AND TYPEDEFS

   public:

     template<class FT, class FT1, class FT2>
     class MaterialStateField;

     // FRIENDS
    
   private:

     template <class U>
     template <class FT, class FT1, class FT2>
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


     // returns the units to be used for mass, length, time, 
     // and temperature
     const XTM::Units &getUnits() const
     {
	 return spumcmp->getUnits();
     }

     // returns the number of energy groups in the problem
     int getNumGroups(int materialId) const
     {
	 return spumcmp->getNumGroups(materialId);
     }

     // returns the maximum Pn scattering order to be found in
     // all the materials in the problem.
     int getMaxScatteringPnOrder(int materialId) const
     {
	 return spumcmp->getMaxScatteringPnOrder(materialId);
     }

     //------------------------------------------------------------------------//
     // getMaterialState:
     //    Returns a material state field from density, temperature,
     //    volume, and material id fields. The inputs arguments are
     //    expected to be "containers of containers". e.g. density_ is
     //    a container of size ncells, with each element itself being a
     //    container of size nmat where  nmat may be different for each cell.
     //    FT2 is provided in the anticipation that the material Id's will
     //    be integer based, whereas the other arguments (FT1) are 
     //    anticipated to be based on double.
     //
     //    density_  - material density (mass/length**3). Should be based on 
     //    the volume of the individual material within the cell, not 
     //    the total cell volume.
     //
     //    electronTemp_ - Electron temperatures (temperature).
     //
     //    ionTemp_ - Ion temperatures (temperature).
     //
     //    volumeFraction - Volume (length**3) or volume fraction (dimensionless)
     //    of each material in the cell. Sum of volumeFraction within a cell 
     //    will be normalized to one.
     //
     //    matId_ - Material Id flags.
     //
     //    FT is the class for output by the various methods which return
     //    a single quantity on a cell, e.g. a cell average result.
     //
     //------------------------------------------------------------------------//

     template<class FT, class FT1, class FT2>
     MaterialStateField<FT, FT1, FT2> 
     getMaterialState(const FT1 &density_,
		      const FT1 &electronTemp_,
		      const FT1 &ionTemp_,
		      const FT1 &volumeFraction_,
		      const FT2 &matId_) const
     {
	 return MaterialStateField<FT,FT1,FT2>(*this, density_, electronTemp_,
					       ionTemp_, volumeFraction_, 
					       matId_);
     }
 };


 //===========================================================================//
 // class MaterialStateField
 //===========================================================================//

 template <class UMCMP>
 template <class FT, class FT1, class FT2>
 class MultiMatCellMatProps<UMCMP>::MaterialStateField
 {

     // NESTED TYPEDEFS

     typedef typename FT::value_type value_type;
     typedef typename FT1::value_type::value_type value_type1;
     typedef std::vector<value_type1> VV1;
     typedef typename FT2::value_type::value_type value_type2;
     typedef std::vector<value_type2> VV2;
     typedef typename UMCMP:: template MaterialStateField<VV1,VV2> UmcMsf;

     // FRIENDS
    
     friend class MultiMatCellMatProps;

     // DATA
    
   private:

     const MultiMatCellMatProps<UMCMP> *pMatProps;

     //  std::vector<UMCMP::MaterialStateField<VV1> >  matState;
     std::vector<UmcMsf>  matState;
    
     std::vector<std::vector<value_type1> > volumeFraction;
    
     // matState.size() = volumeFraction.size() = ncells
     // matState[i].size() = volumeFraction[i].size() = nmat in cell i

     // CREATORS

   private:

     MaterialStateField( const MultiMatCellMatProps<UMCMP> &matprops_,
			 const FT1 &density_, const FT1 &electronTemp_,
			 const FT1 &ionTemp_, const FT1 &volumeFraction_,
			 const FT2 &matId_);
   public:
    
     ~MaterialStateField()
     {
	 // empty
     }


     // ACCESSORS

   public:

     // returns the units to be used for mass, length, time, 
     // and temperature
     const XTM::Units &getUnits() const { return pMatProps->getUnits(); }

     // For all the methods found below, the size of
     // arguments of type FT1 and FT (i.e. ncell) is  expected to be the
     // same as the size of the type FT supplied to the constructor.
     // If the size of the individual containers within arguments of
     // type FT1 (i.e. nmats) is not the same as the size of the individual
     // containers found in the FT1 supplied to the constructor, the
     // individual containers will be resized to agree with the
     // size supplied to the constructor.


     // *** Pre and Post processing utilities ***


     // Duplicates a single value on a cell for each material in
     // the cell
     void cpCell2Mats(FT1 &results, const FT &cellValue) const;

     // Maps a cell average electron temperature to each of the materials
     // found in the cell.
     void mapAvgElectronTemp(FT1 &results, const FT &avgElectronTemp) const;

     // Maps a cell average ion temperature to each of the materials
     // found in the cell.
     void mapAvgIonTemp(FT1 &results, const FT &avgIonTemp) const;

     // Given a new cell average electron temperature (temperature),
     // and a cell volume (length**3), returns the energy depositied in 
     // the electrons of each material in the cell (energy).
     void getElectronEnergyDepbyMat(FT1 &results, 
				    const FT &newElectonTemp,
				    const FT &cellVolume) const;

     // Given a new cell average ion temperature (temperature), and a cell
     // volume (lenght**3), returns the energy depositied in the
     // ions of each material in the cell (energy).
     void getIonEnergyDepbyMat(FT1 &results, 
			       const FT &newIonTemp,
			       const FT &cellVolume) const;


     // *** Cell Average properties. ***
    

     // returns absorption opacity (1/length)
     void getSigmaAbsorption(const int group, FT &results) const;
     
     // returns scattering opacity (1/length)
     void getSigmaScattering(const int group, FT &results) const;
     
     // returns total opacity (1/length)
     void getSigmaTotal(const int group, FT &results) const;

     // returns emission opacity (1/length)
     void getSigmaEmission(const int group, FT &results) const;

     // returns electron-ion coupling coefficient 
     // (energy/length**3-time-temperature)
     void getElectronIonCoupling(FT &results) const;

     // returns electron thermal conduction coefficient 
     // (energy/length-time-temperature)
     void getElectronConductionCoeff(FT &results) const;

     // returns ion thermal conduction coefficient 
     // (energy/length-time-temperature)
     void getIonConductionCoeff(FT &results) const;

     // returns electron specific heat
     // (energy/length**3-temperature)
     void getElectronSpecificHeat(FT &results) const;

     // returns ion specific heat
     // (energy/length**3-temperature)
     void getIonSpecificHeat(FT &results) const;

     // returns electron temperature (temperature)
     void getElectronTemperature(FT &results) const;

     // returns ion temperature (temperature)
     void getIonTemperature(FT &results) const;
    

     // *** Properties by material ***


     // returns electron specific heat
     // (energy/length**3-temperature)
     void getElectronSpecificHeatByMat(FT1 &results) const;

     // returns ion specific heat
     // (energy/length**3-temperature)
     void getIonSpecificHeatByMat(FT1 &results) const;

     // returns electron temperature (temperature)
     void getElectronTempByMat(FT1 &results) const;

     // returns ion temperature (temperature)
     void getIonTempByMat(FT1 &results) const;

     // returns density supplied at construction (mass/length**3)
     // density is based on the individual material volumes,
     // not the total cell volume.
     void getDensity(FT1 &results) const;

     // Returns material volume fraction (dimensionless) in each cell.
     // The sum of the volume fraction equals one in each cell.
     void getVolumeFraction(FT1 &results) const;

     // returns material Id flags supplied at construction
     void getMatId(FT2 &results) const;

 };

} // end of rtt_matprops namespace

#endif                          // __matprops_MultiMatCellMatProps_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/MultiMatCellMatProps.hh
//---------------------------------------------------------------------------//
