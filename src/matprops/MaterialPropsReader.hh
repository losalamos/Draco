//----------------------------------*-C++-*----------------------------------//
// MaterialPropsReader.hh
// Randy M. Roberts
// Mon Apr 20 10:27:42 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_MaterialPropsReader_hh__
#define __matprops_MaterialPropsReader_hh__

#include "ds++/Mat.hh"

#include <string>
#include <vector>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

//===========================================================================//
// class MaterialPropsReader - 
//    Abstract base class for reading material property data,
//    and giving that data to a Material Props object.
//===========================================================================//

class MaterialPropsReader
{

  public:

    virtual bool getNextMaterial(int &materialId_, std::string &name) const = 0;

    virtual void getTemperatureGrid(int materialId,
				    std::vector<double> &tempGrid_) const = 0;
    
    virtual void getDensityGrid(int materialId,
				std::vector<double> &densityGrid_) const = 0;

    virtual void getNumGroups(int materialId, int &numGroups) const = 0;


    virtual void getEnergyUpperbounds(int materialId, int group,
				      double &energyUpperbounds_) const = 0;

    virtual void getEnergyLowerbounds(int materialId, int group,
				      double &energyLowerbounds_) const = 0;

    virtual void getSigmaTotal(int materialId, int group,
			       dsxx::Mat2<double> &data) const = 0;

    virtual void getSigmaAbsorption(int materialId, int group,
				    dsxx::Mat2<double> &data) const = 0;

    virtual void getSigmaEmission(int materialId, int group,
				  dsxx::Mat2<double> &data) const = 0;

    virtual void getElectronIonCoupling(int materialId,
					dsxx::Mat2<double> &data) const = 0;
	
    virtual void getElectronConductionCoeff(int materialId,
					    dsxx::Mat2<double> &data) const = 0;
	
    virtual void getIonConductionCoeff(int materialId,
				       dsxx::Mat2<double> &data) const = 0;
	
    virtual void getElectronSpecificHeat(int materialId,
					 dsxx::Mat2<double> &data) const = 0;
	
    virtual void getIonSpecificHeat(int materialId,
				    dsxx::Mat2<double> &data) const = 0;

};

#endif                          // __matprops_MaterialPropsReader_hh__

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of matprops/MaterialPropsReader.hh
//---------------------------------------------------------------------------//
