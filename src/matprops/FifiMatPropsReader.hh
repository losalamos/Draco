//----------------------------------*-C++-*----------------------------------//
// FifiMatPropsReader.hh
// Randy M. Roberts
// Mon Apr 20 13:44:11 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_FifiMatPropsReader_hh__
#define __matprops_FifiMatPropsReader_hh__

#include "matprops/MaterialPropsReader.hh"

#include "ds++/Mat.hh"

#include <string>
#include <vector>
#include <iosfwd>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class FifiMatPropsReader - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class FifiMatPropsReader : public MaterialPropsReader
{
    // NESTED CLASSES AND TYPEDEFS

    // DATA

    // DISALLOWED DEFAULT METHODS
    
  private:
    
    FifiMatPropsReader(const FifiMatPropsReader &rhs);
    FifiMatPropsReader& operator=(const FifiMatPropsReader &rhs);

  public:

    // CREATORS
    
    FifiMatPropsReader(const Units &units_)
	: MaterialPropsReader(units_)
    {
	// ** empty **
    }
	
    
    // MANIPULATORS
    
    // ** none **
    
    // ACCESSORS

    virtual bool getNextMaterial(int materialId_, std::string &name);

    virtual void getTemperatureGrid(int materialId,
				    std::vector<double> &tempGrid_);
    
    virtual void getDensityGrid(int materialId,
				std::vector<double> &densityGrid_);

    virtual void getNumGroups(int materialId, int &numGroups);


    virtual void getEnergyUpperbounds(int materialId, int group,
				      double &energyUpperbounds_);

    virtual void getEnergyLowerbounds(int materialId, int group,
				      double &energyLowerbounds_);

    virtual void getSigmaTotal(int materialId, int group,
			       dsxx::Mat2<double> &data);

    virtual void getSigmaAbsorption(int materialId, int group,
				    dsxx::Mat2<double> &data);

    virtual void getSigmaEmission(int materialId, int group,
				  dsxx::Mat2<double> &data);

    virtual void getElectronIonCoupling(int materialId,
					dsxx::Mat2<double> &data);
	
    virtual void getElectronConductionCoeff(int materialId,
					    dsxx::Mat2<double> &data);
	
    virtual void getIonConductionCoeff(int materialId,
				       dsxx::Mat2<double> &data);
	
    virtual void getElectronSpecificHeat(int materialId,
					 dsxx::Mat2<double> &data);
	
    virtual void getIonSpecificHeat(int materialId,
				    dsxx::Mat2<double> &data);

  private:
    
    // IMPLEMENTATION
};


END_NS_XTM  // namespace XTM

#endif                          // __matprops_FifiMatPropsReader_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/FifiMatPropsReader.hh
//---------------------------------------------------------------------------//
