//----------------------------------*-C++-*----------------------------------//
// TestMatPropReader.hh
// Randy M. Roberts
// Mon Apr 20 15:55:22 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_test_TestMatPropReader_hh__
#define __matprops_test_TestMatPropReader_hh__

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
// class TestMatPropReader - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestMatPropReader : public MaterialPropsReader
{
    // NESTED CLASSES AND TYPEDEFS

    // DATA

    // DISALLOWED DEFAULT METHODS
    
  private:
    
    TestMatPropReader(const TestMatPropReader &rhs);
    TestMatPropReader& operator=(const TestMatPropReader &rhs);

  public:

    // CREATORS
    
    TestMatPropReader(const Units &units_)
	: MaterialPropsReader(units_)
    {
	// ** empty **
    }
	
    
    // MANIPULATORS
    
    std::istream &read(std::istream &is);
    
    // ACCESSORS

    virtual bool getNextMaterial(int &materialId_, std::string &name) const;

    virtual void getTemperatureGrid(int materialId,
				    std::vector<double> &tempGrid_) const;
    
    virtual void getDensityGrid(int materialId,
				std::vector<double> &densityGrid_) const;

    virtual void getNumGroups(int materialId, int &numGroups) const;


    virtual void getEnergyUpperbounds(int materialId, int group,
				      double &energyUpperbounds_) const;

    virtual void getEnergyLowerbounds(int materialId, int group,
				      double &energyLowerbounds_) const;

    virtual void getSigmaTotal(int materialId, int group,
			       dsxx::Mat2<double> &data) const;

    virtual void getSigmaAbsorption(int materialId, int group,
				    dsxx::Mat2<double> &data) const;

    virtual void getSigmaEmission(int materialId, int group,
				  dsxx::Mat2<double> &data) const;

    virtual void getElectronIonCoupling(int materialId,
					dsxx::Mat2<double> &data) const;
	
    virtual void getElectronConductionCoeff(int materialId,
					    dsxx::Mat2<double> &data) const;
	
    virtual void getIonConductionCoeff(int materialId,
				       dsxx::Mat2<double> &data) const;
	
    virtual void getElectronSpecificHeat(int materialId,
					 dsxx::Mat2<double> &data) const;
	
    virtual void getIonSpecificHeat(int materialId,
				    dsxx::Mat2<double> &data) const;

  private:
    
    // IMPLEMENTATION
};


// Convenience function.

inline std::istream &operator>>(std::istream &is, TestMatPropReader &rhs)
{
    return rhs.read(is);
}


END_NS_XTM  // namespace XTM

#endif                          // __matprops_test_TestMatPropReader_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/test/TestMatPropReader.hh
//---------------------------------------------------------------------------//
