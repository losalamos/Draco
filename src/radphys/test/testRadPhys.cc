//----------------------------------*-C++-*----------------------------------//
// testRadPhys.cc
// Randy M. Roberts
// Thu May 14 11:07:34 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "../RadiationPhysics.hh"

#include <iostream>

#include <vector>
using std::vector;

using namespace XTM;

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(*x))

void testRadPhys(int testnum)
{
    Units units = Units::getAstroPhysUnits();
    
    double Telectron[] = {1.0e-3, 4.0e-3, 1.0};
    vector<double> vTelectron(Telectron, Telectron+ARRAY_SIZE(Telectron));

    vector<double> planckian(vTelectron.size());

    double density[] = {2.0e-02, 2.1, 2.2};
    vector<double> vdensity(density, density+ARRAY_SIZE(density));

    RadiationPhysics radphys(units);

    switch (testnum)
    {
    case 0:
    {
	double sb = radphys.getStefanBoltzmann();
	std::cout << "Stefan-Boltzmann: " << sb << std::endl;

	double a = radphys.getRadConstant();
	std::cout << "Radiation constant: " << a << std::endl;

	std::cout << "Planckian:" << std::endl;
	radphys.getPlanck(vTelectron, planckian);
	for (int i=0; i<planckian.size(); i++)
	    std::cout << " " << planckian[i];
    }
    break;
    case 1:
	radphys.getPlanckTemperatureDerivative(vTelectron, planckian);
	std::cout << "Planckian temperature derivateve:" << std::endl;
	for (int i=0; i<planckian.size(); i++)
	    std::cout << " " << planckian[i];
	break;
    case 2:
	std::cout << "Electron-Ion Coupling (scalars):" << std::endl;
	for (int i=0; i<vdensity.size(); i++)
	{
	    for (int j=0; j<vTelectron.size(); j++)
	    {
		double z = 5.0e-02;
		double abar = 25.0; // BeO
		double electIonCoupling = -999.9;
		
		radphys.getElectIonCoupling(vdensity[i], vTelectron[j],
					    z, abar, electIonCoupling);
		std::cout << " " << electIonCoupling;
	    }
	}
	break;
    case 3:
	std::cout << "Electron-Ion Coupling (vectors):" << std::endl;
	{
	    vector<double> z(vdensity.size(), 5.0e-2);
	    double abar = 25.0;
	    vector<double> electIonCoupling(vdensity.size());
	    radphys.getElectIonCoupling(vdensity, vTelectron,
					z, abar, electIonCoupling);
	    for (int i=0; i<vdensity.size(); i++)
		std::cout << " " << electIonCoupling[i];
	}
	break;
    case 4:
	std::cout << "Electron-Ion Coulomb Log (scalars):" << std::endl;
	for (int i=0; i<vdensity.size(); i++)
	{
	    for (int j=0; j<vTelectron.size(); j++)
	    {
		double z = 5.0e-02;
		double abar = 25.0; // BeO
		double electIonCoulombLog = -999.9;
		
		radphys.getElectIonCoulombLog(vdensity[i], vTelectron[j],
					    z, abar, electIonCoulombLog);
		std::cout << " " << electIonCoulombLog;
	    }
	}
	break;
    case 5:
	std::cout << "Electron-Ion Coulomb Log (vectors):" << std::endl;
	{
	    vector<double> z(vdensity.size(), 5.0e-2);
	    double abar = 25.0;
	    vector<double> electIonCoulombLog(vdensity.size());
	    radphys.getElectIonCoulombLog(vdensity, vTelectron,
					  z, abar, electIonCoulombLog);
	    for (int i=0; i<vdensity.size(); i++)
		std::cout << " " << electIonCoulombLog[i];
	}
	break;
    }
    std::cout << std::endl;
}

using std::cout;
using std::endl;

int main()
{
    const int nTests = 6;
    for (int i=0; i<nTests; i++)
    {
	try
	{
	    testRadPhys(i);
	}
	catch (const dsxx::assertion &ass)
	{
	    std::cerr << ass.what() << std::endl;
	}

	cout << endl;
    }
}

//---------------------------------------------------------------------------//
//                              end of testRadPhys.cc
//---------------------------------------------------------------------------//6
