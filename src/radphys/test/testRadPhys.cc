//----------------------------------*-C++-*----------------------------------//
// testRadPhys.cc
// Randy M. Roberts
// Thu May 14 11:07:34 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "../RadiationPhysics.hh"
#include "test_utils.hh"

#include <iostream>

#include <vector>
using std::vector;

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(*x))

namespace rtt_radphys
{

 using namespace XTM;

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
	 
	 std::cout << "Stefan-Boltzmann: " << sb
		   << " ---> " << testMsg(compare_reals(sb, 1.02828, 5))
		   << std::endl;

	 double a = radphys.getRadConstant();
	 std::cout << "Radiation constant: " << a
		   << " ---> " << testMsg(compare_reals(a, 0.01372, 5))
		   << std::endl;

	 std::cout << "Planckian (vector):" ;
	 radphys.getPlanck(vTelectron, planckian);
	 
	 double known_planckian[] = {3.27313e-13, 8.37922e-11, 0.327313};

	 bool testPassed = true;
	 for (int i=0; i<planckian.size(); i++)
	 {
	     testPassed = testPassed && compare_reals(planckian[i],
						      known_planckian[i], 5);
	     std::cout << " " << planckian[i];
	 }
	 std::cout << " ---> " << testMsg(testPassed);
     }
     break;
     case 1:
     {
	 radphys.getPlanckTemperatureDerivative(vTelectron, planckian);
	 double answer[] = {1.30925e-09, 8.37922e-08, 1.30925};
	 std::cout << "Planckian temperature derivative (vector):" ;
	 bool testPassed = true;
	 for (int i=0; i<planckian.size(); i++)
	 {
	     testPassed = testPassed && compare_reals(planckian[i],
						      answer[i], 5);
	     std::cout << " " << planckian[i];
	 }
	 std::cout << " ---> " << testMsg(testPassed);
	 break;
     }
     case 2:
     {
	 double answer[] = {0.000749646, 0.00017076, 1.20839e-07, 0.0466896,
			    0.0058362, 9.25243e-06, 0.0489129, 0.00611412,
			     9.65704e-06};
	 double *ans = answer;
	 std::cout << "Electron-Ion Coupling (scalars):" ;
	 bool testPassed = true;
	 for (int i=0; i<vdensity.size(); i++)
	 {
	     for (int j=0; j<vTelectron.size(); j++)
	     {
		 double z = 5.0e-02;
		 double abar = 25.0; // BeO
		 double electIonCoupling = -999.9;
		
		 radphys.getElectIonCoupling(vdensity[i], vTelectron[j],
					     z, abar, electIonCoupling);
		 testPassed = testPassed && compare_reals(electIonCoupling,
							  *ans++, 5);
		 std::cout << " " << electIonCoupling;
	     }
	 }
	 std::cout << " ---> " << testMsg(testPassed);
	 break;
     }
     case 3:
     {
	 double answer[] = {.000749646, 0.0058362, 9.65704e-06};
	 double *ans = answer;
	 bool testPassed = true;
	 std::cout << "Electron-Ion Coupling (vectors):" ;
	 {
	     vector<double> z(vdensity.size(), 5.0e-2);
	     double abar = 25.0;
	     vector<double> electIonCoupling(vdensity.size());
	     radphys.getElectIonCoupling(vdensity, vTelectron,
					 z, abar, electIonCoupling);
	     for (int i=0; i<vdensity.size(); i++)
	     {
		 testPassed = testPassed && compare_reals(electIonCoupling[i],
							  *ans++, 5);
		 std::cout << " " << electIonCoupling[i];
	     }
	 }
	 std::cout << " ---> " << testMsg(testPassed);
	 break;
     }
     case 4:
     {
	 double answer[] = {1.68587, 3.07217, 8.59363, 1, 1, 6.26665,
			    1, 1, 6.24339};
	 double *ans = answer;
	 bool testPassed = true;
	 std::cout << "Electron-Ion Coulomb Log (scalars):" ;
	 for (int i=0; i<vdensity.size(); i++)
	 {
	     for (int j=0; j<vTelectron.size(); j++)
	     {
		 double z = 5.0e-02;
		 double abar = 25.0; // BeO
		 double electIonCoulombLog = -999.9;
		
		 radphys.getElectIonCoulombLog(vdensity[i], vTelectron[j],
					       z, abar, electIonCoulombLog);
		 testPassed = testPassed && compare_reals(electIonCoulombLog,
							  *ans++, 5);
		 std::cout << " " << electIonCoulombLog;
	     }
	 }
	 std::cout << " ---> " << testMsg(testPassed);
	 break;
     }
     case 5:
     {
	 double answer[] = {1.68587, 1, 6.24339};
	 double *ans = answer;
	 bool testPassed = true;
	 std::cout << "Electron-Ion Coulomb Log (vectors):" ;
	 {
	     vector<double> z(vdensity.size(), 5.0e-2);
	     double abar = 25.0;
	     vector<double> electIonCoulombLog(vdensity.size());
	     radphys.getElectIonCoulombLog(vdensity, vTelectron,
					   z, abar, electIonCoulombLog);
	     for (int i=0; i<vdensity.size(); i++)
	     {
		 testPassed = testPassed && compare_reals(electIonCoulombLog[i],
							  *ans++, 5);
		 std::cout << " " << electIonCoulombLog[i];
	     }
	 }
	 std::cout << " ---> " << testMsg(testPassed);
	 break;
     }
     }
     std::cout << std::endl;
 }

}

using std::cout;
using std::endl;

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main(int argc, char *argv[])
{
    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    return 0;
	}
    }

    const int nTests = 6;
    for (int i=0; i<nTests; i++)
    {
	try
	{
	    rtt_radphys::testRadPhys(i);
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
