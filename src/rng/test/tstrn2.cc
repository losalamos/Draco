//----------------------------------*-C++-*----------------------------------//
// testrn2.cc
// Thomas M. Evans
// Thu Apr 30 13:33:20 1998
//---------------------------------------------------------------------------//
// @> Test of SPRNG and Sprng random number class.
//---------------------------------------------------------------------------//

#include "rng/Random.hh"
#include <vector>
#include <iostream>

using namespace rtt_rng;
using namespace std;

void t1()
{
    cout.setf(ios::floatfield, ios::fixed);
    cout.precision(4);
    cout << "<< TEST 1 >>" << endl << endl;

    Rnd_Control rcon(493875348);
    vector<Sprng> bucket;

    {
	Sprng ran  = rcon.get_rn();
	
	cout << rcon.get_num() << " " << ran.get_num() << endl;
	for (int i = 0; i < 5; i++)
	    cout << ran.ran() << endl;

	bucket.push_back(ran);
    }

    cout << "\nRandom Numbers from Bank" << endl << endl;

    Sprng ran2 = bucket.front();

    bucket.pop_back();

    for (int i = 2; i < 5; i++)
	cout << ran2.ran() << endl;

    ran2.print();

    cout << endl << "<< END TEST1 >>" << endl << endl;
}

void t2()
{
    cout << "<< TEST 2 >>" << endl << endl;

    Rnd_Control rcon(493875348);

    Sprng r1 = rcon.get_rn();
    Sprng r2 = rcon.get_rn();
    Sprng r3 = r1;

    r2 = r1;

    cout << r1.ran() << " " << r2.ran() << " " << r3.ran() << endl;

    cout << endl << "<< END TEST2 >>" << endl << endl;
}

void t3()
{
    cout << "<< TEST 3>>" << endl << endl;
    
    Rnd_Control rcon(493875348);

    Sprng r1 = rcon.get_rn();
    Sprng r2 = rcon.get_rn();
    Sprng rs = rcon.spawn(r1);
    Sprng rt = rcon.spawn(r1);
    
    for (int i = 0; i < 5; i++)
	cout << r1.ran() << "\t" << rs.ran() << "\t" << rt.ran() 
             << endl;

    for (int i = 0; i < 463; i++)
        Sprng rspawn = rcon.spawn(r1);
    cout << "Spawn 1 complete" << endl;
    for (int i = 0; i < 465; i++)
	Sprng rspawn = rcon.spawn(r2);
    cout << "Spawn 2 complete" << endl;
    for (int i = 0; i < 464; i++)
	Sprng rspawn = rcon.spawn(rs);
    cout << "Spawn 3 complete" << endl;
    cout << "If no warnings then TEST 3 succeeded" << endl;
}

void t4()
{
    cout << "<< TEST 4>>" << endl << endl;
    Rnd_Control rcon(493875348);
    cout << rcon.get_size() << endl;
}
	
main()
{
    t1();
    t2();
    t3();
    t4();
}

//---------------------------------------------------------------------------//
//                              end of testrn2.cc
//---------------------------------------------------------------------------//
