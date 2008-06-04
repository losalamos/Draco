//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstSend_Receive.cc
 * \author Mike Buksas
 * \date   Tue Jun  3 14:19:33 2008
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "../Release.hh"
#include "c4_test.hh"

#include "ds++/Packing_Utils.hh"
#include "../Send_Receive.hh"

using namespace std;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// Implementation classes
//---------------------------------------------------------------------------//

struct Send_Double_Vector : public Sender
{

    Send_Double_Vector(int node) : Sender(node) {  }

    void send(const vector<double>& v)
    {
        Sender::send(v.size(), &v[0]);
    }

};

struct Receive_Double_Vector : public Receiver
{

    Receive_Double_Vector(int node) : Receiver(node) {  }

    vector<double> receive()
    {
        int size = receive_size();
        vector<double> v(size);
        receive_data(size, &v[0]);
        return v;
    }



};

    




//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void auto_communication_test()
{

    Check(nodes() == 1);

    Send_Double_Vector sdv(0);
    Receive_Double_Vector rdv(0);

    vector<double> v(3);
    v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;

    sdv.send(v);
    vector<double> v2 = rdv.receive();

    if (v2.size() != 3) ITFAILS;
    if (v2[0] != v[0]) ITFAILS;
    if (v2[1] != v[1]) ITFAILS;
    if (v2[2] != v[2]) ITFAILS;

    if (rtt_c4_test::passed)
        PASSMSG("Passed auto communication test");

}

    
void single_comm_test()
{

    Check(nodes() == 2);

    if (node() == 0)
    {
        Send_Double_Vector sdv(1);

        vector<double> v(3);
        v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;

        sdv.send(v);

        sdv.wait();

    }

    if (node() == 1)
    {
        Receive_Double_Vector sdv(0);

        vector<double> v = sdv.receive();

        if (v.size() != 3) ITFAILS;
        if (v[0] != 1.0) ITFAILS;

    }

    if (rtt_c4_test::passed)
        PASSMSG("Passed single communication test");
}


void double_comm_test()
{

    Check(nodes() == 2);

    // Make data vectors of differeing size:
    vector<double> data((node()+1)*2);
    for (int i = 0; i < data.size(); ++i)
        data[i] = static_cast<double>(i);
    
    // Make a sender and receiver on each node to/from the other node.
    Send_Double_Vector sdv(1-node());
    Receive_Double_Vector rdv(1-node());
    
    // Send, receive and wait.
    sdv.send(data);
    vector<double> r = rdv.receive();

    if (r.size() != (2-node())*2) ITFAILS;

    for (int i = 0; i < r.size(); ++i)
        if (r[i] != static_cast<double>(i)) ITFAILS;

    if (rtt_c4_test::passed)
        PASSMSG("Passed double communication test");
}
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (std::string(argv[arg]) == "--version")
        {
            if (rtt_c4::node() == 0)
                cout << argv[0] << ": version " 
                     << rtt_c4::release() 
                     << endl;
            rtt_c4::finalize();
            return 0;
        }

    try
    {
        // >>> UNIT TESTS
        if (rtt_c4::nodes() == 1)
            auto_communication_test();
        
        if (rtt_c4::nodes() == 2)
        {
            single_comm_test();
            double_comm_test();
        }
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstSend_Receive, " 
                  << err.what()
                  << std::endl;
        rtt_c4::finalize();
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstSend_Receive, " 
                  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << std::endl;
        rtt_c4::finalize();
        return 1;
    }

    {
        rtt_c4::HTSyncSpinLock slock;

        // status of test
        std::cout << std::endl;
        std::cout <<     "*********************************************" 
                  << std::endl;
        if (rtt_c4_test::passed) 
        {
            std::cout << "**** tstSend_Receive Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
        }
        std::cout <<     "*********************************************" 
                  << std::endl;
        std::cout << std::endl;
    }
    
    rtt_c4::global_barrier();

    std::cout << "Done testing tstSend_Receive on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstSend_Receive.cc
//---------------------------------------------------------------------------//
