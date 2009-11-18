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
#include "../global.hh"
#include "../SpinLock.hh"
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

    void wait()
    {
        Sender::wait();
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

struct Receive_Double_Vector_Autosize : public Receiver
{

    Receive_Double_Vector_Autosize(int node) : Receiver(node) {  }

    vector<double> receive()
    {
        double *v;
        unsigned const size = Receiver::receive(v);
        return vector<double>(v, v+size);
    }

};


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*
 * A single node communicates with itself.
 */
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

    
//---------------------------------------------------------------------------//
/* 
 * One way communication from node 0 to 1.
 */
void single_comm_test()
{

    Check(nodes() == 2);

    if (node() == 0)
    {
        vector<double> v(3);
        v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;

        Send_Double_Vector sdv(1);
        sdv.send(v);

        // Check zero length branch
        
        sdv.wait();
        v.clear();
        sdv.send(v);
    }

    if (node() == 1)
    {
        Receive_Double_Vector sdv(0);
        vector<double> v = sdv.receive();

        if (v.size() != 3) ITFAILS;
        if (v[0] != 1.0) ITFAILS;
        if (v[1] != 2.0) ITFAILS;
        if (v[2] != 3.0) ITFAILS;

        // Check zero length branch
        v = sdv.receive();

        if (v.size() != 0) ITFAILS;
    }

    if (rtt_c4_test::passed)
        PASSMSG("Passed single communication test");
}

    
//---------------------------------------------------------------------------//
/* 
 * One way communication from node 0 to 1, autosized receive.
 */
void single_comm_autosize_test()
{

    Check(nodes() == 2);

    if (node() == 0)
    {
        vector<double> v(3);
        v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;

        Send_Double_Vector sdv(1);
        sdv.send(v);

        // Check zero length branch
        
        sdv.wait();
        v.clear();
        sdv.send(v);
    }

    if (node() == 1)
    {
        Receive_Double_Vector_Autosize sdv(0);
        vector<double> v = sdv.receive();

        if (v.size() != 3) ITFAILS;
        if (v[1] != 2.0) ITFAILS;
        if (v[2] != 3.0) ITFAILS;

        // Check zero length branch
        v = sdv.receive();

        if (v.size() != 0) ITFAILS;
    }

    if (rtt_c4_test::passed)
        PASSMSG("Passed single communication test");
}


//---------------------------------------------------------------------------//
/*
 * Nodes 0 and 1 both send and receive data from the other.
 * 
 */
void double_comm_test()
{

    Check(nodes() == 2);

    const int other = 1-node();

    // Assign sizes and contents of the two vectors.
    const int sizes[2] = {4,7};
    vector<double> data(sizes[node()]);
    for (int i = 0; i < data.size(); ++i) data[i] = static_cast<double>(i);
    
    // Make a sender and receiver on each node to/from the other node.
    Send_Double_Vector sdv(other);
    Receive_Double_Vector rdv(other);
    
    // Send and receive.
    sdv.send(data);
    vector<double> r = rdv.receive();

    // Check the size and contents of the received vector.
    if (r.size() != sizes[other]) ITFAILS;

    for (int i = 0; i < r.size(); ++i)
        if (r[i] != static_cast<double>(i)) ITFAILS;

    if (rtt_c4_test::passed)
        PASSMSG("Passed double communication test");
}



//---------------------------------------------------------------------------//
/*
 * Four nodes pass data to the right.
 * 
 */
void ring_test()
{

    Check(nodes() == 4);

    const int to_node   = (node()+1) % nodes();
    const int from_node = (node()-1+nodes()) % nodes();

    const int sizes[] = {1, 4, 7, 10};
    vector<double> data(sizes[node()]);

    for (int i=0; i<data.size(); ++i) data[i] = static_cast<double>(i*i);

    Send_Double_Vector sender(to_node);
    Receive_Double_Vector receiver(from_node);

    sender.send(data);
    vector<double> r = receiver.receive();

    if (r.size() != sizes[from_node]) ITFAILS;
    for (int i=0; i<r.size(); ++i)
        if (r[i] != static_cast<double>(i*i)) ITFAILS;

    if (rtt_c4_test::passed)
        PASSMSG("Passed ring communication test");

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
#ifndef C4_SCALAR
        if (rtt_c4::nodes() == 1)
            auto_communication_test();
#endif
        
        if (rtt_c4::nodes() == 2)
        {
            single_comm_test();
            single_comm_autosize_test();
            double_comm_test();
        }

        if (rtt_c4::nodes() == 4)
        {
            ring_test();
        }
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstSend_Receive, " 
                  << err.what()
                  << std::endl;
        rtt_c4::abort();
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstSend_Receive, " 
                  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << std::endl;
        rtt_c4::abort();
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
