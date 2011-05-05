//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstSend_Receive.cc
 * \author Mike Buksas
 * \date   Tue Jun  3 14:19:33 2008
 * \brief  
 * \note   Copyright (C) 2006-2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Send_Receive.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "ds++/Release.hh"
#include "c4_test.hh"

#include "ds++/Packing_Utils.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// Implementation classes
//---------------------------------------------------------------------------//

class Send_Double_Vector : public Sender
{
  public:
  
    explicit Send_Double_Vector(int node) 
      : Sender(node) 
    { 
      return; 
    }

    /*! \bug should this be const data?  If so the base class Sender must also
     *       mark the data as const. */
    void send(vector<double> const & v)
    {
        Sender::send(v.size(), (v.size()>0 ? &v[0] : NULL));
    }

    void wait()
    {
        Sender::wait();
    }

};

class Receive_Double_Vector : public Receiver
{
  public:
    explicit Receive_Double_Vector(int node) 
      : Receiver(node)
    { return; }

    vector<double> receive()
    {
        int size = Receiver::receive_size();
        vector<double> v(size,0.0);
        if( size > 0 )
            receive_data(size, &v[0]);

        return v;
    }

};

class Receive_Double_Vector_Autosize : public Receiver
{
  public:
  
    explicit Receive_Double_Vector_Autosize(int node) 
      : Receiver(node) 
    {
      return; 
    }

    vector<double> receive()
    {
        double *v(NULL);
        unsigned const size = Receiver::receive(v);
        vector<double> retVal(v, v+size);
        delete [] v;
        return retVal;
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

    Require(nodes() == 1);
    vector<double> v(3);
    v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;

    Send_Double_Vector sdv(0);
    sdv.send(v);

    Receive_Double_Vector rdv(0);
    vector<double> v2 = rdv.receive();
    Check( v2.size() == 3 );

    if (v2.size() != 3) ITFAILS;
    if (v2[0] != v[0]) ITFAILS;
    if (v2[1] != v[1]) ITFAILS;
    if (v2[2] != v[2]) ITFAILS;

    if (rtt_c4_test::passed)
        PASSMSG("Passed auto communication test");
    return;
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
        sdv.wait();
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
        sdv.wait();
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

    size_t const other = 1-node();

    // Assign sizes and contents of the two vectors.
    size_t const sizes[2] = {4,7};
    vector<double> data(sizes[node()]);
    for( size_t i = 0; i < data.size(); ++i )
        data[i] = static_cast<double>(i);
    
    // Make a sender and receiver on each node to/from the other node.
    Send_Double_Vector sdv(other);
    Receive_Double_Vector rdv(other);
    
    // Send and receive.
    sdv.send(data);
    vector<double> r = rdv.receive();

    // Check the size and contents of the received vector.
    if (r.size() != sizes[other]) ITFAILS;

    for (size_t i = 0; i < r.size(); ++i)
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

    size_t const sizes[] = {1, 4, 7, 10};
    vector<double> data(sizes[node()]);

    for (size_t i=0; i<data.size(); ++i)
        data[i] = static_cast<double>(i*i);

    Send_Double_Vector sender(to_node);
    Receive_Double_Vector receiver(from_node);

    sender.send(data);
    vector<double> r = receiver.receive();

    if (r.size() != sizes[from_node]) ITFAILS;
    for (size_t i=0; i<r.size(); ++i)
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
                     << rtt_dsxx::release() 
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
