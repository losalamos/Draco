//----------------------------------*-C++-*----------------------------------//
// tgsync.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//
// @> Test program for C4_gsync().
//---------------------------------------------------------------------------//

#include "../global.hh"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include <unistd.h>

using std::cout;
using std::endl;

// if(mype == 0) {
//
//     for(i=1;i<npes;i++) get message from other process
//    
//     n = 0;
//     while(n<ntimes) {
//        iprobe for additional messages
//        if(iprobe returns true) abort();
//        n++;
//     }
//
//     MPI_Barrier;
//
//     for(i=1;i<npes;i++) {
//       get message from other process
//       if(!second message) abort();
//     }
//
// } else {
//
//    send first message to process 0
//    MPI_Barrier
//    send second message to process 0
// }

// unnamed namespace to keep from clutering up the global namespace

namespace
{

std::string FailReason;
    
void version(const std::string &progname)
{
    if (C4::node() == 0)
    {
	std::string version = "1.0.0";
	cout << progname << ": version " << version << endl;
    }
}

inline bool everyProcReceived(int nodes, int expectedValue)
{
    // everyProcReceived posts a receive for each non-zero process.
    // It then makes sure that all have them have completed.
    // If not, it frees the posted requests.
    
    int nReceived = 0;
    std::vector<bool> received(nodes-1, false);
    std::vector<C4::C4_Req> reqs(nodes-1);
    
    int buf;

    // Post the receives.
    
    for (int ip=1; ip<nodes; ip++)
	reqs[ip-1] = C4::RecvAsync<int>(&buf, 1, ip);

    // Check to make sure all of the receives have completed with the
    // expected value.
    
    while (nReceived < nodes-1)
    {
	for (int ip=1; ip<nodes; ip++)
	{
	    if (!received[ip-1] && reqs[ip-1].complete())
	    {
		reqs[ip-1].wait();
		nReceived++;
		cout << "Recieved node: " << ip << " value: " << buf << endl;
		cout << "nReceived: " << nReceived << endl;
		if (buf != expectedValue)
		{
		    cout << "Recieved value not expected." << endl;
		    FailReason = "Recieved value not expected.";
		    return false;
		}
		received[ip-1] = true;
	    }
	}
    }

    cout << "Out of while loop." << endl;

    // Test to make sure each non-zero process received at least once.

    bool testPassed = true;
    
    for (int ip=1; ip<nodes; ip++)
    {
	if (!received[ip-1])
	{
	    reqs[ip-1].free();
	    FailReason = "Not All Messages Received";
	    testPassed = false;
	}
    }
    
    return testPassed;
}

inline bool anyMoreReceives(int nodes)
{
    // See if any extra receives come in.
    
    bool moreReceives = false;
    
    int buf;
    std::vector<bool> received(nodes-1, false);
    std::vector<C4::C4_Req> reqs(nodes-1);

    for (int ip=1; ip<nodes; ip++)
	reqs[ip-1] = C4::RecvAsync<int>(&buf, 1, ip);
    
    const int nTimes = 1000;

    cout << "Begin checking for extra receives "
	 << nTimes << " times." << endl;
    
    for (int it=0; it<nTimes; it++)
    {
	// cout << "Check #" << it << endl;
	
	for (int ip=1; ip<nodes; ip++)
	{
	    if (!received[ip-1] && reqs[ip-1].complete())
	    {
		received[ip-1] = true;
		moreReceives = true;
	    }
	}
    }

    if (moreReceives)
	cout << "There are extra receives." << endl;
    else
	cout << "No extra receives." << endl;

    // Free up the requests.
    
    for (int ip=1; ip<nodes; ip++)
	reqs[ip-1].free();
    
    return moreReceives;
}

// Semi-random values to send/receive

const int FirstVal = 99834;
const int SecondVal = 224773;

inline bool zeroNodeDoitOK(int nodes)
{

    // Make sure every process has been received with the FirstVal.
    
    if (!everyProcReceived(nodes, FirstVal))
	return false;
    
    // We should not get any more receives

    if (anyMoreReceives(nodes))
    {
	FailReason = "Too Many Receives.";
	return false;
    }

    // This is what we are testing... gsync()!!!!
    
    C4::gsync();

    // Make sure every process has been received with the SecondVal.
    
    if (!everyProcReceived(nodes, SecondVal))
	return false;

    return true;

} // end unamed namespace

inline void otherNodesDoit(int node)
{
    cout << "Node: " << node << " sending " << FirstVal << "." << endl;
    
    int buf = FirstVal;
    C4::C4_Req r1 = C4::SendAsync<int>(&buf, 1, 0);

    // This is what we are testing... gsync()!!!!
    
    C4::gsync();

    cout << "Node: " << node << " sending " << SecondVal << "." << endl;

    buf = SecondVal;
    C4::C4_Req r2 = C4::SendAsync<int>(&buf, 1, 0);

    r1.wait();
    r2.wait();
}

} // end unnamed namespace

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    C4::Finalize();
	    return 0;
	}
    }

    int node = C4::node();
    int nodes = C4::nodes();

    cout << "Begining Test." << endl;
    
    if (node == 0)
    {
	if (zeroNodeDoitOK(nodes))
	    cout << "tgsync: Test: Passed." << endl;
	else
	    cout << "tgsync: Test: Failed: " << FailReason << "." << endl;
    }
    else
    {
	// sleep(5);
	otherNodesDoit(node);
    }

    cout << "Test Completed." << endl;
    
    C4::Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tgsync.cc
//---------------------------------------------------------------------------//
