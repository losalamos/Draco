//----------------------------------*-C++-*----------------------------------//
// testasync.cc
// Thomas M. Evans
// Sat Aug  1 11:58:27 1998
//---------------------------------------------------------------------------//
// @> Test of async com
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include <iostream>

using namespace std;
using namespace C4;

void nonblock()
{
    if (!node())
	cout << endl << ">>> NONBLOCKING TEST <<<" << endl;
  // make C4 objects
    C4_Req recv;
    C4_Req send;

    int r[4];
    bool receiving = true;

  // post a receive on node 1
    if (node())
	RecvAsync (recv, &r[0], 4, 0, 10);

  // make data on node 0
    if (!node())
    {
	int x[4] = {1,2,3,4};
	int y[4] = {1,2,4,6};
	cout << "node 0 is iterating" << endl;
	for (int i = 0; i < 1000000; i++);
	cout << "node 0 is sending" << endl;
	SendAsync(send, &x[0], 4, 1, 10);
	send.wait();
	cout << "0 has sent x" << endl;
	SendAsync(send, &y[0], 4, 1, 10);
	send.wait();
	cout << "0 has sent y" << endl;
	int f[4] = {-1};
	SendAsync(send, &f[0], 4, 1, 10);
	send.wait();
	cout << "0 is done" << endl;
    }

  // recieve data on node 1
    if (node())
    {
	int ctr;
	while (receiving)
	{
	    ctr++;
	    if (recv.complete())
	    {
		if (r[0] >= 0)
		{
		    cout << "1 has received data" << endl;
		    for (int i = 0; i < 4; i++)
			cout << "r["  << i << "]" << " = " << r[i] << endl;
		    RecvAsync (recv, &r[0], 4, 0, 10);
		}
		else if (r[0] < 0)
		{
		    cout << "1 has received a terminating message" << endl;
		    receiving = false;
		}
	    }	
	}
	cout << "1 waited " << ctr << " iterations" << endl;
    }

  // lets get together at end of this test
    gsync();
    cout << node() << " is done with async communication!" << endl;
    gsync();
}

void block()
{
    if (!node())
	cout << endl << "<<< BLOCKING TEST >>>" << endl;

    int x;
    if (!node())
    {
	int s = 5;
	Send (s, 1, 100);
	cout << "0 has sent a blocking send of " << s << endl;
    }
    if (node())
    {
	Recv (x, 0, 100);
	cout << "1 has received a blocking send x = " << x << endl;
    }
    
  // lets get together at end of this test
    gsync();
    cout << node() << " is done with sync communication!" << endl;
    gsync();
}

void reorder()
{
    if (!node())
	cout << endl << "<<< SWITCH ORDER TEST >>>" << endl;

    int x;
    int y[4];

    if (!node())
    {
	x = 5;
	y[0] = 1;
	y[1] = 2;
	y[2] = 3;
	y[3] = 4;
	Send (x, 1, 500);
	cout << "0 has sent x = " << x << " to 1" << endl;
	Send (&y[0], 4, 1, 501);
	cout << "0 has sent y = {" << y[0] << "," << y[1] << ","
	     << y[2] << "," << y[3] << "} to 1" << endl;
    }
    
    if (node())
    {
	Recv (&y[0], 4, 0, 501);
	cout << "1 has received y = {" << y[0] << "," << y[1] << ","
	     << y[2] << "," << y[3] << "}" << endl;
	Recv (x, 0, 500);
	cout << "1 has received x = " << x << endl;
    }

  // lets get together
    gsync();
    cout << node() << " is done with reordered communication!" << endl;
    gsync();
}

void mix()
{    
    if (!node())
	cout << endl << "<<< MIXED TEST >>>" << endl;

  // get async handles
    C4_Req send;
    C4_Req recv;

  // data
    int r[] = {1,2,3,4};

  // post a receive on 1
    if (node())
	RecvAsync(recv, &r[0], 4, 0, 500);

  // mangle data on 0
    if (!node())
    {
	int x[4];
	int y[4];
	int message;
	SendAsync(send, &x[0], 4, 1, 500);
	send.wait();
	SendAsync(send, &y[0], 4, 1, 500);
	send.wait();
	cout << "0 has async sent all its data" << endl;
	Recv (message, 1, 501);
	cout << "0 has sync received the message " << message << endl;
    }

  // get data on 1 
    if (node())
    {
	recv.wait();
	cout << "1 has async received from 0" << endl;
	RecvAsync(recv, &r[0], 4, 0, 500);
	recv.wait();
	cout << "1 has async received from 0" << endl;
	int message = 10;
	Send (message, 0, 501);
	cout << "1 has sync sent a message to 0" << endl;
    }

  // lets get together at end of this test
    gsync();
    cout << node() << " is done with mix communication!" << endl;
    gsync();
}

void reduce_scalar()
{
   if (!node())
       cout << endl << "<<< REDUCTION TEST 1 >>>" << endl;

   gsync();

 // value on each processor
   int x = 10 * (node() + 1);
   cout << "x on " << node() << " is " << x << endl;
   gsum(x);
   gsync();
   cout << "x on " << node() << " after reduction is " << x << endl;
   gsync();
}

void reduce_array()
{
   if (!node())
       cout << endl << "<<< REDUCTION TEST 2 >>>" << endl;

 // value on each processor
   int x[] = {1*(node()+1), 10*(node()+1)};
   cout << "x on " << node() << " is [" << x[0] << "," << x[1] << "]"
	<< endl;
   if (node() == 0)
       for (int i = 0; i < 10000000; i++);
   gsum(x,2);
   gsync();
   cout << "x on " << node() << " after reduction is [" << x[0] << "," 
	<< x[1] << "]" << endl;
   gsync();
}

int main(int argc, char *argv[])
{
  // C4 initialization
    Init(argc, argv);

  // lets try a nonblocking send receive
    nonblock();

  // lets try a blocking send receive
    block();

  // lets mix non-blocking and blocking sends
    mix();

  // lets try a reduction
    reduce_scalar();
    reduce_array();

  // lets mix up the order on blocking sends
    reorder();

  // C4 final
    Finalize();
}

//---------------------------------------------------------------------------//
//                              end of testasync.cc
//---------------------------------------------------------------------------//
