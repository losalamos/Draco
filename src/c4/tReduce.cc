//----------------------------------*-C++-*----------------------------------//
// tReduce.cc
// Maurice LeBrun
// Mon Feb  6 15:16:26 1995
//---------------------------------------------------------------------------//
// @> Test program for reduction operations
//---------------------------------------------------------------------------//

#include "c4/global.hh"

#include <iostream.h>

int node, nodes;

//---------------------------------------------------------------------------//

void tReduce_sum_int()
{
    int local_sum = 0, global_sum = 0;

    for (int i = 0; i < nodes; i++)
	local_sum += i;

    global_sum = node;
    C4_gsum(global_sum);

    if (node == 0) {
	if (local_sum == global_sum)
	    cout << "gsum for type int correct: " << global_sum << endl;
	else
	    cout << "gsum for type int bogus: " << global_sum << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_sum_long()
{
    long local_sum = 0, global_sum = 0;

    for (int i = 0; i < nodes; i++)
	local_sum += i;

    global_sum = node;
    C4_gsum(global_sum);

    if (node == 0) {
	if (local_sum == global_sum)
	    cout << "gsum for type long correct: " << global_sum << endl;
	else
	    cout << "gsum for type long bogus: " << global_sum << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_sum_float()
{
    float local_sum = 0, global_sum = 0;

    for (int i = 0; i < nodes; i++)
	local_sum += i;

    global_sum = node;
    C4_gsum(global_sum);

    if (node == 0) {
	if (local_sum == global_sum)
	    cout << "gsum for type float correct: " << global_sum << endl;
	else
	    cout << "gsum for type float bogus: " << global_sum << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_sum_double()
{
    double local_sum = 0, global_sum = 0;

    for (int i = 0; i < nodes; i++)
	local_sum += i;

    global_sum = node;
    C4_gsum(global_sum);

    if (node == 0) {
	if (local_sum == global_sum)
	    cout << "gsum for type double correct: " << global_sum << endl;
	else
	    cout << "gsum for type double bogus: " << global_sum << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_min_int()
{
    int local_min = -nodes+1, global_min;

    global_min = -node;
    C4_gmin(global_min);

    if (node == 0) {
	if (local_min == global_min)
	    cout << "gmin for type int correct: " << global_min << endl;
	else
	    cout << "gmin for type int bogus: " << global_min << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_min_long()
{
    long local_min = -nodes+1, global_min;

    global_min = -node;
    C4_gmin(global_min);

    if (node == 0) {
	if (local_min == global_min)
	    cout << "gmin for type long correct: " << global_min << endl;
	else
	    cout << "gmin for type long bogus: " << global_min << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_min_float()
{
    float local_min = -nodes+1, global_min;

    global_min = -node;
    C4_gmin(global_min);

    if (node == 0) {
	if (local_min == global_min)
	    cout << "gmin for type float correct: " << global_min << endl;
	else
	    cout << "gmin for type float bogus: " << global_min << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_min_double()
{
    double local_min = -nodes+1, global_min;

    global_min = -node;
    C4_gmin(global_min);

    if (node == 0) {
	if (local_min == global_min)
	    cout << "gmin for type double correct: " << global_min << endl;
	else
	    cout << "gmin for type double bogus: " << global_min << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_max_int()
{
    int local_max = nodes-1, global_max;

    global_max = node;
    C4_gmax(global_max);

    if (node == 0) {
	if (local_max == global_max)
	    cout << "gmax for type int correct: " << global_max << endl;
	else
	    cout << "gmax for type int bogus: " << global_max << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_max_long()
{
    long local_max = nodes-1, global_max;

    global_max = node;
    C4_gmax(global_max);

    if (node == 0) {
	if (local_max == global_max)
	    cout << "gmax for type long correct: " << global_max << endl;
	else
	    cout << "gmax for type long bogus: " << global_max << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_max_float()
{
    float local_max = nodes-1, global_max;

    global_max = node;
    C4_gmax(global_max);

    if (node == 0) {
	if (local_max == global_max)
	    cout << "gmax for type float correct: " << global_max << endl;
	else
	    cout << "gmax for type float bogus: " << global_max << endl;
    }
}

//---------------------------------------------------------------------------//

void tReduce_max_double()
{
    double local_max = nodes-1, global_max;

    global_max = node;
    C4_gmax(global_max);

    if (node == 0) {
	if (local_max == global_max)
	    cout << "gmax for type double correct: " << global_max << endl;
	else
	    cout << "gmax for type double bogus: " << global_max << endl;
    }
}

//---------------------------------------------------------------------------//
// Verify that global reduction operations work.
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    C4_Init( argc, argv );

    node = C4_node();
    nodes = C4_nodes();

    tReduce_sum_int();
    tReduce_sum_long();
    tReduce_sum_float();
    tReduce_sum_double();

    tReduce_min_int();
    tReduce_min_long();
    tReduce_min_float();
    tReduce_min_double();

    tReduce_max_int();
    tReduce_max_long();
    tReduce_max_float();
    tReduce_max_double();

    C4_Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tReduce.cc
//---------------------------------------------------------------------------//
