//---------------------------------------------------------------------------//
// @> Test of async communication with vectors
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include <iostream>
#include <vector>

using namespace std;
using namespace C4;

void from_master()
{
  // master
    if (!node())
    {
	vector<int> send_vector;
	send_vector.resize(nodes()-1);

	for (int i = 1; i < nodes(); i++)
	    send_vector[i-1] = 12*i;

	for (int i = 1; i < nodes(); i++)
	    Send(&send_vector[i-1], 1, i, 100);

	cout << " Master has sent the elements to the lesser nodes." 
	     << endl;
    }

    else if (node())
    {
	C4_Req r;
	int recv_element = -1;
	int rcount = 0;

	r = RecvAsync(&recv_element, 1, 0, 100);

	while (!r.complete())
	    rcount++;

	if (recv_element == 12*node())
	    cout << " Node " << node() 
		 << " received the proper element value with "
		 << rcount << " tests." << endl;
	else 
	    cout << " recv_element = " << recv_element << endl;
    }
}

void to_master()
{
  // master
    if (!node())
    {
	vector<int> recv_vector;
	recv_vector.resize(nodes()-1);

	vector<C4_Req> rmsg;
	rmsg.resize(nodes()-1);

	vector<int> received;
	received.resize(nodes()-1);

      // init recv vectors
	for (int i = 1; i < nodes(); i++)
	{
	    recv_vector[i-1] = -1;
	    received[i-1] = 0;
	}

      // post async recvs
	for (int i = 1; i < nodes(); i++)
	    rmsg[i-1] = RecvAsync(&recv_vector[i-1], 1, i, 100);
	
      // initialize counter
	int num_waiting_for = nodes()-1; 
	int loop_count = 0;

	while (num_waiting_for)
	{
	    loop_count++;
	    if (!(loop_count % 10000))
		cout << " master's loop count = " << loop_count 
		     << " received[1]= " << received[1] 
		     << " num_waiting_for = " << num_waiting_for 
		     << " rmsg[1] = " << rmsg[1].complete() << endl;
	    for (int i = 1; i < nodes(); i++)
	    {
		if (!received[i-1])
		{
		    if (rmsg[i-1].complete())
		    {
			num_waiting_for--;
			received[i-1] = 1;
		    }
		}
	    }
	}

	int check = 0;
	int calcd = 0;
	for (int i = 1; i < nodes(); i++)
	{
	    check += i*12;
	    calcd += recv_vector[i-1];
	}


	if (calcd == check)
	    cout << " MASTER received all " << nodes()-1 
		 << " nodes' values, requiring " << loop_count
		 << " loops over all nodes." << endl;

	else
	{
	    cout << "MASTER required " << loop_count
		 << " loops over all nodes." << endl;
	    cout << " Test failed.  Master's values follow:" << endl;
	    for (int i = 1; i < nodes(); i++)
		cout << "    Node " << i << ": " << recv_vector[i-1] 
		     << endl;
	}
    }

  // lesser nodes
    else if (node())
    {
	int send_data = 12*node();
	Send(&send_data, 1, 0, 100);
	cout << " I, node " << node() 
	     << ", have sent my data to the master." << endl;
    }
}

int main(int argc, char *argv[])
{
  // C4 initialization
    Init(argc, argv);

  // send from master out to other nodes
    if (!node())
    {
	cout << endl << endl;
	cout << " FIRST TEST: master sending out to all other nodes " 
	     << endl << endl;
    }
    from_master();

  // everyone regroup for the next test
    gsync();

  // all lesser nodes send to the master
    if (!node())
    {
	cout << endl << endl;
	cout << " SECOND TEST: all nodes sending to master " 
	     << endl << endl;
    }
    to_master();
    
  // C4 final
    Finalize();
}

