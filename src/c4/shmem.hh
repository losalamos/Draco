//----------------------------------*-C++-*----------------------------------//
// shmem.hh
// Geoffrey Furnish
// 25 February 1998
//---------------------------------------------------------------------------//
// @> Header for implementation details of the SHMEM layer.
//---------------------------------------------------------------------------//

#ifndef __c4_shmem_hh__
#define __c4_shmem_hh__

#include <pthread.h>

#include "c4/tags.hh"

#include "ds++/Assert.hh"
#include "ds++/DynArray.hh"

C4_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Debugging stuff.

void C4_shm_dbg_1();

//---------------------------------------------------------------------------//
// C4_shmem API.

void C4_Wait( pthread_t& tid );

//---------------------------------------------------------------------------//
// Miscelaneous stuff.

//---------------------------------------------------------------------------//
// Utility classes, needed by the shmem interface.

class Msg_DB {
// Anybody dumb enough to mess with these deserves what they get :-).
  public:
    int source_pe;
    int tag;
    char *buf;
    int buf_size;
    int msg_size;

  public:
    Msg_DB( int _bfsz )
	: buf_size(_bfsz)
    {
	buf = new char[ buf_size ];
	msg_size = 0;
    }
    ~Msg_DB() { delete[] buf; }

    int matches( int t ) const { return t == C4_Any_Tag || t == tag; }
};

class Msg_Queue {

    int nmsgs;
    DynArray<Msg_DB *> msgs;

    int nfree;
    DynArray<Msg_DB *> freelist;

  public:
    Msg_Queue()
	: nmsgs(0), nfree(0)
    {}

    ~Msg_Queue()
    {
	if (nmsgs) {
	// Hmmm.  Program ending with queued messages?  Sounds fishy.

	// Delete all messages.

	    for( int i=0; i < nmsgs; i++ )
		delete msgs[i];
	}

	if (nfree) {

	    for( int i=0; i < nfree; i++ )
		delete freelist[i];
	}
    }

    int is_empty() const { return (nmsgs == 0); }
    int num_msgs() const { return nmsgs; }

    Msg_DB *msg( int i )
    {
	Assert( i >= 0 && i < nmsgs );
	return msgs[i];
    }

    void enqueue( Msg_DB *pm )
    {
	msgs[nmsgs++] = pm;
    }

    Msg_DB *new_msg( int msglen )
    {
	Msg_DB *pm = NULL;

    // First check to see if we have a suitably sized Msg_DB on the free
    // list.  If so, return it, and shorten the free list.

	for( int i=0; i < nfree; i++ ) {
	    if (freelist[i]->buf_size >= msglen) {
		pm = freelist[i];

		for( int j=i; j < nfree-1; j++ )
		    freelist[j] = freelist[j+1];

		    nfree--;
		    break;
	    }
	}

    // If not, make a new Msg_DB and return that.

	if (!pm)
	    pm = new Msg_DB( msglen );

    // Set the size as requested.

	pm->msg_size = msglen;

	return pm;
    }

    void discard( int i )
    {
	Assert( i >= 0 && i < nmsgs );

    // Move message i to the free list.
	freelist[nfree++] = msgs[i];

    // Contract the message list.
	for( int j=i; j < nmsgs-1; j++ )
	    msgs[j] = msgs[j+1];

	nmsgs--;
    }
};

enum AsyncMsgType {
    Async_Send,
    Async_Recv
};

enum AsyncMsgState {
    Inactive,
    Setting_Up,
    Pending,
    No_Longer_Pending,		// So says the remote node.
    Complete			// So says the local node.
};

struct Async_DB {
    AsyncMsgState state;
    void *data;
    int tag;
    int size;
};

extern int C4_shm_mynode;
extern int C4_shm_nodes;

class C4_Req;

int SHM_Send( void *buf, int size, int dest, int tag, int group );
int SHM_Recv( void *buf, int size, int dest, int tag, int group );

C4_Req SHM_SendAsync( void *buf, int size, int dest, int tag, int group );
C4_Req SHM_RecvAsync( void *buf, int size, int source, int tag, int group );

C4_NAMESPACE_END

#endif                          // __c4_shmem_hh__

//---------------------------------------------------------------------------//
//                              end of c4/shmem.hh
//---------------------------------------------------------------------------//
