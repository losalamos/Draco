//----------------------------------*-C++-*----------------------------------//
// shmem_sr.cc
// Geoffrey M. Furnish
// Fri Mar  6 10:32:46 1998
//---------------------------------------------------------------------------//
// @> Implement SHMEM send/receive functionality.
//---------------------------------------------------------------------------//

#include "config.hh"
#include "global.hh"
#include "C4_Req.hh"

#include "ds++/Assert.hh"

#include <iostream>
using namespace std;

#include <time.h>
#include <mpp/shmem.h>
#include <pthread.h>

// #define SHMDBG 1

C4_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// C4 SHMEM DATA
//
// The C4 messaging layer over SHMEM uses a variety of data structures to
// orchestrate reliable point to point messaging.  This section introduces
// the data definitions used by the various messaging functions in this
// file. 
//---------------------------------------------------------------------------//

int C4_shm_mynode;		// id of this node.
int C4_shm_nodes;		// # of nodes in the virtual multicomputer.

// pe_msg_waiting[node] indicates if there is an inbound message sitting in
// the inbox, waiting to be recived from node.

static int *pe_msg_waiting;

// pe_ready[node] indicates whether node is ready to receive a message.

static int *pe_ready;

// C4_max_buf_sz is the size of the largest message which can be stuffed into 
// the inbox.

// const int C4_max_buf_sz = 1024;
const int C4_max_buf_sz = 4*1024;

// pe_recv_buf_s is the struct which represents the inbox.

struct pe_recv_buf_s {
    int length;
    int tag;
    char data[ C4_max_buf_sz ];
};

pthread_mutex_t *send_lock;
pthread_mutex_t *recv_lock;

// pe_recv_buf[node] is the inbox for receiving messages from node.

static pe_recv_buf_s *pe_recv_buf;

// msg_queue is where we put messages that come in via the inbox, but aren't
// ready to be received yet.

static Msg_Queue *msg_queue;

// We sometimes need to pass the address of a 1 to shmem in order to set some 
// flag on the destination node.  

int one = 1;

static int c4verbose = 0;

//---------------------------------------------------------------------------//

// defined in shmem_reduce.cc
void C4_shm_init_scalar_work_arrays();

static void C4_shm_init_pt2pt();

//---------------------------------------------------------------------------//
// Initialize the data structures needed by the comm layer.
//---------------------------------------------------------------------------//

void C4_shm_init_pt2pt()
{
    C4_shm_mynode = _my_pe();
    C4_shm_nodes = _num_pes();

    int nodes = C4_shm_nodes;

    pe_msg_waiting = (int *) shmalloc( nodes * sizeof(int) );
    Assert( pe_msg_waiting );

    pe_ready = (int *) shmalloc( nodes * sizeof(int) );
    Assert( pe_ready );

    pe_recv_buf = (pe_recv_buf_s *)
	shmalloc( nodes * sizeof(struct pe_recv_buf_s) );
    Assert( pe_recv_buf );

// Initialize the flag arrays.

    for( int i=0; i < nodes; i++ ) {
	pe_msg_waiting[i] = 0;
	pe_ready[i] = 1;
    }

// Now construct the message queues so we have somewhere to store the
// inbound messages.

    msg_queue = new Msg_Queue[ C4_shm_nodes ];

// Build the thread synchronization support for point to point messaging.

    send_lock = new pthread_mutex_t[ C4_shm_nodes ];
    recv_lock = new pthread_mutex_t[ C4_shm_nodes ];

    for( int i=0; i < nodes; i++ ) {
	Insist( !pthread_mutex_init( &send_lock[i], NULL ), 
		"Unable to initialize mutex" );
	Insist( !pthread_mutex_init( &recv_lock[i], NULL ), 
		"Unable to initialize mutex" );
    }

    if (c4verbose)
        cout << "pt2pt arrays are ready on node " << node() << endl;
    gsync();
}


//---------------------------------------------------------------------------//
// Utility functions, not part of the public interface.

//---------------------------------------------------------------------------//
// Extract data from the inbox and store into a buffer.
//---------------------------------------------------------------------------//

static void pull_msg_from_inbox( int source, void *buf, int msglen )
{
    if ( msglen > C4_max_buf_sz ) {
	throw( "receive of packetized inbounds not implemented." );
    } else {
	memcpy( buf, pe_recv_buf[source].data, msglen );
    }
}

//---------------------------------------------------------------------------//
// Update book keeping data to reflect that the inbox is clear and source is
// free to send another message.
//---------------------------------------------------------------------------//

void mark_inbox_as_clear( int source )
{
// Clear my own msg_waiting flag, and tell soruce I'm ready to receive again.

// Must do it in this order to avoid a race.

    pe_msg_waiting[source] = 0;

    shmem_int_put( &pe_ready[C4_shm_mynode], &one, 1, source );
}

//---------------------------------------------------------------------------//
// Orchestrate the receive side of message delivery.
//---------------------------------------------------------------------------//

void pull_msg_from_node( int source, void *buf, int msglen )
{
// Wait till there is a message coming through the inbox.
    while( !pe_msg_waiting[source] )
	;

    char *cbuf = static_cast<char *>( buf );
    int len = pe_recv_buf[source].length;
    Insist( len <= msglen, "Receive buffer not big enough." );

    int base=0, recvlen=0;

    do {
    // Update base from prior packet receive.
	base += recvlen;

    // Calculate length of next packet.
	recvlen = (base + C4_max_buf_sz < len ?
		   C4_max_buf_sz : len - base);

    // Wait for data to be ready.
	while( !pe_msg_waiting[source] )
	    ;

    // Receive the next packet's worth of data.
	memcpy( cbuf+base, pe_recv_buf[source].data, recvlen );

    // Indicate readiness to receive next packet.
	mark_inbox_as_clear( source );
    } 
    while( base + recvlen < len );
    
}

//---------------------------------------------------------------------------//
// Orchestrate the send side of message delivery.
//---------------------------------------------------------------------------//

void push_msg_to_node( int dest, void *buf, int msglen, int tag )
{
// Build header and send.

    struct {
	int length;
	int tag;
    } hdr;

    hdr.length = msglen; hdr.tag = tag;

    shmem_int_put( (int *) &pe_recv_buf[C4_shm_mynode],
		   (int *) &hdr, 2, dest );

    int base=0, xmitlen=0;
    char *cbuf = static_cast<char *>( buf );

// Now send the messages through one packet at a time.
    do {
    // Update base from prior packet send.
	base += xmitlen;

    // Caclulate size of next packet.
	xmitlen = (base + C4_max_buf_sz < msglen ?
		   C4_max_buf_sz : msglen - base);

    // Wait for dest to be ready for next packet
	while( !pe_ready[dest] )
	    ;

    // Send the next packet's worth of data.
	shmem_putmem( (void *) &pe_recv_buf[C4_shm_mynode].data,
		      cbuf+base, xmitlen, dest );

    // indicate dest is not ready
	pe_ready[dest] = 0;

    // set msg_waiting flag in dest.
	shmem_int_put( &pe_msg_waiting[C4_shm_mynode],
		       &one, 1, dest );
    }
    while( base + xmitlen < msglen );
}

//---------------------------------------------------------------------------//
// Pull a message out of the inbox and stuff it onto the message queue.
//---------------------------------------------------------------------------//

static void empty_inbox_to_msg_queue( int source )
{
    if ( !pe_msg_waiting[source] ) return;

    int msgtag = pe_recv_buf[source].tag;
    int msglen = pe_recv_buf[source].length;

// Build the message to be put on our local queue.

    Msg_DB *pm = msg_queue[source].new_msg( msglen );
    pm->tag = msgtag;

    pull_msg_from_inbox( source, pm->buf, msglen );

// Okay, msg is now in pm, so clear the inbox.

    mark_inbox_as_clear( source );

// Now queue the message.

    msg_queue[source].enqueue(pm);
}


//---------------------------------------------------------------------------//
// The Init() function initializes virtual nodes, and the data structures
// used by the C4 SHMEM communication layer.
//---------------------------------------------------------------------------//

void Init( int& iargc, char **& iargv )
{
    int npes = 0;

// Search through the arglist, looking for C4 options.

    int argc = iargc;
    char **argv = iargv;

    argc--, argv++;		// Skip program name.

    while( argc )
    {
	string opt = *argv;

	if (opt == "-npes")
	{
	    Insist( argc > 1, "Syntax: -npes #" );
	    npes = atoi( argv[1] );
	    argc -= 2, argv += 2;
	    continue;
	}

        if (opt == "c4v")
        {
            c4verbose = 1;
            argc--, argv++;
            continue;
        }

	cout << "Unrecognized option: " << *argv << endl;
	argc--, argv++;
    }

    if (c4verbose)
        cout << "Starting C4\n";

// Now fire off the pe's.
    start_pes( npes );

    if (c4verbose)
        cout << "Virtual processor is alive.\n";

    C4_shm_init_scalar_work_arrays();
    C4_shm_init_pt2pt();

// Why is this not the default?
    shmem_set_cache_inv();
}

void Finalize()
{
// Whaddya think?  Should I shmem_free the symmetric data buffer?  heh heh.

// Check that the shmalloc books haven't been hosed.

#ifdef _CRAYMPP
    Assert( !shmalloc_check(0) );
#endif

    shfree( pe_msg_waiting );
    shfree( pe_ready );
    shfree( pe_recv_buf );

    delete[] msg_queue;

// Check that the shmalloc books haven't been hosed.  Again, just to be sure.

#ifdef _CRAYMPP
    Assert( !shmalloc_check(0) );
#endif
}

int node()
{
    return C4_shm_mynode;
}

int nodes()
{
    return C4_shm_nodes;
}

int group()
{
    int group = 0;
    return group;
}

void gsync()
{
    shmem_barrier_all();
}

struct timespec tsclock;

double Wtime()
{
    int clock =
#if defined(CLOCK_SGI_CYCLE)
        CLOCK_SGI_CYCLE;
#else
    CLOCK_REALTIME;
#endif

    clock_gettime( clock, &tsclock );
    return tsclock.tv_sec + tsclock.tv_nsec*1.0e-9;
}

double Wtick() 
{
    int clock =
#if defined(CLOCK_SGI_CYCLE)
        CLOCK_SGI_CYCLE;
#else
    CLOCK_REALTIME;
#endif

    clock_getres( clock, &tsclock );
    return tsclock.tv_sec + tsclock.tv_nsec*1.0e-9;
}

//---------------------------------------------------------------------------//
// MPI send/receive calls (basic set)
//
// Synchronous:
//	MPI_Send(void* buf, int count, MPI_Datatype datatype,
//		 int dest, int tag, MPI_Comm comm);
//	MPI_Recv(void* buf, int count, MPI_Datatype datatype,
//		 int source, int tag, MPI_Comm comm, MPI_Status *status);
//
// Asynchronous:
//	MPI_Isend(void* buf, int count, MPI_Datatype datatype,
//		  int dest, int tag, MPI_Comm comm, MPI_Request *request);
//	MPI_Irecv(void* buf, int count, MPI_Datatype datatype,
//		  int source, int tag, MPI_Comm comm, MPI_Request *request);
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Perform a normal (blocking) send.
//---------------------------------------------------------------------------//

int SHM_Send( void *buf, int size, int dest, int tag, int group )
{
// Acquire access lock to the send data structures.
    Insist( !pthread_mutex_lock( &send_lock[dest] ),
	    "Error aequiring send lock." );

// Now busy wait till dest is ready.

    while( !pe_ready[dest] )
    // Do not have to flush cache, b/c C4_Init sets it for messages to
    // automatically invalidate the cache.
	;

    push_msg_to_node( dest, buf, size, tag );

// Release the send lock so other threads can send data too.
    Insist( !pthread_mutex_unlock( &send_lock[dest] ),
	    "Error releasing send lock." );

    return C4_SUCCESS;
}

//---------------------------------------------------------------------------//
// Perform a normal (blocking) receive.
//---------------------------------------------------------------------------//

int SHM_Recv( void *buf, int size, int source, int tag, int group )
{
    int msgtag, msglen = 0;

    if (source == C4_Any_Source)
    {
	throw "bogus";
    }
    else
    {
    // Receive only from specified node 'source'.

	Insist( !pthread_mutex_lock(&recv_lock[source]),
		"Erorr acquiring recv_lock." );

	int bytes_recvd = 0;

    // This is the "easy" case.  Look for a message from source which matches
    // tag.

	Insist( source >= 0 && source < C4_shm_nodes,
		"C4_Recv: source node out of range." );

#ifdef SHMDBG
	printf( "%d trying to receive a msg.\n", C4_shm_mynode );
#endif

    // First, check to see if the message queue is empty.  If it's not, the
    // message we want could be in there.

	if ( !msg_queue[source].is_empty() ) {

	// search msg queue. If find matching tag, copy to user buffer, and
	// remember to empty the inbox before returning.

	    bool match = false;
	    for( int i=0; i < msg_queue[source].num_msgs(); i++ ) {

		Msg_DB *pm = msg_queue[source].msg(i);

		if ( pm->matches(tag) ) {
		// Dequeue this message, copy to user buffer.

		    if (pm->msg_size > size) {
			throw( "recv buffer not large enough for message." );
		    }

		    bytes_recvd = pm->msg_size;
		    memcpy( buf, pm->buf, bytes_recvd );

		    msg_queue[source].discard(i);

		    match = true;
		    break;
		}
	    }

	    if (match) {
		Insist( !pthread_mutex_unlock(&recv_lock[source]),
			"Unable to release recv_lock." );
		return bytes_recvd;
	    }
	}

    // It wasn't already in the queue.  So, ... let's just start pulling
    // messages in from the inbox.  If we find what we are looking for,
    // return it.  If not, queue it.

#ifdef SHMDBG
	printf( "%d, it wasn't in the queue, looking in inbox.\n", C4_shm_mynode );
#endif

	bool match = false;
	do {
	// Now busy wait till one shows up for me.

	    while( !pe_msg_waiting[source] )
		;

#ifdef SHMDBG
	    printf( "%d now there is a message waiting.\n", C4_shm_mynode );
#endif

	// Okay, something is finally here, now dispense.

	    msgtag = pe_recv_buf[source].tag;
	    msglen = pe_recv_buf[source].length;

#ifdef SHMDBG
	    printf( "%d, tag=%d, msgtag=%d, msglen=%d\n",
		    C4_shm_mynode, tag, msgtag, msglen );
#endif

	    if (tag == C4_Any_Tag || tag == msgtag) {
#ifdef SHMDBG
		printf( "%d, match!\n", C4_shm_mynode );
#endif

	    // Receive and exit this spinlock.
		pull_msg_from_node( source, buf, size );

		Insist( !pthread_mutex_unlock(&recv_lock[source]),
			"Unable to release recv_lock." );
		return msglen;
	    }
	    else {
		empty_inbox_to_msg_queue( source );

		msglen = 0;	// clear this, just to be safe.
	    }

	} while( !match );
    }

    throw "SHM_Recv: fell off the end!";
}

struct async_msg_stuff {
    void *buf;
    int size;
    int node;
    int tag;
    int group;
  public:
    async_msg_stuff( void *b, int s, int n, int t, int g )
	: buf(b), size(s), node(n), tag(t), group(g)
    {}
};

void *SHM_thread_send_async( void *arg )
{
    async_msg_stuff *mdb = (async_msg_stuff *) arg;

    SHM_Send( mdb->buf, mdb->size, mdb->node, mdb->tag, mdb->group );

    delete mdb;
    return NULL;
}

//---------------------------------------------------------------------------//
// Perform a non blocking send.
//---------------------------------------------------------------------------//

C4_Req SHM_SendAsync( void *buf, int size, int dest, int tag, int group )
{
    async_msg_stuff *mdb = new async_msg_stuff( buf, size, dest,
						tag, group );
    C4_Req r;

    if (pthread_create( &r.thread, NULL, SHM_thread_send_async, mdb ))
	printf( "%d unable to spawn send thread.", C4_shm_mynode );

    r.type = Async_Send;
    r.set();

    return r;
}

void *SHM_thread_recv_async( void *arg )
{
    async_msg_stuff *mdb = (async_msg_stuff *) arg;

    SHM_Recv( mdb->buf, mdb->size, mdb->node, mdb->tag, mdb->group );

    delete mdb;
    return NULL;
}

//---------------------------------------------------------------------------//
// Perform a non blocking receive.
//---------------------------------------------------------------------------//

C4_Req SHM_RecvAsync( void *buf, int size, int source, int tag, int group )
{
    async_msg_stuff *mdb = new async_msg_stuff( buf, size, source,
						tag, group );
    C4_Req r;

    if (pthread_create( &r.thread, NULL, SHM_thread_recv_async, mdb ))
	printf( "%d unable to spawn recv thread.\n", C4_shm_mynode );

    r.type = Async_Recv;
    r.set();

    return r;
}

//---------------------------------------------------------------------------//
// "Optimized" forms of the above.  These are supposed to avoid spurios
// object creation and copy.  However, for the purpose of getting going,
// these are implemented by just calling the ones above.  Later we can switch
// the implementations so that these go faster.
//---------------------------------------------------------------------------//

void SHM_SendAsync( C4_Req& r, void *buf, int size, int dest, int tag, 
		    int group /*=0*/ )
{
// Not checking that r is not in use, which is of course a concern...
    r = SHM_SendAsync( buf, size, dest, tag, group );
}

void SHM_RecvAsync( C4_Req& r, void *buf, int size, int source, int tag, 
		    int group /*=0*/ )
{
// Not checking that r is not in use, which is of course a concern...
    r = SHM_RecvAsync( buf, size, source, tag, group );
}

//---------------------------------------------------------------------------//
// Wait on an asynchronous message to complete.
//---------------------------------------------------------------------------//

void C4_Wait( pthread_t& tid )
{
    pthread_join( tid, NULL );
}

C4_NAMESPACE_END

#include "ds++/DynArray.cc"

template class dsxx::DynArray<C4::Msg_DB *>;

//---------------------------------------------------------------------------//
//                              end of shmem_sr.cc
//---------------------------------------------------------------------------//
