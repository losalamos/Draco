//----------------------------------*-C++-*----------------------------------//
// Banded_Matrix.cc
// Geoffrey M. Furnish
// Thu Mar 12 12:51:51 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "linalg/Banded_Matrix.hh"

#include "c4/SpinLock.hh"

#include <algorithm>
#include <functional>
#include <list>

//---------------------------------------------------------------------------//
// Compute and store the info which describes the messaging pattern required
// to perform mathematically correct matrix-vector multiplies.
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::compute_message_routing()
{
// Build up an "occupancy vector" which shows which column positions are
// filled in any row owned by this processor.

    Mat1<int> occupancy( nrt );
    occupancy = 0;

// Now figure out which columns we use.

    for( int n=0; n < N; n++ )
    {
    // Look at diagonal number n.

	for( int c=0; c < nrp; c++ ) {
	    int col = nro + diag_offset[n] + c;
	    if ( 0 <= col && col < nrt )
		occupancy( col ) = 1;
	}
    }

    if (verbose) emit_occupancy_vectors( occupancy );

// Now count up the occupied cells so we know how big the data buffer needs
// to be.

    int ncells = std::count( occupancy.begin(), occupancy.end(), 1 );

    SPINLOCK( cout << node() << " uses " << ncells << " cells." << endl );

// Now allocate the space for receiving data.

    cvdata.redim( ncells );

// Set the indexes to an invalid value, so that any screwups will be caught.

    cvindex = -1;

// Calculate the local index (in cvdata) where each occupied element will
// reside.

    Mat1<int> local_index( nrt ), global_index( ncells );
    global_index = -1;

    int current_index = 0;
    for( int i=0; i < nrt; i++ )
	if (occupancy(i) == 1) {
	    local_index(i) = current_index;
	    global_index(current_index) = i;
	    current_index++;
	}

    Check( current_index == ncells );

// Need to figure out which processor owns each of these elements.  I guess
// we should start by figuring out which processors own which rows.

    Mat1<int> start_row( nodes() ), nrows( nodes() );

    start_row = 0;
    nrows = 0;

    start_row( node() ) = nro;
    nrows( node() ) = nrp;

    gsum( start_row.begin(), nodes() );
    gsum( nrows.begin(), nodes() );

    Mat1<int> cell_owner( ncells );
    Mat1<int> active_cells_owned_by( nodes() );
    active_cells_owned_by = 0;

    for( int c=0; c < ncells; c++ )
    {
	int gi = global_index(c);
	int p=0;
	for( ; p < nodes(); p++ )
	    if (start_row(p) <= gi && gi < start_row(p) + nrows(p))
		break;
	Check( p < nodes() );
	cell_owner( c ) = p;
	active_cells_owned_by(p)++;
    }

    {
	HTSyncSpinLock h;

	cout << node() << " active_cells_owned_by: ";
	for( int i=0; i < nodes(); i++ )
	    cout << active_cells_owned_by(i) << "  ";
	cout << endl;
    }

// Calculate nsenders.

    senders = std::count_if( active_cells_owned_by.begin(),
			     active_cells_owned_by.end(),
			     std::bind2nd( std::greater<int>(), 0 ) );

    SPINLOCK( cout << node() << " has " << senders << " senders.\n" );

// Resize some things that require one element per sender.

    sender_node.redim( senders );
    sender_nels.redim( senders );
    sender_head.redim( senders );

    rreq.redim( senders );

// Now walk through the senders, and record each sender's node id, the number 
// of elements they'll be sending us, and the location in the cvdata buffer
// where that sender's first element will go.

    int s=0;
    for( int p=0; p < nodes(); p++ )
	if (active_cells_owned_by(p) > 0)
	{
	    sender_node(s) = p;
	    sender_nels(s) = active_cells_owned_by(p);
	    sender_head(s) = ( s == 0 ?
			       0 :
			       sender_head(s-1) + sender_nels(s-1) );
	    s++;
	}

    Check( s == senders );
    Check( senders == 0 ||
	   sender_head(senders-1) + sender_nels(senders-1) == ncells );

// Now check to see if we are sending data to ourself.  If so, annotate this
// correctly so that the send/receive code can use copies instead of
// messaging for sending data to self.

    self_sender_id = -1;

    for( int s=0; s < senders; s++ )
	if (sender_node(s) == node()) self_sender_id = s;

    cout << node() << " set self_sender_id\n" << flush;

// Okay, the above basically represents the determination of who is sending
// us data.  Now we need the opposite, who we are sending data to, or who
// receives data from us.

// We're gonna need some buffer space of size nrp.  But nrp isn't always the
// same on all nodes, so fetch the max.

    int max_nrp = nrp;
    gmax( max_nrp );

    Mat1<int> global_index_buf( max_nrp );

    receivers = 0;
    int max_noc = 0;

    std::list<rcvr_tag *> rcvr_list;

// Now loop through all the nodes.  Each one in trun listens to the
// requirements articulated by the others.

    for( int p=0; p < nodes(); p++ )
    {
	gsync();

	if (p == node())
	{
	// Okay, the other nodes are going to tell us what they need from
	// us.
	    cout << p << " ready to listen." << endl;

	    for( int n=0; n < nodes(); n++ )
	    {
		if (n == node()) {
		// Fill up struct for ourself.

		    int noc = std::count( &occupancy(start_row(p)),
					  &occupancy(start_row(p)) + nrows(p),
					  1 );

		    cout << p << " needs to send " << noc << " to self." << endl;

		    Check ( 0 <= noc && noc <= max_nrp );

		    if (noc > 0)
		    {
			receivers++;
			max_noc = std::max( noc, max_noc );
			
		    // Build the list of indices we need from ourself.

			int ndx = 0;
			for( int i = nro; i < nro + nrp; i++ )
			    if (occupancy(i))
				global_index_buf(ndx++) = i;

			Check( ndx == noc );

			cout << node() << " pushing new rcvr_tag( " << n
			     << ", " << noc << ", "
			     << global_index_buf(0) << "  "
			     << global_index_buf(1) << "  " << endl;

			rcvr_list.push_back(
			    new rcvr_tag( n, noc, global_index_buf.begin() )
			    );
		    }
		}
		else {
		// Tell the next processor that we are ready to hear his
		// requirements.  This may seem kind of dumb, but experience
		// has shown it helps a lot to use handshaking in situations
		// like this.

		    int dummy=0;
		    Send( dummy, n );

		    int noc;
		    Recv( noc, n );

		    Check( noc >= 0 );

		    cout << p << " needs to send " << noc << " to " << n << endl;

		    if (noc > 0)
		    {
			receivers++;
			max_noc = std::max( noc, max_noc );

		    // He needs data from us,  Now he's gonna tell us which
		    // elements exactly it is that he needs from us.

			int nn = Recv( global_index_buf.begin(), noc, n );
			Check( nn == noc );

			rcvr_list.push_back(
			    new rcvr_tag( n, noc, global_index_buf.begin() )
			    );

			for( int k=0; k < nn; k++ )
			    cout << node() << " needs to send index "
				 << global_index_buf(k) << " to node "
				 << n << endl;

		    }
		}
	    }
	}
	else {
	    if (p == node()) continue; // Nothing to do here.

	// Wait for processor p to say he's listening to us.  This seems
	// dumb, but prior experience with message passing implementations
	// indicates that we can easily swamp the inbox if we don't do a
	// little handshaking in situations like this.

	    int dummy;
	    Recv( dummy, p );

	// Tell processor p how many elements we need to get from him.

	    int noc = std::count( &occupancy(start_row(p)),
				  &occupancy(start_row(p)) + nrows(p), 1 );

	    Check ( 0 <= noc && noc <= max_nrp );

	    Send( noc, p );

	    if (noc > 0)
	    {
	    // Build the list of indices we need from him ...

		int n = 0;
		for( int i = start_row(p);
		     i < start_row(p) + nrows(p); i++ )
		    if (occupancy(i))
			global_index_buf(n++) = i;

		Check( n == noc );

		cout << node() << " telling " << p << " that we need "
		     << noc << " els with indices "
		     << global_index_buf(0) << "  "
		     << global_index_buf(1) << "  " << endl;
// 		     << global_index_buf(2) << "  "
// 		     << global_index_buf(3) << "  " << endl;

	    // and send it over.

		Send( global_index_buf.begin(), noc, p );
	    }
	}
    }

    if (verbose)
	SPINLOCK( cout << node() << " has " << receivers
		  << " receivers, with max_noc=" << max_noc << endl );

// Okay, by this point, rcvr_list has the prescription of what data has to
// get shipped where.  Now we just need to build up the receiver structs
// based on this info.

    receiver_node.redim( receivers );
    receiver_nels.redim( receivers );
    receiver_ndxs = Mat2<int>( max_noc, receivers );
    receiver_data = Mat2<T>( max_noc, receivers );
    sreq.redim( receivers );

    int rcvr = 0;
    for( std::list<rcvr_tag *>::iterator x = rcvr_list.begin();
	 x != rcvr_list.end(); x++ )
    {
	cout << node() << " rcvr=" << rcvr << " x: node=" << (*x)->node
	     << " nels=" << (*x)->nels << "  "
	     << (*x)->indexes[0] << "  "
	     << (*x)->indexes[1] << "  "
	     << endl;

    // Stuff the data for this recevier into our receiver structs.
	receiver_node(rcvr) = (*x)->node;
	receiver_nels(rcvr) = (*x)->nels;
	for( int i=0; i < receiver_nels(rcvr); i++ )
	    receiver_ndxs(i,rcvr) = (*x)->indexes[i];

	delete *x;		// We're done, wipe it away.

	rcvr++;
    }

    if (verbose) emit_receiver_structs();

// Finally, check to see if we are receiving data from ourself.  If so,
// annotate this correctly so that the send/receive code can use copies
// instead of messaging for sending data to self.

    self_receiver_id = -1;

    for( int r=0; r < receivers; r++ )
	if (receiver_node(r) == node()) self_receiver_id = r;

    cout << node() << " set self_receiver_id\n" << flush;

// Now we need to work on setting up the start_bands, nbands, and cvindex
// arrays. 

    nbands = 0;

// Loop over rows.

    for( int r=0; r < nrp; r++ )
    {
	int sb, eb;

    // Find the starting band, that is, the first diagnoal which is on
    // processor for this row.

	for( sb=0; sb < N; sb++ )
	{
	    int c = nro + diag_offset[sb] + r;
	    if (0 <= c) {
		start_band(r) = sb;
		break;
	    }
	}

    // Now count the bands which are on processor.

	for( eb=sb; eb < N; eb++ )
	{
	    int c = nro + diag_offset[eb] + r;
	    if (c >= nrt)
		break;
	}

	nbands(r) = eb - sb;

	cout << node() << " row=" << nro+r
	     << " start_band=" << start_band(r)
	     << " nbands=" << nbands(r) << endl;
    }

// Okay, we have the active bands per row calculated.  Now we should be able
// to calculate the indexes into the cvdata buffer, which will be needed by
// each row.

    for( int r=0; r < nrp; r++ )
	for( int ib=0; ib < nbands(r); ib++ ) {
	    int b = start_band(r) + ib;
	    int c = nro + diag_offset[b] + r;
	    cvindex(b,r) = local_index(c);
	}

// Let's print it out, just for fun.

    {
	cout << flush;
	HTSyncSpinLock h;

	cout << node() << " dumping cvindex:" << endl;
	for( int r=0; r < nrp; r++ ) {
	    for( int b=0; b < N; b++ )
		cout << cvindex(b,r) << " ";
	    cout << endl;
	}
    }

#if 0
// Now we need to deduce which processors we will need to send elements from
// our column vector to, as well as the inverse--which processors will be
// sending us data from their column vector.  

// Let's start by getting the occupancy vector for each node communicated to
// each other node.

// This is a horribly inefficient way to calculate this, but we have to start 
// somewhere...  REPLACE THIS LATER!!!

    Mat2<int> ov( nrt, nodes() );
    ov = 0;
    {
	Mat1<C4_Req> rreq( nodes() );
	for( int i=0; i < nodes(); i++ ) {
	    if (i == node()) continue;
	    RecvAsync( rreq[i], &ov(0,i), nrt, i );
	}
	for( int i=0; i < nodes(); i++ ) {
	    if (i == node()) continue;
	    Send( &occupancy(0), nrt, i );
	}
    // Now copy ours into place.
	for( int i=0; i < nrt; i++ )
	    ov( i, node() ) = occupancy(i);
    }

    {
	HTSyncSpinLock h;
	
	cout << "Ocupancy vectors for node " << node() << endl;

	for( int n=0; n < nodes(); n++ ) {
	    for( int i=0; i < nrt; i++ )
		cout << ov(i,n) << " ";
	    cout << endl;
	}
    }

    gsync();

// Okay, we have the occupancy vectors for all nodes now, so we should be
// able to start calculating the required messaging patterns.

// Consider our portion of the column vector.  Each of our elements could
// have to be sent to as few as zero or as many as all of the other
// processors. 

    Mat1<bool> recipient( nodes() );
    recipient = false;

    SPINLOCK( cout << "Ready to calculate recipient list.\n" );

    for( int n=0; n < nodes(); n++ )
	for( int i=nro; i < nro+nrp; i++ )
	    if ( ov(i,n) )
		recipient(n) = true;

    int nrecipients = std::count( recipient.begin(), recipient.end(), true );

    SPINLOCK( cout << node() << " has " << nrecipients << " recipients.\n" );

    {
	HTSyncSpinLock h;

	cout << node() << " : ";
	for( int i=0; i < nodes(); i++ )
	    cout << recipient(i) << " ";
	cout << endl;
	cout << flush;
    }

    Check( nrecipients <= nodes() );

    gsync();

// Now we know how many recipients there are, so we can allocate the send
// buffers for them.  For simplicity sake, we will make the buffers be as big 
// as our share of elements is.  This should only be a small waste, since we
// are a distributed matrix anyway.

    receivers = nrecipients;
    receiver_node.redim( receivers );
    receiver_nels.redim( receivers );
//    receiver_ndxs = Mat2<int>( nrp, receivers );
//    receiver_data = Mat2<T>( nrp, receivers );
    receiver_ndxs = new int*[ receivers ];
    receiver_data = new T*[ receivers ];

    sreq.redim( receivers );

    for( int n=0,		// node # of i-th receiver node.
	     i=0;
	 i < receivers; i++ )
    {
    // Figure out which is the next receiver.
	while( !recipient(n) ) n++;

    // Keep data on this node.
	receiver_node(i) = n;
	receiver_nels(i) = 0;

	for( int j=0; j < nrows(n); j++ )
	    if ( ov(start_row(n)+j,n) )
		receiver_nels(i)++;

	receiver_ndxs[i] = new int[ receiver_nels(i) ];
	receiver_data[i] = new T[ receiver_nels(i) ];

	int nels=0;
	for( int j=0; j < nrows(n); j++ )
	    if ( ov(start_row(n)+j,n) )
		receiver_ndxs[i][nels++] = j;

	Check( nels == receiver_nels(i) );

    // Now point the index to the next node for consideration.
	n++;
    }

    cout << node() << " done building receiver structs." << endl;
#endif
//     if (verbose) emit_receiver_structs();

// Okay, now we need the inverse map.  In the above we calculated what others 
// need to receive from us.  Now we need to calculate what others will be
// sending us.

// I guess we should start by determining which processor owns each row.
// Then we can map our occupancy vector against this to see which processors
// will be sending us data.
#if 0
    Mat1<int> owner( nrt );
    owner = 0;
    for( int i=0; i < nrp; i++ )
	owner( nro+i ) = node();

    C4::gsum( &owner(0), nrt );

// Now we can look at our occupancy vector, and determine who will be sending
// us data.

    int nsenders = 0,
	last_discovered_sender_node = -1;

    for( int i=0; i < nrt; i++ )
	if (occupancy(i) == 1) {
	    if (owner(i) != last_discovered_sender_node) {
		last_discovered_sender_node = owner(i);
		nsenders++;
	    }
	}

// Okay we now know how many nodes will be sending us data, so we can make
// places to receive the data now.

    senders = nsenders;
    sender_ndxs = Mat2<int>( nrp, senders );
    sender_data = Mat2<T>( nrp, senders );

    sender_nels = 0;

// Now start building up these data structures.

    last_discovered_sender_node = -1;
    int nsender = -1;
    for( int i=0; i < nrt; i++ )
	if (occupancy(i)) {
	    if (owner(i) != last_discovered_sender_node) {
		last_discovered_sender_node = owner(i);
		nsender++;
	    }
	    sender_node(nsender) = owner(i);
	    sender_ndxs( sender_nels(nsender)++, nsender ) = i;
	}

    if (verbose) emit_sender_structs();
#endif
}

//---------------------------------------------------------------------------//
// Package up our data and send it off to the processors that need it in
// order to do their part of the work.
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::initiate_sends( const T *pd )
{
    SPINLOCK( cout << node() << " initiating send's." << endl );

    for( int r=0; r < receivers; r++ )
    {
    // Pack data into send buffers.

    // Note that receiver_ndxs has /global/ indexes, but pd is the local data 
    // vector, so we have to convert the global indexes to local indexes.

	for( int i=0; i < receiver_nels(r); i++ )
	    receiver_data(i,r) = pd[ receiver_ndxs(i,r) - nro ];

    // But don't send message to ourselves.

	if (receiver_node(r) == node())
	    continue;		// Copy data when completing receives.

    // Post async send to other processor.
	SendAsync( sreq[r],
		   &receiver_data(0,r), receiver_nels(r),
		   receiver_node(r) );
    }
}

//---------------------------------------------------------------------------//
// Make sure the async sends are all done, in preparation for reusing the
// send buffers.
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::complete_sends()
{
    SPINLOCK( cout << node() << " completing send's." << endl );

    for( int r=0; r < receivers; r++ )
	sreq[r].wait();

}

//---------------------------------------------------------------------------//
// Post receives for the data which is inbound to us.
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::initiate_receives()
{
    SPINLOCK( cout << node() << " initiating recv's." << endl );

    for( int s=0; s < senders; s++ )
    {
	if (sender_node(s) == node()) continue;

	cout << node() << " posting recv for sender " << s
	     << " s_node=" << sender_node(s) << endl;
	cout << "s_head=" << sender_head(s)
	     << " s_nels=" << sender_nels(s) << endl;

	RecvAsync( rreq(s),
		   &cvdata( sender_head(s) ), sender_nels(s),
		   sender_node(s) );
    }
}

//---------------------------------------------------------------------------//
// Make sure all data has been received and/or moved into place.
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::complete_receives()
{
    SPINLOCK( cout << node() << " completing recv's." << endl );

    for( int s=0; s < senders; s++ )
	if (sender_node(s) == node())
	{
	    cout << node() << " copying own data into place." << endl;

	    Check( self_receiver_id >= 0 && self_receiver_id < receivers );
	    Check( self_sender_id >= 0 && self_sender_id < senders );

	// Copy our own data into place.
	    for( int i=0; i < receiver_nels( self_receiver_id ); i++ )
		cvdata( sender_head( self_sender_id ) + i ) =
// 		    receiver_data[ self_receiver_id ][i];
		    receiver_data(i,self_receiver_id);

	    cout << node() << " done copying own data into place." << endl;
	}
	else
	    rreq[s].wait();
}

//---------------------------------------------------------------------------//
// Ta Da!  Perform matrix vector multiply...
//---------------------------------------------------------------------------//

template<class T, int N> template<class ColumnVector>
void Banded_Matrix<T,N>::multiply( const ColumnVector& x,
				   ColumnVector& b )
{
#if 0
    cout << node() << " doing the matvec thing on a ColumnVector of type\n"
	 << typeid(ColumnVector).name()
	 << endl;

    Mat1<C4_Req> rreq(senders);

// Start by posting receives for the data.
    for( int i=0; i < senders; i++ ) {

    // Skip ourselves, we will just copy the data into place.
	if (sender_node(i) == node()) continue;

    // Post receives for data coming from other nodes.
	RecvAsync( rreq(i), &sender_data(0,i), sender_nels(i),
		   sender_node(i) );
    }

// Now send our data to wherever it's going.
    for( int i=0; i < receivers; i++ ) {

    // Don't use messaging to send data to self.
	if (receiver_node(i) == node()) {
	// Just copy into place.
	    for( int j=0; j < receiver_nels(i); j++ )
		sender_data(j,i) = x( receiver_ndxs(j,i) );

	// Go on to next receiver.
	    continue;
	}

    // Copy the data we're sending into the buffers for transmission to
    // receivers. 
	for( int j=0; j < receiver_nels(i); j++ )
	    receiver_data(j,i) = x( receiver_ndxs(j,i) );

	Send( &receiver_data(0,i), receiver_nels(i), receiver_node(i) );
    }

    cout << node() << " done sending data." << endl;

    for( int i=0; i < senders; i++ )
	rreq(i).wait();

    cout << node() << " done receiving data." << endl;

// Need a data structure which tells me for each row that I won, which
// elements of the bands are to be multiplied agasint which elements of the
// inbound data sent to me by the senders.

    Mat1<int> nels( nrp );
    nels = 0;
    Mat2<int> band(N,nrp), index(N,nrp);
    band = 0;
    index = 0;
    Mat1<T> xdata( &sender_data(0,0), nrp*senders );

    for( int i=0; i < nrp; i++ )
    {
	b(i) = 0;

    // Now multiply the elements against the data which has been sent to us
    // and sum it all up.

	for( int j=0; j < nels(i); j++ )
	    b(i) += data( i, band(j,i) ) * xdata( index(j,i) );

    }
#endif
    complete_sends();		// From the previous call...

    initiate_receives();
    initiate_sends( x.begin() );
    complete_receives();

    gsync();

    if (verbose)
    {
	cout << flush;
	HTSyncSpinLock h;

	cout << "Node " << node() << " dumping sparse data:\n";
	cout << "cvdata: size=" << cvdata.size() << endl;
	for( int i=0; i < cvdata.size(); i++ )
	    cout << cvdata(i) << " ";
	cout << endl;
    }

// Multiply the data.

    for( int r=0; r < nrp; r++ )
    {
    // Initialize the receiver array.
	b(r) = 0;

    // Figure out which bands we're working with.
	int sb = start_band(r), nb = nbands(r);

	cout << node() << " processing row " << r
	     << " sb=" << sb << " nb=" << nb << endl;

    // Now loop over the particpating bands, and accumulate the product.
	for( int ib = sb; ib < sb+nb; ib++ )
	    b(r) += data(ib,r) * cvdata( cvindex(ib,r) );
    }

// You might think we would complete sends here, but there is no reason to
// stall this processor here waiting on other processors to complete their
// receives.  So, we just let 'er rip, and we won't actually block on these
// sends until we have to reuse the send buffers on the next time this method 
// is called.
}

//---------------------------------------------------------------------------//
// For each processor, print out the occupancy vector.  This should look a
// little bit like the fill pattern for the matrix being represented, except
// that diagonals get turned into blocks based on which processors own which
// rows. 
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::emit_occupancy_vectors( const Mat1<int>& ov ) const
{
    cout << flush;
    HTSyncSpinLock h;

    for( int i=0; i < nrp; i++ ) {
	cout << node() << "  " << nro+i << ":  ";
	for( int j=0; j < nrt; j++ )
	    cout << ov(j) << "  ";
	cout << endl;
    }
}

//---------------------------------------------------------------------------//
// Dump the receiver structs for visual inspection.
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::emit_receiver_structs() const
{
    cout << flush;
    HTSyncSpinLock h;

    cout << "\nReceiver structs for node " << node() << endl;
    cout << "receivers = " << receivers << endl;
    for( int r=0; r < receivers; r++ )
    {
	cout << "The " << r << "th receiver is node "
	     << receiver_node(r) << endl;
	cout << "which needs to obtain " << receiver_nels(r)
	     << " with indexes ";
	for( int n=0; n < receiver_nels(r); n++ )
	    cout << receiver_ndxs(n,r) << "  ";
	cout << endl;
    }
}

//---------------------------------------------------------------------------//
// Dump the sender structs for visual inspection.
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::emit_sender_structs() const
{
    cout << flush;
    HTSyncSpinLock h;

    cout << "\nSender structs for node " << node() << endl;
    cout << "senders = " << senders << endl;
    for( int s=0; s < senders; s++ )
    {
	cout << "The " << s << "th sender is node "
	     << sender_node(s) << endl;
	cout << "which sends " << sender_nels(s)
	     << " with indexes ";
	for( int n=0; n < sender_nels(s); n++ )
	    cout << sender_ndxs(n,s) << "  ";
	cout << endl;
    }
}

//---------------------------------------------------------------------------//
//                              end of Banded_Matrix.cc
//---------------------------------------------------------------------------//
