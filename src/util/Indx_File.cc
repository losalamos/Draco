//----------------------------------*-C++-*----------------------------------//
// Indx_File.cc
// Geoffrey Furnish
// 3 March 1994
//---------------------------------------------------------------------------//
// @> A class for managing a file containing indexed entries of unspecified
// @> form and content.
//---------------------------------------------------------------------------//

#include "util/Indx_File.hh"
#include "ds++/Assert.hh"

#if defined(__linux) || defined(_CRAYMPP) || defined(__sgi)
#define MAP_SHARED_IS_BUSTED 1
#endif

#if defined(_CRAYMPP)
#define NO_MMAP_AT_ALL
#endif

#if defined(__sgi) || defined(__sun)
#define MAP_FILE 0
#endif

#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#ifndef NO_MMAP_AT_ALL
extern "C" {
#include <sys/mman.h>
}
#endif

// Check for compilation on braindead SUNOS systems.  grrrrrrrrr.

#ifndef SEEK_SET
#define SEEK_SET 0
#endif

//---------------------------------------------------------------------------//
// Construct a TOC_entry from its member/arguments.
//---------------------------------------------------------------------------//

template<class Key>
TOC_entry<Key>::TOC_entry( const Key& k, int _base, int _len /*=0*/ )
    : key(k), base(_base), len(_len)
{
}

//---------------------------------------------------------------------------//
// Compare two TOC_entry's for equality.
//---------------------------------------------------------------------------//

template<class Key>
int TOC_entry<Key>::operator==( const TOC_entry<Key>& e ) const {
    if ( !(key == e.key) ) return 0;

// Should we check len?  Are we comparing the key's, the table of
// contents entry, or what?  Hmm.  Maybe we should insist they are all
// equall.  Grrr, have to think about this.

    if ( !(len == e.len) ) return 0;

// What about base?  Should it matter where the data is in the file?
// Let's say it doesn't, for now.

    return 1;
}

//---------------------------------------------------------------------------//
// Constructor.  Open the file, do initial setup, etc.
//---------------------------------------------------------------------------//

template<class Key>
Indx_File<Key>::Indx_File( String _fname, int _amode, 
			   int _max_toc_entries /*=0*/,
			   int _active /*=1*/ )
    : fname(_fname), access_mode(_amode),
      max_toc_entries(_max_toc_entries),
      toc(max_toc_entries), active(_active)
{
    if (!active) return;

    int i;

    ce = toc_entries = 0;
    search_dir = 1;
    write_in_progress = 0;

    header_size = get_header_size();

    toc_base = 0 + header_size;

// Now see about openning the file.

    switch(access_mode) {

    case READ:
	fp = fopen( fname, "r" );

	Insist( fp, "Unable to open Indx_File<T> for reading.\n" );

	read_header();

    // Do some sanity checks.

	if ( header.entry_size != get_toc_entry_size() )
	    cerr << "bogus file.\n";

	max_toc_entries = header.max_entries;

    // Read in the index table.

	read_toc();
	toc_entries = header.entries;

	break;

    case WRITE:
	fp = fopen( fname, "w" );
	if (!fp)
	    cerr << "Indx_File<T>: unable to open " << fname
		 << " for writing.\n\n";

    // Write header, indicating size of TOC_Entry, # thereof, base of data
    // area, etc.

	header.entry_size = get_toc_entry_size();
	if (max_toc_entries)
	    header.max_entries = max_toc_entries;
	else {
	    cerr << "max_entries needs to be specified.  Defautling.\n";
	    header.max_entries = max_toc_entries = 20;
	}
	header.entries = 0;
	header.file_length = 0;
	write_header();

    // Could think about leaving a hole, maybe even a big one, to make
    // expansion of the toc easier.

	data_base = toc_base + max_toc_entries * sizeof(TOC_entry<Key>);
	break;

    case RW_MMAP:

    // This is used to allow read/write access to the file via mmap.  Such
    // use does not allow extending the file, but is useful for "rearranging"
    // data in an existing file.  The currently envisioned use for this is to
    // allow particle sorting in particle history files.

	if ( (fd = open( fname, O_RDWR )) < 0 ) {
	    throw( "Can't open file for read/write." );
	}
	fp = fdopen( fd, "rb+" );

	struct stat statbuf;
	fstat( fd, &statbuf );
	mm_size = (int) statbuf.st_size;

#ifdef MAP_SHARED_IS_BUSTED
    // Just one problem.  Linux (perhaps other unices?) does not support the
    // combination of PROT_WRITE and MAP_SHARED.  In our case we don't really
    // need to share with other concurrent processes.  We just want to use
    // MAP_SHARED in order to get any changes we make communicated back to
    // the disk file.  If you use MAP_PRIVATE, the changes don't go back to
    // the file.  So, we fake it.  We do this by reading the file into a
    // memory buffer, handle the access functions normally, and then be sure
    // to write the whole file back to disk at the end.

	ma = (void *) new char[ mm_size ];
	{
	    int nr = fread( ma, sizeof(char), mm_size, fp );

	    if (nr != mm_size) {
		throw( "Can't fake mmap on this file/system." );
	    }
	}
#else
	ma = mmap( 0, (size_t) mm_size, PROT_READ | PROT_WRITE,
		   MAP_FILE | MAP_SHARED, fd, 0 );

	if ( (caddr_t) ma == (caddr_t) -1) {
	    throw( "Can't mmap file." );
	}
#endif

    // Now we need to reconstruct the book keeping info.

	toc_entries     = ((Indx_File_Header *) ma)->entries;
	max_toc_entries = ((Indx_File_Header *) ma)->max_entries;

    // Rebuild the table of contents.

	for( i=0; i < toc_entries; i++ ) {
	    TOC_entry<Key> ent;
	    memcpy( &ent, (char *) ma + sizeof(Indx_File_Header) +
		    i * sizeof(TOC_entry<Key>), sizeof(TOC_entry<Key>) );

	// fread( (char *) &ent, sizeof(TOC_entry<Key>), 1, fp );
	    toc[i] = ent;
	}
	break;

    default:
	cerr << "\n\nERROR in Indx_File<T>: unrecognized file mode.\n\n";
	break;
    }
}

//---------------------------------------------------------------------------//
// ~Indx_File()

// I think this is bogus.  Should only write the header back if this
// file was opened for writing.  Needs work.
//---------------------------------------------------------------------------//

template<class Key>
Indx_File<Key>::~Indx_File()
{
    if (!active) return;

    Assert(!write_in_progress);	// Could toss this in favor of the
				// following check, or could throw(foo)

// Check that we're not in the middle of making an entry.

    if (write_in_progress) {
	cout << "Indx_File<T> dtor, error\n" <<flush;
	cerr << "Indx_File<T>: A write is in progress, this is not the\n";
	cerr << "right time to be closing the file.  Behavior undefined.\n";
	cerr << flush;
    }

    if (access_mode == WRITE || access_mode == READ_WRITE) {
    // Dump the header.

	header.entries = toc_entries;
	if (header.entries)
	    header.file_length = toc[toc_entries-1].base +
		toc[toc_entries-1].len;

	fseek( fp, 0, SEEK_SET );

	write_header();
	write_toc();
    }

// close the file.

    if (access_mode == RW_MMAP) {
#ifdef MAP_SHARED_IS_BUSTED
    // rewind file and write contents out.
	rewind( fp );
	int nw = fwrite( ma, sizeof(char), mm_size, fp );

	if (nw != mm_size) {
	    throw( "Unable to properly close the faked mmap file." );
	}

	delete[] (char *) ma;
#else
	if (munmap( (caddr_t) ma, mm_size )) {
	    throw( "Error when attempting to unmap file." );
	}
#endif
    }

    fflush(fp);
    fclose(fp);
}

//---------------------------------------------------------------------------//
// Methods for writing the file.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Initiate the writing of a new record to the file.
//---------------------------------------------------------------------------//

template<class Key>
int Indx_File<Key>::Start_new_entry( const Key& key )
{
    Insist( active, "This node is not participating in file I/O." );
    Insist( access_mode != RW_MMAP,
	    "RW_MMAP is not for extending existing files." );

    if (toc_entries >= max_toc_entries) {
	cerr << "Index is full.  Consider expanding the table.\n";
	return 0;			// Instead of throw(overflow).
    }

    ce = toc_entries++;

    toc[ce].key = key;

    if (ce == 0)
	toc[ce].base = data_base;
    else
	toc[ce].base = toc[ce-1].base + toc[ce-1].len;

    fseek( fp, toc[ce].base, SEEK_SET );

    write_in_progress = 1;
    bytes_written = 0;

    return 1;
}

//---------------------------------------------------------------------------//
// Indicate that the current record is now completely written.
//---------------------------------------------------------------------------//

template<class Key>
void Indx_File<Key>::End_of_entry()
{
    Insist( active, "This node is not participating in file I/O." );
    Insist( write_in_progress,
	    "Indx_File<Key>::End_of_entry(), Out of sequence." );

    fflush( fp );

    toc[ce].len = bytes_written;
    write_in_progress = 0;
    bytes_written = 0;
}

// Should this next function return len, or the number actually
// written, which might be less if file error?  Probably the latter,
// have to think on it some more.

//---------------------------------------------------------------------------//
// Write data to file.
//---------------------------------------------------------------------------//

template<class Key>
int Indx_File<Key>::write( void *p, int len, int n /*=1*/ )
{
    Insist( active, "This node is not participating in file I/O." );
    Insist( write_in_progress, "Out of sequence." );

// write len bytes starting from p, to file.

    int nn = fwrite( (char *) p, len, n, fp );
    if ( nn != n ) {
	cerr << "Unable to write entire record.\n";
	cerr << nn << " records written, " << n << " expected.\n\n";
    }

    bytes_written += nn * len;

    return bytes_written;
}

//---------------------------------------------------------------------------//
// Working with the index.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Return 1 if able to set as user requests, 0 if not.  In the later
// case, further read operations on the file are undefined.  User
// should pay attention to this result value.
//---------------------------------------------------------------------------//

template<class Key>
int Indx_File<Key>::Set_TOC_pointer( int point /*=0*/,
				     int direction /*=1*/ )
{
    Insist( active, "This node is not participating in file I/O." );

// Should probably convert this to be an exception. ...

    if ( point < 0 || point >= toc_entries )
      return 0;

    ce = point;
    search_dir = (direction >=0) ? 1 : -1;
    fseek( fp, toc[ce].base, SEEK_SET );
    return 1;
}

//---------------------------------------------------------------------------//
// Go to first record in the file.  Sets both the key and the file pointer.
//---------------------------------------------------------------------------//

template<class Key>
int Indx_File<Key>::First( Key& key )
{
    Insist( active, "This node is not participating in file I/O." );

    ce = 0;
    search_dir = 1;

    key = toc[ce].key;
    fseek( fp, toc[ce].base, SEEK_SET );
    return ce;
}

//---------------------------------------------------------------------------//
// Go to the next record in the file.  Can step in either direction.
//---------------------------------------------------------------------------//

template<class Key>
int Indx_File<Key>::Next( Key& key )
{
    Insist( active, "This node is not participating in file I/O." );

    ce += search_dir;
    if (ce < 0 || ce >= toc_entries)
	return 0;

    key = toc[ce].key;
    fseek( fp, toc[ce].base, SEEK_SET );
    return ce;
}

//---------------------------------------------------------------------------//
// Search in the previously specified (or default) direction for an
// element with the given key.  If wrap is enabled (disabled by
// default), then keep going till full circle.
//---------------------------------------------------------------------------//

template<class Key>
int Indx_File<Key>::Locate( const Key& key, int wrap /*=0*/  )
{
    Insist( active, "This node is not participating in file I/O." );

    int pe = ce;
    int done = 0;
    int idx = -1;		// The failure index.

    if (wrap) {
	while( !done ) {
	// Increment "current" element immediately so that we don't "find
	// ourselves", and thus stall.

	    ce += search_dir;

	    if ( toc[ce].key == key ) {
		fseek( fp, toc[ce].base, SEEK_SET );
		idx = ce;
		return idx;
	    }

	    if (ce < 0) {
		ce = toc_entries-1;
		continue;
	    }
	    if (ce >= toc_entries) {
		ce = -1;	// Will be incremented to zero as first
				// operation of the next iteration.
		continue;
	    }
	    if (done)
		return idx;

	// We've already adjusted ce, so if we're back to pe,
	// we've wrapped without finding the key, so reset and
	// exit. (code below).

	    if (ce == pe)
		break;
	}

    } else {
	for( ; ce >= 0 && ce < toc_entries; /*ce += search_dir*/ ) {
	// Increment "current" element first, so we don't "find ourself".

	    ce += search_dir;

	    if ( toc[ce].key == key ) {
		fseek( fp, toc[ce].base, SEEK_SET );
		idx = ce;
		return idx;
	    }
	}
    }

    ce = pe;			// Don't move ce if no match.
    return idx;			// Indicate failure.
}

//---------------------------------------------------------------------------//
// Search for the entry with the specified key, and return it's memory
// address. Only for use with files opened with mmap.
//---------------------------------------------------------------------------//

template<class Key>
void *Indx_File<Key>::mm_Locate( const Key& key, int& l, int wrap /*=0*/ )
{
    Insist( active, "This node is not participating in file I/O." );
    Assert( access_mode == IDX::RW_MMAP );

    int idx = Locate( key, wrap );

    return mm_loc_indx( idx, l );
}

//---------------------------------------------------------------------------//
// Return memory address of an entry, specified by index.
//---------------------------------------------------------------------------//

template<class Key>
void *Indx_File<Key>::mm_loc_indx( int e, int& l )
{
    Insist( active, "This node is not participating in file I/O." );
    Assert( access_mode == IDX::RW_MMAP );

    l = len( e );
    return (void *) ( (char *) ma + toc[e].base );
}

//---------------------------------------------------------------------------//
// Methods to read the file.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Block read data from the file.
//---------------------------------------------------------------------------//

template<class Key>
int Indx_File<Key>::read( void *p, int len )
{
    Insist( active, "This node is not participating in file I/O." );

    int n = fread( (char *) p, 1, len, fp );

    if (n != len ) cerr << "partial read !!!\n";

    return n;
}

template<class Key>
void Indx_File<Key>::write_header()
{
    int n = fwrite( (char *) &header, sizeof(header), 1, fp );

    if (n != 1)
	cerr << "Unable to correctly write header to disk.\n";
}

template<class Key>
void Indx_File<Key>::write_toc()
{
    for( int i=0; i < toc_entries; i++ ) {
	int n = fwrite( (char *) &toc[i], sizeof(TOC_entry<Key>), 1, fp );
	if (n != 1)
	    cerr << "Unable to write toc entry " << i << " to disk.\n"
		 << flush;
    }
}

template<class Key>
int Indx_File<Key>::get_header_size()
{
    return sizeof(header);
}

template<class Key>
int Indx_File<Key>::get_toc_entry_size()
{
    return sizeof(TOC_entry<Key>);
}

template<class Key>
void Indx_File<Key>::read_header()
{
    fread( (char *) &header, sizeof(header), 1, fp );
}

template<class Key>
void Indx_File<Key>::read_toc()
{
    for( int i=0; i < header.entries; i++ ) {
	TOC_entry<Key> ent;
	fread( (char *) &ent, sizeof(TOC_entry<Key>), 1, fp );
	toc[i] = ent;
    }
}

//===========================================================================//
// Need a special derived class so we can ensure the header gets written out
// in a machine architecture independent way.  Grrr.
//===========================================================================//

// Wooops!  The pathetic HP compiler is failing to call virtual functions
// correctly when Indx_File<int> is subclassed by Indx_File_int.  So, at
// tremendous cost in grossness, I am copying Indx_File<int> to
// Indx_File_int, and will perform horrid hacky hacks from now (10/27/95)
// until after graduation.  Man, I hate pathetic compilers...

//---------------------------------------------------------------------------//
// Constructor.  Open the file, do initial setup, etc.
//---------------------------------------------------------------------------//

Indx_File_int::Indx_File_int( const String& _fname, int _amode, 
			      int _max_toc_entries /*=0*/,
			      int _active /*=1*/ )
    : fname(_fname), access_mode(_amode),
      max_toc_entries(_max_toc_entries),
      toc(max_toc_entries), active(_active)
{
    if (!active) return;

    ce = toc_entries = 0;
    search_dir = 1;
    write_in_progress = 0;

    header_size = get_header_size();

    toc_base = 0 + header_size;

// Now see about openning the file.

    switch(access_mode) {

    case READ:
	fp = fopen( fname, "r" );

	Insist( fp, "Unable to open Indx_File<T> for reading.\n" );

	read_header();

    // Do some sanity checks.

	if ( header.entry_size != get_toc_entry_size() )
	    cerr << "bogus file.\n";

	max_toc_entries = header.max_entries;

    // Read in the index table.

	read_toc();
	toc_entries = header.entries;

	break;

    case WRITE:
	fp = fopen( fname, "w" );
	if (!fp)
	    cerr << "Indx_File<T>: unable to open " << fname
		 << " for writing.\n\n";

    // Write header, indicating size of TOC_Entry, # thereof, base of data
    // area, etc.

	header.entry_size = get_toc_entry_size();
	if (max_toc_entries)
	    header.max_entries = max_toc_entries;
	else {
	    cerr << "max_entries needs to be specified.  Defautling.\n";
	    header.max_entries = max_toc_entries = 20;
	}
	header.entries = 0;
	header.file_length = 0;
	write_header();

    // Could think about leaving a hole, maybe even a big one, to make
    // expansion of the toc easier.

	data_base = toc_base + max_toc_entries * get_toc_entry_size();
	break;

    case RW_MMAP:
    // This is used to allow read/write access to the file via mmap.  Such
    // use does not allow extending the file, but is useful for "rearranging"
    // data in an existing file.  The currently envisioned use for this is to
    // allow particle sorting in particle history files.

	if ( (fd = open( fname, O_RDWR )) < 0 ) {
	    throw( "Can't open file for read/write." );
	}
	fp = fdopen( fd, "rb+" );

	struct stat statbuf;
	fstat( fd, &statbuf );
	mm_size = (int) statbuf.st_size;

#ifdef MAP_SHARED_IS_BUSTED
    // Just one problem.  Linux (perhaps other unices?) does not support the
    // combination of PROT_WRITE and MAP_SHARED.  In our case we don't really
    // need to share with other concurrent processes.  We just want to use
    // MAP_SHARED in order to get any changes we make communicated back to
    // the disk file.  If you use MAP_PRIVATE, the changes don't go back to
    // the file.  So, we fake it.  We do this by reading the file into a
    // memory buffer, handle the access functions normally, and then be sure
    // to write the whole file back to disk at the end.

	ma = (void *) new char[ mm_size ];
	{
	    int nr = fread( ma, sizeof(char), mm_size, fp );

	    if (nr != mm_size) {
		throw( "Can't fake mmap on this file/system." );
	    }
	}
#else
	ma = mmap( 0, (size_t) mm_size, PROT_READ | PROT_WRITE,
		   MAP_FILE | MAP_SHARED, fd, 0 );

	if ( (caddr_t) ma == (caddr_t) -1) {
	    throw( "Can't mmap file." );
	}
#endif

    // Now we need to reconstruct the book keeping info.

	rewind(fp);
	read_header();

	max_toc_entries = header.max_entries;
	read_toc();
	toc_entries = header.entries;
	break;

    default:
	cerr << "\n\nERROR in Indx_File_int: unrecognized file mode.\n\n";
	break;
    }
}

//---------------------------------------------------------------------------//
// ~Indx_File_int()

// I think this is bogus.  Should only write the header back if this
// file was opened for writing.  Needs work.
//---------------------------------------------------------------------------//

Indx_File_int::~Indx_File_int()
{
    if (!active) return;

    Assert(!write_in_progress);	// Could toss this in favor of the
				// following check, or could throw(foo)

// Check that we're not in the middle of making an entry.

    if (write_in_progress) {
	cout << "Indx_File<T> dtor, error\n" <<flush;
	cerr << "Indx_File<T>: A write is in progress, this is not the\n";
	cerr << "right time to be closing the file.  Behavior undefined.\n";
	cerr << flush;
    }

    if (access_mode == WRITE || access_mode == READ_WRITE) {
    // Dump the header.

	header.entries = toc_entries;
	header.file_length = toc[toc_entries-1].base +
	    toc[toc_entries-1].len;

	fseek( fp, 0, SEEK_SET );

	write_header();
	write_toc();
    }

    if (access_mode == RW_MMAP) {
#ifdef MAP_SHARED_IS_BUSTED
    // rewind file and write contents out.
	rewind( fp );
	int nw = fwrite( ma, sizeof(char), mm_size, fp );

	if (nw != mm_size) {
	    throw( "Unable to properly close the faked mmap file." );
	}

	delete[] (char *) ma;
#else
	if (munmap( (caddr_t) ma, mm_size )) {
	    throw( "Error when attempting to unmap file." );
	}
#endif
    }

    fflush(fp);
    fclose(fp);
}

//---------------------------------------------------------------------------//
// Methods for writing the file.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Initiate the writing of a new record to the file.
//---------------------------------------------------------------------------//

int Indx_File_int::Start_new_entry( int key )
{
    Insist( active, "This node is not participating in file I/O." );
    Insist( access_mode != RW_MMAP,
	    "RW_MMAP is not for extending existing files." );
/*
    if (toc_entries >= max_toc_entries) {
	cerr << "Index is full.  Consider expanding the table.\n";
	return 0;			// Instead of throw(overflow).
    }
    */
    Assert( access_mode == IDX::WRITE );
    Insist( toc_entries < max_toc_entries,
	   "Index is full, you need to request a larger table.\n" );

    ce = toc_entries++;

    toc[ce].key = key;

    if (ce == 0)
	toc[ce].base = data_base;
    else
	toc[ce].base = toc[ce-1].base + toc[ce-1].len;

    fseek( fp, toc[ce].base, SEEK_SET );

    write_in_progress = 1;
    bytes_written = 0;

    return 1;
}

//---------------------------------------------------------------------------//
// Indicate that the current record is now completely written.
//---------------------------------------------------------------------------//

void Indx_File_int::End_of_entry()
{
    Insist( active, "This node is not participating in file I/O." );
    Insist( write_in_progress,
	    "Indx_File_int::End_of_entry(), Out of sequence." );

    fflush( fp );

    toc[ce].len = bytes_written;
    write_in_progress = 0;
    bytes_written = 0;
}

// Should this next function return len, or the number actually
// written, which might be less if file error?  Probably the latter,
// have to think on it some more.

//---------------------------------------------------------------------------//
// Write data to file.
//---------------------------------------------------------------------------//

int Indx_File_int::write( void *p, int len, int n /*=1*/ )
{
    Insist( active, "This node is not participating in file I/O." );
    Insist( write_in_progress, "Out of sequence." );

// write len bytes starting from p, to file.

    int nn = fwrite( (char *) p, len, n, fp );
    if ( nn != n ) {
	cerr << "Unable to write entire record.\n";
	cerr << nn << " records written, " << n << " expected.\n\n";
    }

    bytes_written += nn * len;

    return bytes_written;
}

//---------------------------------------------------------------------------//
// Working with the index.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Return 1 if able to set as user requests, 0 if not.  In the later
// case, further read operations on the file are undefined.  User
// should pay attention to this result value.
//---------------------------------------------------------------------------//

int Indx_File_int::Set_TOC_pointer( int point /*=0*/,
				     int direction /*=1*/ )
{
    Insist( active, "This node is not participating in file I/O." );

// Should probably convert this to be an exception. ...

    if ( point < 0 || point >= toc_entries )
      return 0;

    ce = point;
    search_dir = (direction >=0) ? 1 : -1;
    fseek( fp, toc[ce].base, SEEK_SET );
    return 1;
}

//---------------------------------------------------------------------------//
// Go to first record in the file.  Sets both the key and the file pointer.
//---------------------------------------------------------------------------//

int Indx_File_int::First( int& key )
{
    Insist( active, "This node is not participating in file I/O." );

    ce = 0;
    search_dir = 1;

    key = toc[ce].key;
    fseek( fp, toc[ce].base, SEEK_SET );
    return ce;
}

//---------------------------------------------------------------------------//
// Go to the next record in the file.  Can step in either direction.
//---------------------------------------------------------------------------//

int Indx_File_int::Next( int& key )
{
    Insist( active, "This node is not participating in file I/O." );

    ce += search_dir;
    if (ce < 0 || ce >= toc_entries)
	return 0;

    key = toc[ce].key;
    fseek( fp, toc[ce].base, SEEK_SET );
    return ce;
}

//---------------------------------------------------------------------------//
// Search in the previously specified (or default) direction for an
// element with the given key.  If wrap is enabled (disabled by
// default), then keep going till full circle.
//---------------------------------------------------------------------------//

int Indx_File_int::Locate( int key, int wrap /*=0*/  )
{
    Insist( active, "This node is not participating in file I/O." );

    int pe = ce;
    int done = 0;
    int idx = -1;		// The failure index.

    if (wrap) {
	while( !done ) {
	// Increment "current" element immediately so that we don't "find
	// ourselves", and thus stall.

	    ce += search_dir;

	    if ( toc[ce].key == key ) {
		fseek( fp, toc[ce].base, SEEK_SET );
		idx = ce;
		return idx;
	    }

	    if (ce < 0) {
		ce = toc_entries-1;
		continue;
	    }
	    if (ce >= toc_entries) {
		ce = -1;	// Will be incremented to zero as first
				// operation of the next iteration.
		continue;
	    }
	    if (done)
		return idx;

	// We've already adjusted ce, so if we're back to pe,
	// we've wrapped without finding the key, so reset and
	// exit. (code below).

	    if (ce == pe)
		break;
	}

    } else {
	for( ; ce >= 0 && ce < toc_entries; /*ce += search_dir*/ ) {
	// Increment "current" element first, so we don't "find ourself".

	    ce += search_dir;

	    if ( toc[ce].key == key ) {
		fseek( fp, toc[ce].base, SEEK_SET );
		idx = ce;
		return idx;
	    }
	}
    }

    ce = pe;			// Don't move ce if no match.
    return idx;			// Indicate failure.
}

//---------------------------------------------------------------------------//
// Search for the entry with the specified key, and return it's memory
// address. Only for use with files opened with mmap.
//---------------------------------------------------------------------------//

void *Indx_File_int::mm_Locate( int key, int& l, int wrap /*=0*/ )
{
    Insist( active, "This node is not participating in file I/O." );
    Assert( access_mode == IDX::RW_MMAP );

    int idx = Locate( key, wrap );

    return mm_loc_indx( idx, l );
}

//---------------------------------------------------------------------------//
// Return memory address of an entry, specified by index.
//---------------------------------------------------------------------------//

void *Indx_File_int::mm_loc_indx( int e, int& l )
{
    Insist( active, "This node is not participating in file I/O." );
    Assert( access_mode == IDX::RW_MMAP );

    l = len( e );
    return (void *) ( (char *) ma + toc[e].base );
}

//---------------------------------------------------------------------------//
// Methods to read the file.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Block read data from the file.
//---------------------------------------------------------------------------//

int Indx_File_int::read( void *p, int len )
{
    Insist( active, "This node is not participating in file I/O." );

    int n = fread( (char *) p, 1, len, fp );

    if (n != len ) cerr << "partial read !!!\n";

    return n;
}

void Indx_File_int::write_header()
{
#if 0
    cout << "Entering Indx_File_int::write_header()\n";
    cout << "header.entry_size =" << header.entry_size << endl;
    cout << "header.max_entries =" << header.max_entries << endl;
    cout << "header.entries =" << header.entries << endl;
    cout << "header.file_length =" << header.file_length << endl;
#endif

    mbs.clear();

    mbs.wr( header.entry_size  );
    mbs.wr( header.max_entries );
    mbs.wr( header.entries     );
    mbs.wr( header.file_length );

    Assert( mbs.size() == 20 );

    int n = fwrite( (char *) mbs.buf(), mbs.size(), 1, fp );

    if (n != 1)
	cerr << "Unable to correctly write header to disk.\n";
}

void Indx_File_int::write_toc()
{
    for( int i=0; i < toc_entries; i++ ) {
	mbs.clear();
	mbs.wr( toc[i].key  );
	mbs.wr( toc[i].base );
	mbs.wr( toc[i].len  );
	int n = fwrite( (char *) mbs.buf(), mbs.size(), 1, fp );
	if (n != 1)
	    cerr << "Unable to write toc entry " << i << " to disk.\n"
		 << flush;
    }
}

int Indx_File_int::get_header_size()
{
    int siz_int = mbs.sizeof_int();

    return (siz_int + 1) * 4;
}

int Indx_File_int::get_toc_entry_size()
{
    int siz_int = mbs.sizeof_int();

    return (siz_int + 1) * 3;
}

void Indx_File_int::read_header()
{
    mbs.reset();
    fread( (char *) mbs.buf(), get_header_size(), 1, fp );

    mbs.rd( header.entry_size  );
    mbs.rd( header.max_entries );
    mbs.rd( header.entries     );
    mbs.rd( header.file_length );
}

void Indx_File_int::read_toc()
{
    for( int i=0; i < header.entries; i++ ) {
	TOC_entry<int> ent;
	mbs.clear();
	fread( (char *) mbs.buf(), get_toc_entry_size(), 1, fp );

	mbs.rd( ent.key  );
	mbs.rd( ent.base );
	mbs.rd( ent.len  );

	toc[i] = ent;
    }
}

//---------------------------------------------------------------------------//
//                              end of Indx_File.cc
//---------------------------------------------------------------------------//
