//----------------------------------*-C++-*----------------------------------//
// Indx_File.hh
// Geoffrey Furnish
// 3 March 1994
//---------------------------------------------------------------------------//
// @> A class for managing a file containing indexed entries of unspecified
// @> form and content.
//---------------------------------------------------------------------------//

#ifndef __util_Indx_File_hh__
#define __util_Indx_File_hh__

#include "ds++/String.hh"
#include "ds++/DynArray.hh"

#include <stdio.h>

template<class Key>
class TOC_entry {
  public:
    Key key;
    int base;
    int len;

  public:
    TOC_entry( int bogus =0 ) { bogus++; }
    TOC_entry( const Key& k, int _base, int _len =0 );

    // We don't need copy ctor, dtor, or operator= for this simple
    // class, since the only (potential) aggregate is key, and we are
    // specifically mandating that it must be "flat".  So the default
    // versions should work just fine.

    int operator==( const TOC_entry<Key>& e ) const;
    int operator!=( const TOC_entry<Key>& e ) const {
	return !(*this == e);
    }
};

struct Indx_File_Header {
    int entry_size;
    int max_entries;
    int entries;
    int file_length;
};

class IDX {
  public:
    enum { READ, WRITE, READ_WRITE, RW_MMAP };
};

//===========================================================================//
// class Indx_File<Key> - Access records of a file using a Key.

// The purpose of this class is to provide indexed access to a file.
// For example, you want to store the phi data to disk every n
// timesteps, and then come back later and read out phi at timestep t.
// This class is supposed to be able to keep an index, or a table of
// contents, etc, of the file, allowing you to go directly to the
// point in the file holding an indexed item of a given number.
//
// The class is templated with a Key, allowing you to associate some
// sort of user defined data with each entry in the toc.  The Key must
// be some sort of a "flat" data type, as it will be copied to the
// file using fwrite or the like.  So, no embedded pointers, etc.
// Also, the key need not be unique--functions are provided to search
// for keys, and to repeat search.
//===========================================================================//

template<class Key>
class Indx_File : public IDX {

    Indx_File( const Indx_File<Key>& );
    Indx_File<Key>& operator=( const Indx_File<Key>& );

  protected:
    Indx_File_Header header;

    int header_size;
    int toc_base;
    int data_base;

    String fname;
    int access_mode;
    int toc_entries, max_toc_entries; // current #, max # of entries.
    int ce;			// current entry
    int search_dir;

    int write_in_progress, bytes_written;

    DynArray< TOC_entry<Key> > toc;

    int fd;
    FILE *fp;
    void *ma;			// memory address of file.
    int mm_size;

    int active;

  public:
    
    Indx_File( String _fname, int _amode, int _max_toc_entries =0,
	       int _active =1 );
    virtual ~Indx_File();

// Methods for writing the file.

    int Start_new_entry( const Key& key );
    void End_of_entry();
    int write( void *p, int len, int n =1 );

// Working with the index.

// (You should be able to work either using the actual integer index,
// or with the keys.)

// These deal with indexes.
    int entries() const { return toc_entries; }
    int Set_TOC_pointer( int point =0, int direction =1 );
    int Set( int point ) { return Set_TOC_pointer( point ); }

// These deal with keys.

    int First( Key& key );
    int Next( Key& key );
    int Locate( const Key& key, int wrap =0 ); // returns index of entry, or
				               // -1 if not found. 
    Key key() const { return key(ce); }
    int len() const { return len(ce); }

    Key key( int e ) const { return toc[e].key; }
    int len( int e ) const { return toc[e].len; }

    void *mm_Locate( const Key& key, int& l, int wrap =0 );
    void *mm_loc_indx( int e, int& l );

//    operator void* () { return (void *) (ce >= 0 && ce < toc_entries); }
    operator int () { return (ce >= 0 && ce < toc_entries); }

// Methods to read the file.

    int read( void *p, int len );

// Methods to concretize this class

    virtual void write_header();
    virtual void write_toc();
    virtual int get_header_size();
    virtual int get_toc_entry_size();

    virtual void read_header();
    virtual void read_toc();
};

#include "spdf/spdf_stream.hh"

class Indx_File_int : public IDX
{
    mbuf_spdf_stream mbs;

    Indx_File_int( const Indx_File_int& );
    Indx_File_int& operator=( const Indx_File_int& );

  protected:
    Indx_File_Header header;

    int header_size;
    int toc_base;
    int data_base;

    String fname;
    int access_mode;
    int toc_entries, max_toc_entries; // current #, max # of entries.
    int ce;			// current entry
    int search_dir;

    int write_in_progress, bytes_written;

    DynArray< TOC_entry<int> > toc;

    int fd;
    FILE *fp;
    void *ma;			// memory address of file.
    int mm_size;

    int active;

  public:
    
    Indx_File_int( const String& _fname, int _amode, int _max_toc_entries =0,
		   int _active =1 );
    virtual ~Indx_File_int();

// Methods for writing the file.

    int Start_new_entry( int key );
    void End_of_entry();
    int write( void *p, int len, int n =1 );

// Working with the index.

// (You should be able to work either using the actual integer index,
// or with the keys.)

// These deal with indexes.
    int entries() const { return toc_entries; }
    int Set_TOC_pointer( int point =0, int direction =1 );
    int Set( int point ) { return Set_TOC_pointer( point ); }

// These deal with keys.

    int First( int& key );
    int Next( int& key );
    int Locate( int key, int wrap =0 ); // returns index of entry, or
				               // -1 if not found. 
    int key() const { return key(ce); }
    int len() const { return len(ce); }

    int key( int e ) const { return toc[e].key; }
    int len( int e ) const { return toc[e].len; }

    void *mm_Locate( int key, int& l, int wrap =0 );
    void *mm_loc_indx( int e, int& l );

//    operator void* () { return (void *) (ce >= 0 && ce < toc_entries); }
    operator int () { return (ce >= 0 && ce < toc_entries); }

// Methods to read the file.

    int read( void *p, int len );

    void write_header();
    void write_toc();
    int get_header_size();
    int get_toc_entry_size();

    void read_header();
    void read_toc();
};

#endif                          // __util_Indx_File_hh__

//---------------------------------------------------------------------------//
//                              end of util/Indx_File.hh
//---------------------------------------------------------------------------//
