//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/path.hh
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright © 2010 Los Alamos National Security, LLC
 *         All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#ifndef __dsxx_path_hh__
#define __dsxx_path_hh__

/*!
 * \bug Consider replacing path.cc and path.hh with Boost FileSystem.
 */

#include "ds++/config.h"
#include "SystemCall.hh"  // rtt_dsxx::dirSep
#include "Assert.hh"
#include <string>
#include <iostream>
#ifdef UNIX
#include <dirent.h> // struct DIR
#include <sys/stat.h> // struct stat; S_ISDIR
#endif

namespace rtt_dsxx
{

enum FilenameComponent
{
    FC_PATH,       //!< Extract path portion of fully qualified filename
    FC_ABSOLUTE,
    FC_NAME,       //!< Extract filename portion (w/o path).
    FC_EXT,
    FC_NAME_WE,
    FC_REALPATH,
    FC_LASTVALUE
};

//---------------------------------------------------------------------------//
//! Get a specific component of a full filename.
DLL_PUBLIC std::string getFilenameComponent( std::string const & fqName,
                                             FilenameComponent   fc );

//---------------------------------------------------------------------------//
//! Does the file exist?
DLL_PUBLIC bool fileExists(  std::string const & filename );
DLL_PUBLIC bool isDirectory( std::string const & path );

//---------------------------------------------------------------------------//
//! Functor for printing all items in a directory tree
class wdtOpPrint 
{
  public:
    void operator()( std::string const & dirpath ) const {
        std::cout << dirpath << std::endl; }
};
//---------------------------------------------------------------------------//
//! Functor for removing all items in a direcotry tree
class wdtOpRemove 
{
  public:
    void operator()( std::string const & dirpath ) const
    {
        std::cout << "Removing \"" << dirpath << "\"" << std::endl;
        remove( dirpath.c_str() );
        Ensure( ! fileExists( dirpath.c_str() ) );
    }
};  

//---------------------------------------------------------------------------//
/*!
 * \brief Walk a directory tree structure, perform myOperator() action on
 *  each entry.
 * \arg dirname String representing the top node of the directory to be parsed.
 * \arg myOperator Functor that defines action to be taken on each entry in
 *  the directory. Recommend using wdtOpPrint or wdtOpRemove
 * \return void
 *
 * \sa draco_remove Helper function to recurively delete a directory and all
 * its contents.
 * \sa draco_dir_print Helper function that will print a directory and all
 * its contents.
  *
 * Sample implementation for Win32 (uses Win API which I don't want to do)
 * http://stackoverflow.com/questions/1468774/recursive-directory-deletion-with-win32
 * http://msdn.microsoft.com/en-us/library/aa365488%28VS.85%29.aspx
 *
 * Sample implementation for Unix
 * http://www.linuxquestions.org/questions/programming-9/deleting-a-directory-using-c-in-linux-248696/
 *
 * Consider using Boost.FileSystem
 * \c boost::filesystem::remove_all(path);
 *
 * \code{.cpp}
#include "boost/filesystem.hpp"
#include <iostream>
using namespace boost::filesystem;
int main()
{
  path current_dir("."); //
  for (recursive_directory_iterator iter(current_dir), end; iter != end; ++iter)
    std::cout << iter->path() << "\n";
  return 0;
}
 * \endcode
 */
template< typename T > 
void draco_walk_directory_tree( std::string const & dirname, T const & myOperator )
{
    // If file does not exist, report and continue.
    if( ! fileExists( dirname ) )
    {
        std::cout << "File/directory \"" << dirname
                  << "\"does not exist.  Continuing..." << std::endl;
        return;
    }

    // If this is not a directory, no recursion is needed:
    if( isDirectory( dirname ) )
    {
        DIR *dir;
        struct dirent *entry;
        struct stat statbuf;
        
        dir = opendir( dirname.c_str() );
        Insist(dir != NULL, "Error opendir()");

        // Loop over all entries in the directory.
        while( (entry = readdir(dir)) != NULL )
        {
            std::string d_name( entry->d_name );
            // Don't include "." or ".." entries.
            if( d_name != std::string(".") && d_name != std::string("..") )
            {
                std::string itemPath;
                if( dirname[dirname.length()-1] == UnixDirSep )
                    itemPath = dirname + d_name;
                else
                    itemPath = dirname + UnixDirSep + d_name;
                
                // if the entry is a directory, recursively delete it,
                // otherwise, delete the file

                // This implementation fails on the lightwight compute kernels
                // on CI/CT, so we will use a different technique that relies
                // on stat.
//                 if (entry->d_type == DT_DIR)
//                     draco_walk_directory_tree( itemPath, myOperator );
//                 else
//                     myOperator( itemPath );
                int error = stat( itemPath.c_str(), &statbuf );
                Check( error != -1 );
                if( S_ISDIR( statbuf.st_mode ) )
                    draco_walk_directory_tree(itemPath, myOperator);
                else
                    myOperator( itemPath );
            }
        }
        closedir(dir);
    }
    
    // Perform action on the top level entry
    myOperator(dirname);

    return;
}

//---------------------------------------------------------------------------//
//! Recursively remove a directory.
void draco_remove( std::string const & path );
//! Recursively print a directory tree.
void draco_dir_print( std::string const & path );

} // end namespace rtt_dsxx

#endif // __dsxx_path_hh__

//---------------------------------------------------------------------------//
// end of path.hh
//---------------------------------------------------------------------------//
