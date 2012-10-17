//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/path.cc
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright © 2011 Los Alamos National Security, LLC
 *         All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#include "path.hh"
#include "SystemCall.hh"
#include "Assert.hh"

#include <cerrno>       // errno
#include <cstring>      // strerror
#include <cstdlib>      // realpath
#include <cstdio>       // remove
#include <sstream>
#include <iostream>
#include <sys/stat.h>   // stat
#ifdef UNIX
#include <dirent.h> // struct DIR
#endif

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*!
 * \brief Helper function to extract path information from a filename.
 * \param fqName A fully qualified filename (/path/to/the/unit/test)
 * \return filename only, or path to file only.
 *
 * This function expects a fully qualfied name of a unit test (e.g.:
 * argv[0]).  It strips off the path and returns the name of the unit test.
 *
 * Options:
 *    FC_PATH        Return path portion only for fqName
 *    FC_ABSOLUTE    not implemented.
 *    FC_NAME        Return filename w/o path for fqName
 *    FC_EXT         not implemented.
 *    FC_NAME_WE     not implemented.
 *    FC_REALPATH    resolve all symlinks
 *    FC_LASTVALUE
 */
std::string getFilenameComponent( std::string const & fqName,
                                  FilenameComponent fc )
{
    using std::string;
    string retVal;
    string::size_type idx;
    std::string fullName( fqName );
    
    switch( fc )
    {
        case FC_PATH :
            // if fqName is a directory and ends with "/", trim the trailing
            // dirSep.
            if( fqName.rfind( rtt_dsxx::UnixDirSep ) == fqName.length()-1 )
                fullName=fqName.substr(0,fqName.length()-1);
            
            idx=fullName.rfind( rtt_dsxx::UnixDirSep );
            if( idx == string::npos ) 
            {
                // Didn't find directory separator, as 2nd chance look for Windows
                // directory separator. 
                idx=fullName.rfind( rtt_dsxx::WinDirSep );
            }
            // If we still cannot find a path separator, return "./"
            if( idx == string::npos )
                retVal = string( string(".") + rtt_dsxx::dirSep );
            else
                retVal = fullName.substr(0,idx+1); 
            break;
            
        case FC_NAME :
            idx=fqName.rfind( UnixDirSep );
            if( idx == string::npos )
            {
                // Didn't find directory separator, as 2nd chance look for Windows
                // directory separator.
                idx=fqName.rfind( WinDirSep );
            }
            // If we still cannot find a path separator, return the whole
            // string as the test name.
            if( idx == string::npos )
                retVal = fqName;
            else
                retVal = fqName.substr(idx+1);
            break;

        case FC_REALPATH :
            {
                std::string path( getFilenameComponent( fqName, FC_PATH ) );
                if( draco_getstat(path) )
                {
                    // On error, return empty string.
                    retVal = std::string();
                    // retVal = draco_getcwd();
                }
                else
                {
                    retVal = draco_getrealpath(path);
                }
                break;
            }
        case FC_ABSOLUTE :
            Insist( false, "case for FC_ABSOLUTE not implemented." );
            break;
        case FC_EXT :
            Insist( false, "case for FC_EXT not implemented." );
            break;
        case FC_NAME_WE :
            Insist( false, "case for FC_NAME_WE not implemented." );
            break;

        default:
            std::ostringstream msg;
            msg << "Unknown mode for rtt_dsxx::setName(). fc = " << fc;
            Insist( false, msg.str() );
    }
    return retVal;
}


//---------------------------------------------------------------------------//
/*! \brief Does the file exist?
 */
bool fileExists( std::string const & strFilename )
{
    struct stat stFileInfo;
    bool retVal( false );
    int intStat;

    // Attempt to get the file attributes
    intStat = stat( strFilename.c_str(), &stFileInfo );
    if(intStat == 0)
    {
        // We were able to get the file attributes so the file obviously
        // exists.
        retVal = true;
    }
    else
    {
        // We were not able to get the file attributes.  This may mean that we
        // don't have permission to access the folder which contains this
        // file. If you need to do that level of checking, lookup the return
        // values of stat which will give you more details on why stat failed.
        retVal = false;
    } 
    
    return retVal;
}

bool isDirectory( std::string const & path )
{
    // If the path does not exist, then it cannot be a directory.
    if( ! fileExists(path) ) return false;
    
    struct stat stFileInfo;
    bool retVal( false );

    // Attempt to get the file attributes
    Remember( int intStat = ) stat( path.c_str(), &stFileInfo );
    Check( intStat == 0 );

    if( S_ISDIR( stFileInfo.st_mode) ) retVal=true;

    return retVal;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Recursively remove a directory.
 * \arg   path String representing the top node ofthe directory to be removed.
 * \return void
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
 */
void draco_remove( std::string const & dirname )
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
                if (entry->d_type == DT_DIR)
                    draco_remove(itemPath);
                else
                {
                    std::cout << "Deleting \"" << itemPath << "\"" << std::endl;
                    remove(itemPath.c_str());
                }
            }
        }
        closedir(dir);
    }
    
    // Now the directory is empty, finally delete the directory itself. 
    std::cout << "Deleting \"" << dirname << "\"" << std::endl;
    remove(dirname.c_str());
    return;
}

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of path.hh
//---------------------------------------------------------------------------//
