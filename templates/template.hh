//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/<class>.hh
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_<spkg>_<class>_hh
#define rtt_<spkg>_<class>_hh

namespace rtt_<spkg>
{

//===========================================================================//
/*!
 * \class <class>
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at http://www.doxygen.org.
 *
 * \sa <class>.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*! 
 * \example <pkg>/test/<pkg>_test.cc 
 * 
 * description of example
 */
// revision history:
// -----------------
// 0) (<date>) <user>: original
// 
//===========================================================================//

class <class> 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! default constructors
    <class>();

    //! copy constructor (the long doxygen description is in the .cc file)
    <class>(const <class> &rhs);

    //! destructor
    ~<class>();

    // MANIPULATORS
    
    //! Assignment operator for <class>
    <class>& operator=(const <class> &rhs);

    // ACCESSORS

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA

};

} // end namespace rtt_<spkg>

#endif // rtt_<spkg>_<class>_hh

//---------------------------------------------------------------------------//
//              end of <pkg>/<class>.hh
//---------------------------------------------------------------------------//
