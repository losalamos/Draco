//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/<class>.hh
 * \author <user>
 * \date   <date>
 * \brief  <start>
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
 * commands can be found at:
 *
 * Doxygen tutorial: http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * Doxygen keywords: http://www.stack.nl/~dimitri/doxygen/commands.html
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
// 0) original
// 
//===========================================================================//

class <class> 
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    // ----------------------------------------
    // CREATORS
    // ----------------------------------------
    
    //! constructors
    <class>();
    <class>(const <class> &rhs);

    //! destructor
    ~<class>();

    // ----------------------------------------
    // MANIPULATORS
    // ----------------------------------------
    
    /*!
     * \brief Assignment operator for <class>
     *
     * This will create a copy of the instantiated object.  For now the
     * assignment operator has been disabled because we have a prototype with
     * no definition.  This prevents the compiler from creating a default
     * assignment operator.
     * 
     * \param rhs The object that is being copies
     * \return <class> A copy of the original object.
     */
    <class>& operator=(const <class> &rhs);

    // ----------------------------------------
    // ACCESSORS
    // ----------------------------------------

    /*
     * \brief Short description
     *
     * Long description
     *
     * \param c description goes here.
     * \return description goes here.
     *
     * \example test/mytest.cc
     */
    // int dummy_accessor( int c ) const;

  private:
    
    // ----------------------------------------
    // IMPLEMENTATION
    // ----------------------------------------

    // ----------------------------------------
    // DATA
    // ----------------------------------------

    /*
     * \brief short description
     *
     * Long description here.  
     */
    // int dummy_data;

};

} // end namespace rtt_<spkg>

#endif // rtt_<spkg>_<class>_hh

//---------------------------------------------------------------------------//
//              end of <pkg>/<class>.hh
//---------------------------------------------------------------------------//
