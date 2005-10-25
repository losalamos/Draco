/*!
  \file    Thin_Ptr.hh
  \author  Paul Henning
  \brief   Declaration of class Thin_Ptr
  \note    Copyright &copy; 2005 The Regents of the University of California.
  \version $Id$
*/

#ifndef rtt_dsxx_Thin_Ptr_hh
#define rtt_dsxx_Thin_Ptr_hh


namespace rtt_dsxx
{

//! An unchecked pointer-like class for scalars
/*!
  This provides the same interface as rtt_dsxx::Safe_Ptr, but should compile
  out to be raw pointers.

  Except for strange situations, you should use the \c DBC_Ptr macro in
  DBC_Ptr.hh instead of using this class explicitly.

  \sa Safe_Ptr.hh
  \sa DBC_Ptr.hh
*/
template<class T> class Thin_Ptr
{
    template<class X> friend class Thin_Ptr;
  public:
    //! Default constructor.
    /*! The pointer is initialized to \c NULL. */
    Thin_Ptr():d_ptr(0) {}

    //! Explicit constructor from pointer.
    /*! The safest way to use this constructor is if the \p ptr parameter is
      really a temporary, e.g.
      \code
      Thin_Ptr<Some_Class> a(new Some_Class);
      \endcode
      If it is really a normal raw pointer, make sure not to do anything to
      the pointed-to object through that raw pointer.
    */
    explicit Thin_Ptr(T * const ptr) : d_ptr(ptr) {}


    //! Explicit constructor from pointer of type X
    /*! 
      This is to allow things like 
      \code
      Thin_Ptr<Base_Class> bc(new Derived_Class);
      \endcode
      The comments in the first non-trivial ctor also hold.
    */
    template<class X> explicit Thin_Ptr(X * const ptr) : d_ptr(ptr) {}


    //! Assignment operator for type T *.
    Thin_Ptr<T>& operator=(T * const p_in) 
    {
	d_ptr = p_in;
	return *this;
    }

    //! Natural assignment operator
    /*!
      This was added because you can't dynamic cast PODs, as is done in the
      member template version of op=
    */
    Thin_Ptr<T>& operator=(Thin_Ptr<T> const & ptr) 
    {
	d_ptr = ptr.d_ptr;
	return *this;
    }
    
    //! Assignment from another class
    /*! 
      This is to allow things like 
      \code
      Thin_Ptr<Base_Class> bc;
      Thin_Ptr<Derived_Class> dc(new Derived_Class);
      bc = dc;
      \endcode
      The comments in the first non-trivial ctor also hold.
    */
    template<class X> Thin_Ptr<T>& operator=(Thin_Ptr<X> const & ptr) 
    {
	d_ptr = dynamic_cast<T*>(ptr.d_ptr);
	return *this;
    }

    //! Sets the pointer to \c NULL
    void release_data() { d_ptr = 0; }

    //! Calls \c delete on the pointer
    void delete_data() { delete d_ptr; d_ptr = 0; }

    //! Dereference operator
    T& operator*() const { return *d_ptr; }

    //! Access operator
    T* operator->() const { return d_ptr; }

    //! Boolean conversion operator.
    operator bool() const 
    { return d_ptr != 0; }

    //! Operator not.
    bool operator!() const 
    { return d_ptr == 0; }

    //! Equality operator for T*.
    bool operator==(const T *ptr_in) const 
    { return d_ptr == ptr_in; }

    //! Inequality operator for T*.
    bool operator!=(const T *ptr_in) const 
    { return d_ptr != ptr_in; }

    //! Equality operator for Thin_Ptr<T>.
    bool operator==(const Thin_Ptr<T> &sp_in) const 
    { return d_ptr == sp_in.d_ptr; }

    //! Inequality operator for Thin_Ptr<T>.
    bool operator!=(const Thin_Ptr<T> &sp_in) const 
    { return d_ptr != sp_in.d_ptr; }


  private:
    T* d_ptr;
};

}


#endif
