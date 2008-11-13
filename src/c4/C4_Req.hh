//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Req.hh
 * \author Thomas M. Evans, Geoffrey Furnish
 * \date   Thu Jun  2 09:54:02 2005
 * \brief  C4_Req class definition.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_C4_Req_hh
#define c4_C4_Req_hh

// C4 package configure
#include "c4/config.h"

#ifdef C4_MPI
#include "c4_mpi.h"
#endif

namespace rtt_c4
{
//===========================================================================//
/*!
 * \class C4_ReqRefRep
 * \brief Handle for non-blocking message requests.
 *
 * This class provides an encapsulator for the message requests (MPI) which
 * are produced by non blocking calls.  This class automatically waits for the
 * message to complete when the containing object goes out of scope, thus
 * plugging one of the easiest types of programming errors with non blocking
 * messaging.  Reference counting is used so that these may be passed by value
 * without accidentally triggering a program stall.
 *
 * This is a "work" class. The user interface for requests is provided by
 * rtt_c4::C4_Req.
 */
//===========================================================================//

class C4_ReqRefRep {

    friend class C4_Req;
    
    int n;
    int assigned;
    
#ifdef C4_MPI
    MPI_Status  s;
    MPI_Request r;
#endif

  private:

    // Disallowed methods

    C4_ReqRefRep( const C4_ReqRefRep& rep );
    C4_ReqRefRep& operator=( const C4_ReqRefRep& rep );

    // Private default ctor and dtor for access from C4_Req only.
    
    C4_ReqRefRep();
    ~C4_ReqRefRep();

  public:
    
    void wait();
    void free();

    bool complete();

    unsigned count();

    int inuse() const { return assigned; }

  private:
    
    void set()      { assigned = 1; }
    void clear()    { assigned = 0; }
};

//===========================================================================//
/*!
 * \class C4_Req
 * \brief Non-blocking communication request class.
 *
 * This class provides an encapsulator for the message requests (MPI) which
 * are produced by non blocking calls.  This class automatically waits for the
 * message to complete when the containing object goes out of scope, thus
 * plugging one of the easiest types of programming errors with non blocking
 * messaging.  Reference counting is used so that these may be passed by value
 * without accidentally triggering a program stall.
 *
 * This class provides an interface for non-blocking request handles that
 * should be used by users.  
 */
//===========================================================================//

class C4_Req {

    //! Request handle.
    C4_ReqRefRep *p;

  public:
    
    C4_Req();
    C4_Req( const C4_Req& req );
    ~C4_Req();
    C4_Req& operator=( const C4_Req& req );

    //! \brief Equivalence operator
    bool operator==(const C4_Req& right ) {
        return (p == right.p); }
    bool operator!=(const C4_Req& right ) {
        return (p != right.p ); }

    void wait()     { p->wait(); }
    void free()     { p->free(); }

    bool complete() { return p->complete(); }  // Should be const?

    unsigned count() { return p->count(); }

    int inuse() const { return p->inuse(); }

  private:

    void set()      { p->set(); }

    // Private access to the C4_ReqRefRep internals.
    
#ifdef C4_MPI
// Commented out by KT on 2005 Nov 1    
// It does not appear that this member function serves any purpose.  It is not
// used by C4_Req or friend class C4_ReqReqRef and it is private.
// else can use it.
//    MPI_Status  &s() { return p->s; }
    MPI_Request &r() { return p->r; }
#endif

    // FRIENDSHIP
    
    // Specific friend C4 functions that may need to manipulate the
    // C4_ReqRefRep internals.

#ifdef C4_MPI
    template<class T>
    friend C4_Req send_async(const T *buf, int nels, int dest, int tag);
    
    template<class T>
    friend C4_Req receive_async(T *buf, int nels, int source, int tag);

    template<class T>
    friend void send_async(C4_Req &r, const T *buf, int nels, int dest,
			   int tag);
    template<class T>
    friend void receive_async(C4_Req &r, T *buf, int nels, int source, 
			      int tag);
#endif
};

} // end namespace rtt_c4

#endif // c4_C4_Req_hh

//---------------------------------------------------------------------------//
//              end of c4/C4_Req.hh
//---------------------------------------------------------------------------//
