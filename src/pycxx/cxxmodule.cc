//----------------------------------*-C++-*----------------------------------//
// Copyright 1998 The Regents of the University of California. 
// All rights reserved. See LEGAL.LLNL for full text and disclaimer.
//---------------------------------------------------------------------------//

#include "CXX_Objects.h"
namespace Py {
    Type Object::type () const { 
        return Type (FromAPI(PyObject_Type (p)));
    }
    
    String Object::str () const {
        return String (FromAPI(PyObject_Str (p)));
    }
    
    String Object::repr () const { 
        return String (FromAPI(PyObject_Repr (p)));
    }
    
    std::string Object::as_string() const {
        return static_cast<std::string>(str());
    }
        
    bool Object::isType (const Type& t) const { 
        return type ().ptr() == t.ptr();
    }
    
    Char::operator String() const {
        return String(ptr());
    }

    
    // output
    std::ostream& operator<< (std::ostream& os, const Object& ob) {
        return (os << static_cast<std::string>(ob.str()));
    }  
} // end namespace Py
