//----------------------------------*-C++-*----------------------------------//
// Copyright 1998 The Regents of the University of California. 
// All rights reserved. See LEGAL.LLNL for full text and disclaimer.
//---------------------------------------------------------------------------//

#ifndef __CXX_Objects__h
#define __CXX_Objects__h

#include "Python.h"
#include "CXX_Exception.h"
#include <iostream>
#include <strstream>
#include <string>
#include <utility>
#include <iterator>

namespace Py {
    // Forward declarations
    class Object;
    class Type;
    template<class T> class SeqBase;
    class String;
    template<class T> class MapBase;
    
    // new_reference_to also overloaded below on Object
    inline PyObject* new_reference_to(PyObject* p) {
        Py_XINCREF(p);
        return p;
    } 
    
    // returning Null() from an extension method triggers a 
    // Python exception
    inline PyObject* Null()
    {
        return (static_cast<PyObject*>(0));
    }
    
    // Nothing() is what an extension method returns is 
    // there is no other return value.
    inline PyObject* Nothing()
    {
        Py_INCREF(Py_None);
        return Py_None;
    }
    
    class FromAPI {
        // Python API routines return to you owned pointers.
        // FromAPI is a helper class to let you construct an object from one of these.
        
        // Usage: Object(FromAPI(p)) or Object x = FromAPI(p) 
        // where p is an already-owned pointer such as returned by an API routine.
    private:
        PyObject *p;
        FromAPI& operator=(const FromAPI& other); // no assignment
        FromAPI(const FromAPI& other); // no copy constructor
    public:
        explicit FromAPI(PyObject *powned) { // construct from a pointer you own, only!
            p = powned;
        }
        virtual ~FromAPI() {
            Py_XDECREF(p); // we own it, so kill it
        }
        operator PyObject*() const {return p;}
    };
    
    //===========================================================================//
    // class Object
    // The purpose of this class is to serve as the most general kind of
    // Python object, for the purpose of writing C++ extensions in Python
    // Objects hold a PyObject* which they own. This pointer is always a 
    // valid pointer to a Python object. In children we must maintain this behavior.
    //
    // Instructions on how to make your own class MyType descended from Object:
    // (0) Pick a base class, either Object or perhaps SeqBase<T> or MapBase<T>.
    //     This example assumes Object.
    
    // (1) Write a routine int MyType_Check (PyObject *) modeled after PyInt_Check,
    //     PyFloat_Check, etc.  
    
    // (2) Add method accepts: 
    //     virtual bool accepts (PyObject *pyob) const {
    //         return pyob && MyType_Check (pyob);
    //     }
    
    // (3) Include the following constructor and copy constructor
    // 
    /*
    explicit MyType (PyObject *pyob): Object(pyob) {
       validate();
    }
    MyType(const Object& other): Object(*other) {
       validate();
    }
    */
    // You may wish to add other constructors; see the classes below for examples.
    // Each constructor must use "set" to set the pointer
    // and end by validating the pointer you have created.
    
    // (4) Each class needs at least these two assignment operators:
    /*
    MyType& operator= (const Object& rhs) {
       return (*this = *rhs);
    }
    
      Mytype& operator= (PyObject* rhsp) {
         if(ptr() == rhsp) return *this;
         set(rhsp);
         return *this;
      }
    */  
    // Note on accepts: constructors call the base class 
    // version of a virtual when calling the base class constructor, 
    // so the test has to be done explicitly in a descendent.
    
    // If you are inheriting from PythonExtension<T> to define an object
    // note that it contains PythonExtension<T>::check
    // which you can use in accepts when writing a wrapper class.
    // See Demo/r.h and Demo/r.cxx for an example.
    
    class Object {
    private:
        // the pointer to the Python object 
        // Only Object sets this directly. 
        // The default constructor for Object sets it to Py_None and
        // child classes must use "set" to set it
        //
        PyObject* p;
        
    protected:
        
        void set (PyObject* pyob) {
            release();
            p = pyob;
            Py_XINCREF (p);
            validate();
        }
        
        void release () {
            Py_XDECREF (p);
            p = 0;
        }
        
        void validate() {
            // release pointer if not the right type
            if (! accepts (p)) {
                release ();
                if(PyErr_Occurred()) { // Error message already set
                    throw Exception();
                }
                // Better error message if RTTI available
                // std::string s("Error creating object of type ");
                // s += (typeid (*this)).name();
                
                throw TypeError ("CXX: type error.");
            }
        }
        
    public:
        // Constructors acquire new ownership of pointer
        explicit Object (PyObject* pyob=Py_None): p (pyob) {
            Py_XINCREF (p);
            validate();
        }
        
        // Copy constructor acquires new ownership of pointer
        Object (const Object& ob): p(ob.p) {
            Py_XINCREF (p);
            validate();
        }
        
        // Assignment acquires new ownership of pointer      
        Object& operator= (const Object& rhs) {
            set(rhs.p);
            return *this;
        }
        
        Object& operator= (PyObject* rhsp) {
            if(ptr() == rhsp) return *this;
            set (rhsp);
            return *this;
        }
        
        // Destructor 
        virtual ~Object () {
            release ();
        }
        
        // Loaning the pointer to others, retain ownership
        PyObject* operator* () const { 
            return p;
        }
        
        // Would like to call this pointer() but messes up STL in SeqBase<T>
        PyObject* ptr () const {
            return p;
        }
        
        //
        // Queries	
        //
        
        // Can pyob be used in this object's constructor?
        virtual bool accepts (PyObject *pyob) const {
            return (pyob != 0);
        }
        
        int reference_count () const { // the reference count
            return p ? p->ob_refcnt : 0;
        }
        
        Type type () const; // the type object associated with this one
        
        String str () const; // the str() representation
        
        std::string as_string() const;
        
        String repr () const; // the repr () representation
        
        bool hasAttr (const std::string& s) const {
            return PyObject_HasAttrString (p, const_cast<char*>(s.c_str())) ? true: false;
        }
        
        Object getAttr (const std::string& s) const {
            return Object(FromAPI(
                PyObject_GetAttrString (p, const_cast<char*>(s.c_str()))
                ));
        }
        
        Object getItem (const Object& key) const {
            return Object(FromAPI(PyObject_GetItem(p, *key)));
        }
        
        long hashValue () const {
            return PyObject_Hash (p);
        }
        
        // int print (FILE* fp, int flags=Py_Print_RAW) {
        //     return PyObject_Print (p, fp, flags);
        // }
        bool is(PyObject *pother) const {  // identity test
            return p == pother;
        }
        
        bool is(const Object& other) const { // identity test
            return p == other.p;
        }
        
        bool isCallable () const {
            return PyCallable_Check (p) ? true: false;
        }
        
        bool isDict () const {
            return PyDict_Check (p) ? true: false;
        }
        
        bool isList () const {
            return PyList_Check (p) ? true: false;
        }
        
        bool isMapping () const {
            return PyMapping_Check (p) ? true: false;
        }
        
        bool isNumeric () const {
            return PyNumber_Check (p) ? true: false;
        }
        
        bool isSequence () const {
            return PySequence_Check (p) ? true: false;
        }
        
        bool isTrue () const {
            return PyObject_IsTrue (p) ? true: false;
        }
        
        bool isType (const Type& t) const;
        
        bool isTuple() const {
            return PyTuple_Check(p) ? true: false;
        }
        
        bool isString() const {
            return PyString_Check(p) ? true: false;
        }
        
        // Commands
        void setAttr (const std::string& s, const Object& value) {
            if(PyObject_SetAttrString (p, const_cast<char*>(s.c_str()), *value) == -1) 
                throw AttributeError ("getAttr failed.");
        }
        
        void delAttr (const std::string& s) {
            if(PyObject_DelAttrString (p, const_cast<char*>(s.c_str())) == -1)
                throw AttributeError ("delAttr failed.");
        }
        
        // PyObject_SetItem is too weird to be using from C++
        // so it is intentionally omitted.
        
        void delItem (const Object& key) {
            //if(PyObject_DelItem(p, *key) == -1) 
            // failed to link on Windows?
            throw KeyError("delItem failed.");
        }
        // Equality and comparison use PyObject_Compare
        
        bool operator==(const Object& o2) const {
            int k = PyObject_Compare (p, *o2);
            if (PyErr_Occurred()) throw Exception();
            return k == 0;
        }
        
        bool operator!=(const Object& o2) const {
            int k = PyObject_Compare (p, *o2);
            if (PyErr_Occurred()) throw Exception();
            return k != 0;
            
        }
        
        bool operator>=(const Object& o2) const {
            int k = PyObject_Compare (p, *o2);
            if (PyErr_Occurred()) throw Exception();
            return k >= 0;
        }
        
        bool operator<=(const Object& o2) const {
            int k = PyObject_Compare (p, *o2);
            if (PyErr_Occurred()) throw Exception();
            return k <= 0;
        }
        
        bool operator<(const Object& o2) const {
            int k = PyObject_Compare (p, *o2);
            if (PyErr_Occurred()) throw Exception();
            return k < 0;
        }
        
        bool operator>(const Object& o2) const {
            int k = PyObject_Compare (p, *o2);
            if (PyErr_Occurred()) throw Exception();
            return k > 0;
        }
     };
     // End of class Object
     inline PyObject* new_reference_to(const Object& g) {
         PyObject* p = g.ptr();
         Py_XINCREF(p);
         return p;
     }
     
     std::ostream& operator<< (std::ostream& os, const Object& ob);
     
     // Class Type
     class Type: public Object {
     public:
         explicit Type (PyObject* pyob): Object(pyob) {
             validate();
         }
         
         explicit Type (const Object& ob): Object(*ob) {
             validate();
         }
         
         Type(const Type& t): Object(t) {
             validate();
         }
         
         Type& operator= (const Object& rhs) {
             return (*this = *rhs);
         }
         
         Type& operator= (PyObject* rhsp) {
             if(ptr() == rhsp) return *this;
             set (rhsp);
             return *this;
         }
         virtual bool accepts (PyObject *pyob) const {
             return pyob && PyType_Check (pyob);
         }
     };
     
     // ===============================================
     // class Int
     class Int: public Object {
     public:
         // Constructor
         explicit Int (PyObject *pyob): Object (pyob) {
             validate();
         }
         
         Int (const Int& ob): Object(*ob) {
             validate();
         }
         
         // create from long
         explicit Int (long v = 0L): Object(FromAPI(PyInt_FromLong(v))) {
             validate();
         }
         
         // create from int
         explicit Int (int v) {
             long w = v;
             set(FromAPI(PyInt_FromLong(w)));
             validate();
         }
         
         explicit Int (const Object& ob) {
             set(FromAPI(PyNumber_Int(*ob)));
             validate();
         }
         
         // Assignment acquires new ownership of pointer
         
         Int& operator= (const Object& rhs) {
             return (*this = *rhs);
         }
         
         Int& operator= (PyObject* rhsp) {
             if(ptr() == rhsp) return *this;
             set (FromAPI(PyNumber_Int(rhsp)));
             return *this;
         }
         // Membership
         virtual bool accepts (PyObject *pyob) const {
             return pyob && PyInt_Check (pyob);
         }
         // convert to long
         operator long() const {
             return PyInt_AsLong (ptr());
         }
         // assign from an int
         Int& operator= (int v) {
             *this = FromAPI(PyInt_FromLong (long(v)));
             return *this;
         }
         // assign from long
         Int& operator= (long v) {
             *this = FromAPI(PyInt_FromLong (v));
             return *this;
         }	
     };
     
     // ===============================================
     // class Long
     class Long: public Object {
     public:
         // Constructor
         explicit Long (PyObject *pyob): Object (pyob) {
             validate();
         }
         
         Long (const Long& ob): Object(*ob) {
             validate();
         }
         
         // create from long
         explicit Long (long v = 0L)
             : Object(FromAPI(PyLong_FromLong(v))) 
         {
             validate();
         }
         // create from int
         explicit Long (int v) 
             : Object(FromAPI(PyLong_FromLong(static_cast<long>(v)))) 
         {
             validate();
         }
         
         // try to create from any object
         explicit Long (const Object& ob) 
             : Object(FromAPI(PyNumber_Long(*ob))) 
         {
             validate();
         }
         
         // Assignment acquires new ownership of pointer
         
         Long& operator= (const Object& rhs) {
             return (*this = *rhs);
         }
         
         Long& operator= (PyObject* rhsp) {
             if(ptr() == rhsp) return *this;
             set (FromAPI(PyNumber_Long(rhsp)));
             return *this;
         }
         // Membership
         virtual bool accepts (PyObject *pyob) const {
             return pyob && PyLong_Check (pyob);
         }
         // convert to long
         operator long() const {
             return PyLong_AsLong (ptr());
         }
         operator double() const {
             return PyLong_AsDouble (ptr());
         }
         // assign from an int
         Long& operator= (int v) {
             *this = FromAPI(PyLong_FromLong (long(v)));
             return *this;
         }
         // assign from long
         Long& operator= (long v) {
             *this = FromAPI(PyLong_FromLong (v));
             return *this;
         }	
     };
     
     // ===============================================
     // class Float
     // No constructor from a double since it is bad to overload on a pointer
     // and a numerical type.
     class Float: public Object {
     public:
         // Constructor
         explicit Float (PyObject *pyob): Object(pyob) {
             validate();
         }
         
         Float (const Float& f): Object(f) {
             validate();
         }
         
         // make from double
         explicit Float (double v=0.0)
             :Object(FromAPI(PyFloat_FromDouble (v)))
         {
             validate();
         }
         
         // try to make from any object
         explicit Float (const Object& ob)
             :Object(FromAPI(PyNumber_Float(*ob)))
         {
             validate();
         }
         
         Float& operator= (const Object& rhs) {
             return (*this = *rhs);
         }
         
         Float& operator= (PyObject* rhsp) {
             if(ptr() == rhsp) return *this;
             set (FromAPI(PyNumber_Float(rhsp)));
             return *this;
         }
         // Membership
         virtual bool accepts (PyObject *pyob) const {
             return pyob && PyFloat_Check (pyob);
         }
         // convert to double
         operator double () const {
             return PyFloat_AsDouble (ptr());
         }
         // assign from a double
         Float& operator= (double v) {
             *this = FromAPI(PyFloat_FromDouble (v));
             return *this;
         }
         // assign from an int
         Float& operator= (int v) {
             *this = FromAPI(PyFloat_FromDouble (double(v)));
             return *this;
         }
         // assign from long
         Float& operator= (long v) {
             *this = FromAPI(PyFloat_FromDouble (double(v)));
             return *this;
         }	
         // assign from an Int
         Float& operator= (const Int& iob) {
             *this = FromAPI(PyFloat_FromDouble (double(long(iob))));
             return *this;
         }
     };
     
     // Sequences
     // Sequences are here represented as sequences of items of type T. 
     // The base class SeqBase<T> represents that.
     // In basic Python T is always "Object". 
     
     // seqref<T> is what you get if you get elements from a non-const SeqBase<T>.
     // Note: seqref<T> could probably be a nested class in SeqBase<T> but that might stress
     // some compilers needlessly. Simlarly for mapref later. 
     
     // While this class is not intended for enduser use, it needs some public
     // constructors for the benefit of the STL. 
     
     // See Scott Meyer's More Essential C++ for a description of proxies.
     // This application is even more complicated. We are doing an unusual thing
     // in having a double proxy. If we want the STL to work 
     // properly we have to compromise by storing the rvalue inside. The 
     // entire Object API is repeated so that things like s[i].isList() will 
     // work properly.
     
     template<class T>
         class seqref { 
     protected:
         SeqBase<T>& s; // the sequence
         int offset; // item number
         T the_item; // lvalue
     public:
         seqref (SeqBase<T>& seq, int j)
             : s(seq), offset(j), the_item (s.getItem(j)){}

         seqref (const seqref<T>& r)
             : s(r.s), offset(r.offset), the_item(r.the_item) {}
         
         ~seqref() {}
         
         operator T() const { // rvalue
             return the_item;
         } 
         
         seqref<T>& operator=(const seqref<T>& rhs) { //used as lvalue
             the_item = rhs.the_item;
             s.setItem(offset, the_item);
             return *this;
         }
         
         seqref<T>& operator=(const T& ob){ // used as lvalue
             the_item = ob;
             s.setItem(offset, ob);
             return *this;
         }
         
         // forward everything else to the item           
         PyObject* ptr () const {
             return the_item.ptr();
         }
         
         int reference_count () const { // the reference count
             return the_item.reference_count();
         }
         
         Type type () const {
             return the_item.type();
         }
         
         String str () const {
             return the_item.str();
         }
         
         String repr () const {
             return the_item.repr();
         }
         
         bool hasAttr (const std::string& attr_name) const {
             return the_item.hasAttr(attr_name);
         }
         
         Object getAttr (const std::string& attr_name) const {
             return the_item.getAttr(attr_name);
         }
         
         Object getItem (const Object& key) const {
             return the_item.getItem(key);
         }
         
         long hashValue () const {
             return the_item.hashValue();
         }
         
         bool isCallable () const {
             return the_item.isCallable();
         }
         
         bool isDict () const {
             return the_item.isDict();
         }
         
         bool isList () const {
             return the_item.isList();
         }
         
         bool isMapping () const {
             return the_item.isMapping();
         }
         
         bool isNumeric () const {
             return the_item.isNumeric();
         }
         
         bool isSequence () const {
             return the_item.isSequence();
         }
         
         bool isTrue () const {
             return the_item.isTrue();
         }
         
         bool isType (const Type& t) const {
             return the_item.isType (t);
         }
         
         bool isTuple() const {
             return the_item.isTuple();
         }
         
         bool isString() const {
             return the_item.isString();
         }
         // Commands
         void setAttr (const std::string& attr_name, const Object& value) {
             the_item.setAttr(attr_name, value);
         }
         
         void delAttr (const std::string& attr_name) {
             the_item.delAttr(attr_name);
         }
         
         void delItem (const Object& key) {
             the_item.delItem(key);
         }
         
         bool operator==(const Object& o2) const {
             return the_item == o2;
         }
         
         bool operator!=(const Object& o2) const {
             return the_item != o2;               
         }
         
         bool operator>=(const Object& o2) const {
             return the_item >= o2;
         }
         
         bool operator<=(const Object& o2) const {
             return the_item <= o2;
         }
         
         bool operator<(const Object& o2) const {
             return the_item < o2;
         }
         
         bool operator>(const Object& o2) const {
             return the_item > o2;
         }
     }; // end of seqref
     
     
     // class SeqBase<T>
     // ...the base class for all sequence types
     
     template<class T>
         class SeqBase: public Object {
     public:
         // STL definitions
         typedef int size_type; 
         typedef seqref<T> reference;
         typedef T const_reference;
         typedef seqref<T>* pointer;
         typedef int difference_type;
         
         virtual size_type max_size() const {
             return std::string::npos; // ?
         }
         
         virtual size_type capacity() const {
             return size();
         }
         
         virtual void swap(SeqBase<T>& c) {
             SeqBase<T> temp = c;
             c = ptr();
             set(temp.ptr());
         }
         
         virtual size_type size () const { 
             return PySequence_Length (ptr());
         }
         
         explicit SeqBase<T> ()
             :Object(FromAPI(PyTuple_New(0))) 
         {
             validate();
         }
         
         explicit SeqBase<T> (PyObject* pyob): Object(pyob) {
             validate();
         }
         
         SeqBase<T> (const Object& ob): Object(ob) {
             validate();
         }
         
         // Assignment acquires new ownership of pointer
         
         SeqBase<T>& operator= (const Object& rhs) {
             return (*this = *rhs);
         }
         
         SeqBase<T>& operator= (PyObject* rhsp) {
             if(ptr() == rhsp) return *this;
             set (rhsp);
             return *this;
         }
         
         virtual bool accepts (PyObject *pyob) const {
             return pyob && PySequence_Check (pyob);
         }
         
         size_type length () const {
             return PySequence_Length (ptr());
         }
         
         // Element access 
         const T operator[](size_type index) const {
             return getItem(index);
         }
         
         seqref<T> operator[](size_type index) {
             return seqref<T>(*this, index);
         }
         
         virtual T getItem (size_type i) const {
             return T(FromAPI(PySequence_GetItem (ptr(), i)));
         }
         
         virtual void setItem (size_type i, const T& ob) {
             if (PySequence_SetItem (ptr(), i, *ob) == -1) {
                 throw Exception();
             }
         }
         
         SeqBase<T> repeat (int count) const {
             return SeqBase<T> (FromAPI (PySequence_Repeat (ptr(), count)));
         }
         
         SeqBase<T> concat (const SeqBase<T>& other) const {
             return SeqBase<T> (FromAPI (PySequence_Concat(ptr(), *other)));
         }
         
         // more STL compatability
         const T front () const {
             return getItem(0);
         }
         
         seqref<T> front() {
             return seqref(this, 0);
         }
         
         const T back () const {
             return getItem(size()-1);
         }
         
         seqref<T> back() {
             return seqref(this, size()-1);
         }
         
         void verify_length(size_type required_size) 
         {
             if (size() != required_size) 
                 throw IndexError ("Unexpected SeqBase<T> length.");
         }
         
         void verify_length(size_type min_size, size_type max_size)
         {
             size_type n = size();
             if (n < min_size || n > max_size) 
                 throw IndexError ("Unexpected SeqBase<T> length.");
         }
         
         class iterator
             : public std::iterator<std::random_access_iterator_tag, seqref<T>, int> {
         protected:
             friend class SeqBase<T>;
             SeqBase<T>* seq;
             int count;
             
         public:
             ~iterator () {}
             
             iterator () {
                 seq = 0;
                 count = 0;
             }
             
             iterator (SeqBase<T>* s, int where) 
                 :seq(s), count(where) {} 
             
             iterator (const iterator& other) {
                 seq = other.seq;
                 count = other.count;
             }
             
             bool operator== (const iterator& other) {
                 return (seq == other.seq) && (count == other.count);
             }
             
             bool operator!= (const iterator& other) {
                 return (seq != other.seq) || (count != other.count);
             }
             
             bool operator< (const iterator& other) {
                 return (count < other.count);
             }
             
             bool operator> (const iterator& other) {
                 return (count > other.count);
             }
             
             bool operator<= (const iterator& other) {
                 return (count <= other.count);
             }
             
             bool operator>= (const iterator& other) {
                 return (count >= other.count);
             }
             
             seqref<T> operator*() {
                 return seqref<T>(*seq, count);
             }
             
             seqref<T> operator[] (int i) {
                 return seqref<T>(*seq, count + i);
             }
             
             iterator& operator=(const iterator& other) {
                 if (this == &other) return *this;
                 seq = other.seq;
                 count = other.count;
                 return *this;
             }
             
             iterator operator+(int n) const {
                 return iterator(seq, count + n);
             }
             
             iterator operator-(int n) const {
                 return iterator(seq, count - n);
             }
             
             iterator& operator+=(int n) {
                 count = count + n;
                 return *this;
             }
             
             iterator& operator-=(int n) {
                 count = count - n;
                 return *this;
             }
             
             int operator-(const iterator& other) const {
                 if (seq != other.seq) 
                     throw RuntimeError ("SeqBase<T>::iterator comparison error");
                 return count - other.count;
             }
             
             // prefix ++
             iterator& operator++ () { count++; return *this;}
             // postfix ++
             iterator operator++ (int) { return iterator(seq, count++);}
             // prefix --
             iterator& operator-- () { count--; return *this;}
             // postfix --
             iterator operator-- (int) { return iterator(seq, count--);}
             
             std::string diagnose() const {
                 std::ostrstream oss;
                 oss << "iterator diagnosis " << seq << ", " << count << std::ends;
                 return std::string(oss.str());
             }
         };    // end of class SeqBase<T>::iterator
             
         iterator begin () {
             return iterator(this, 0);
         }
         
         iterator end () {
             return iterator(this, length());
         }
         
         class const_iterator
             : public std::iterator<std::random_access_iterator_tag, const Object,int> {
         protected:
             friend class SeqBase<T>;
             const SeqBase<T>* seq;
             int count;
             
         public:
             ~const_iterator () {}
             
             const_iterator () {
                 seq = 0;
                 count = 0;
             }
             
             const_iterator (const SeqBase<T>* s, int where) 
                 :seq(s), count(where) {} 
             
             const_iterator(const const_iterator& other) {
                 seq = other.seq;
                 count = other.count;
             }
             
             bool operator== (const const_iterator& other) {
                 return (count == other.count);
             }
             
             bool operator!= (const const_iterator& other) {
                 return (count != other.count);
             }
             
             bool operator< (const const_iterator& other) {
                 return (count < other.count);
             }
             
             bool operator> (const const_iterator& other) {
                 return (count > other.count);
             }
             
             bool operator<= (const const_iterator& other) {
                 return (count <= other.count);
             }
             
             bool operator>= (const const_iterator& other) {
                 return (count >= other.count);
             }
             
             const T operator*() const { 
                 return seq->getItem(count);
             }
             
             const T operator[] (int i) const {
                 return seq->getItem(count + i);
             }
             
             const_iterator& operator=(const const_iterator& other) {
                 if (this == &other) return *this;
                 seq = other.seq;
                 count = other.count;
                 return *this;
             }
             
             const_iterator operator+(int n) const {
                 return const_iterator(seq, count + n);
             }
             
             const_iterator operator-(int n) const {
                 return const_iterator(seq, count - n);
             }
             
             const_iterator& operator+=(int n) {
                 count = count + n;
                 return *this;
             }
             
             const_iterator& operator-=(int n) {
                 count = count - n;
                 return *this;
             }
             
             int operator-(const const_iterator& other) const {
                 if (seq != other.seq) 
                     throw RuntimeError ("SeqBase<T>::const_iterator::- error");
                 return count - other.count;
             }
             // prefix ++
             const_iterator& operator++ () { count++; return *this;}
             // postfix ++
             const_iterator operator++ (int) { return const_iterator(seq, count++);}
             // prefix --
             const_iterator& operator-- () { count--; return *this;}
             // postfix --
             const_iterator operator-- (int) { return const_iterator(seq, count--);}
         };    // end of class SeqBase<T>::const_iterator
         
         const_iterator begin () const {
             return const_iterator(this, 0);
         }
         
         const_iterator end () const {
             return const_iterator(this, length());
         }
    };
    
    // Here's an important typedef you might miss if reading too fast...
    typedef SeqBase<Object> Sequence;
    
    // ==================================================
    // class String
    // Python strings return strings as individual elements.
    // I'll try having a class Char which is a String of length 1
    //
    class Char: public Object {
    public:
        explicit Char (PyObject *pyob): Object(pyob) {
            validate();
        }
        
        Char (const Object& ob): Object(ob) {
            validate();
        }
        
        Char (const std::string& v = "")
            :Object(FromAPI(PyString_FromStringAndSize (const_cast<char*>(v.c_str()),1)))
        {
            validate();
        }
        
        Char (char v)
            :Object(FromAPI(PyString_FromStringAndSize (&v, 1))) 
        {
            validate();
        }
        // Assignment acquires new ownership of pointer
        
        Char& operator= (const Object& rhs) {
            return (*this = *rhs);
        }
        
        Char& operator= (PyObject* rhsp) {
            if(ptr() == rhsp) return *this;
            set (rhsp);
            return *this;
        }
        // Membership
        virtual bool accepts (PyObject *pyob) const {
            return pyob && PyString_Check (pyob) && 
                (PySequence_Length (pyob) == 1);
        }
        
        // Assignment from C string
        Char& operator= (const std::string& v) {
            *this = FromAPI(PyString_FromStringAndSize (const_cast<char*>(v.c_str()),1));
            return *this;
        }
        
        Char& operator= (char v) {
            *this = FromAPI(PyString_FromStringAndSize (&v, 1));
            return *this;
        }
        
        // Conversion
        operator String() const;
        
        operator std::string () const {
            return std::string(PyString_AsString (ptr()));
        }           
    };
    
    class String: public SeqBase<Char> {
    public:
        virtual size_type capacity() const {
            return max_size();
        }
        
        explicit String (PyObject *pyob): SeqBase<Char>(pyob) {
            validate();
        }
        
        String (const Object& ob): SeqBase<Char>(ob) {
            validate();
        }
        
        String (const std::string& v = "")
            :SeqBase<Char>(FromAPI(PyString_FromString (const_cast<char*>(v.c_str()))))
        {
            validate();
        }
        
        String (const std::string& v, std::string::size_type vsize)
            :SeqBase<Char>(FromAPI(PyString_FromStringAndSize (const_cast<char*>(v.c_str()), vsize)))
        {
            validate();
        }
        
        String (const char* v)
            :SeqBase<Char>(FromAPI(PyString_FromString (v))) 
        {
            validate();
        }
        // Assignment acquires new ownership of pointer
        
        String& operator= (const Object& rhs) {
            return (*this = *rhs);
        }
        
        String& operator= (PyObject* rhsp) {
            if(ptr() == rhsp) return *this;
            set (rhsp);
            return *this;
        }
        // Membership
        virtual bool accepts (PyObject *pyob) const {
            return pyob && PyString_Check (pyob);
        }
        
        // Assignment from C string
        String& operator= (const std::string& v) {
            *this = FromAPI(PyString_FromString (const_cast<char*>(v.c_str())));
            return *this;
        }
        // Queries
        virtual size_type size () const { 
            return PyString_Size (ptr());
        }
        
        operator std::string () const {
            return std::string(PyString_AsString (ptr()));
        }           
    };
    
    // ==================================================
    // class Tuple
    class Tuple: public Sequence {
    public:
        virtual void setItem (int offset, const Object&ob) {
            // note PyTuple_SetItem is a thief...
            if(PyTuple_SetItem (ptr(), offset, new_reference_to(ob)) == -1) {
                throw Exception();
            }
        }
        
        // Constructor
        explicit Tuple (PyObject *pyob): Sequence (pyob) {
            validate();
        }
        
        Tuple (const Object& ob): Sequence(ob) {
            validate();
        }
        
        // New tuple of a given size	
        explicit Tuple (int size = 0) {
            set(FromAPI(PyTuple_New (size)));
            validate ();
            for (int i=0; i < size; i++) {
                if(PyTuple_SetItem (ptr(), i, new_reference_to(Py_None)) == -1) {
                    throw Exception();
                }
            }
        }
        // Tuple from any sequence
        explicit Tuple (const Sequence& s) { 
            set(FromAPI(PyTuple_New (s.length())));
            validate();
            for(int i=0; i < s.length(); i++) {
                if(PyTuple_SetItem (ptr(), i, new_reference_to(s[i])) == -1) {
                    throw Exception();
                }
            }
        }      
        // Assignment acquires new ownership of pointer
        
        Tuple& operator= (const Object& rhs) {
            return (*this = *rhs);
        }
        
        Tuple& operator= (PyObject* rhsp) {
            if(ptr() == rhsp) return *this;
            set (rhsp);
            return *this;
        }
        // Membership	
        virtual bool accepts (PyObject *pyob) const {
            return pyob && PyTuple_Check (pyob);
        }
        
        Tuple getSlice (int i, int j) const {
            return Tuple (FromAPI(PySequence_GetSlice (ptr(), i, j)));
        }
        
    };
    
    // ==================================================
    // class List
    
    class List: public Sequence { 
    public:
        // Constructor
        explicit List (PyObject *pyob): Sequence(pyob) {
            validate();
        }
        List (const Object& ob): Sequence(ob) {
            validate();
        }
        // Creation at a fixed size	
        List (int size = 0) {
            set(FromAPI(PyList_New (size)));
            validate();
            for (int i=0; i < size; i++) {
                if(PyList_SetItem (ptr(), i, new_reference_to(Py_None)) == -1) {
                    throw Exception();
                }
            }
        }
        
        // List from a sequence
        List (const Sequence& s): Sequence() {
            int n = s.length();
            set(FromAPI(PyList_New (n)));
            validate();
            for (int i=0; i < n; i++) {
                if(PyList_SetItem (ptr(), i, new_reference_to(s[i])) == -1) {
                    throw Exception();
                }
            }
        }
        
        virtual size_type capacity() const {
            return max_size();
        }
        // Assignment acquires new ownership of pointer
        
        List& operator= (const Object& rhs) {
            return (*this = *rhs);
        }
        
        List& operator= (PyObject* rhsp) {
            if(ptr() == rhsp) return *this;
            set (rhsp);
            return *this;
        }
        // Membership
        virtual bool accepts (PyObject *pyob) const {
            return pyob && PyList_Check (pyob);
        }
        
        List getSlice (int i, int j) const {
            return List (FromAPI(PyList_GetSlice (ptr(), i, j)));
        }
        
        void setSlice (int i, int j, const Object& v) {
            if(PyList_SetSlice (ptr(), i, j, *v) == -1) {
                throw Exception();
            }
        }
        
        void append (const Object& ob) {
            if(PyList_Append (ptr(), *ob) == -1) {
                throw Exception();
            }
        }
        
        void insert (int i, const Object& ob) {
            if(PyList_Insert (ptr(), i, *ob) == -1) {
                throw Exception();
            }
        }
        
        void sort () {
            if(PyList_Sort(ptr()) == -1) {
                throw Exception();
            }
        }
        
        void reverse () {
            if(PyList_Reverse(ptr()) == -1) {
                throw Exception();
            }
        }
    };
    // Mappings
    // ==================================================
    template<class T>
        class mapref { 
    protected:
        MapBase<T>& s; // the map
        std::string key; // item number
        T the_item;
        
    public:
        mapref<T> (MapBase<T>& map, const std::string& k)
            : s(map), key(k), the_item() 
        {
            if(map.hasKey(k)) the_item = map.getItem(k);
        }; 
        
        ~mapref<T>() {}
        
        // MapBase<T> stuff
        // lvalue
        mapref<T>& operator=(const mapref<T>& other) {
            if(this == &other) return *this;
            the_item = other.the_item;
            s.setItem(key, other.the_item);
            return *this;
        };
        
        mapref<T>& operator= (const T& ob) {
            the_item = ob;
            s.setItem (key, ob);
            return *this;
        }
        
        // rvalue
        operator T() const {
            return the_item;
        }
        
        // forward everything else to the_item           
        PyObject* ptr () const {
            return the_item.ptr();
        }
        
        int reference_count () const { // the mapref count
            return the_item.reference_count();
        }
        
        Type type () const {
            return the_item.type();
        }
        
        String str () const {
            return the_item.str();
        }
        
        String repr () const {
            return the_item.repr();
        }
        
        bool hasAttr (const std::string& attr_name) const {
            return the_item.hasAttr(attr_name);
        }
        
        Object getAttr (const std::string& attr_name) const {
            return the_item.getAttr(attr_name);
        }
        
        Object getItem (const Object& k) const {
            return the_item.getItem(k);
        }
        
        long hashValue () const {
            return the_item.hashValue();
        }
        
        bool isCallable () const {
            return the_item.isCallable();
        }
        
        bool isList () const {
            return the_item.isList();
        }
        
        bool isMapping () const {
            return the_item.isMapping();
        }
        
        bool isNumeric () const {
            return the_item.isNumeric();
        }
        
        bool isSequence () const {
            return the_item.isSequence();
        }
        
        bool isTrue () const {
            return the_item.isTrue();
        }
        
        bool isType (const Type& t) const {
            return the_item.isType (t);
        }
        
        bool isTuple() const {
            return the_item.isTuple();
        }
        
        bool isString() const {
            return the_item.isString();
        }
        
        // Commands
        void setAttr (const std::string& attr_name, const Object& value) {
            the_item.setAttr(attr_name, value);
        }
        
        void delAttr (const std::string& attr_name) {
            the_item.delAttr(attr_name);
        }
        
        void delItem (const Object& k) {
            the_item.delItem(k);
        }
        }; // end of mapref
        
        template<class T>
            class MapBase: public Object {
    protected:
        explicit MapBase<T>(){}
    public:
        // reference: proxy class for implementing []
        
        // Constructor
        explicit MapBase<T> (PyObject *pyob): Object(pyob) {
            validate();
        }
        
        MapBase<T> (const MapBase<T>& ob): Object(ob) {
            validate();
        }
        // Assignment acquires new ownership of pointer
        
        MapBase<T>& operator= (const Object& rhs) {
            return (*this = *rhs);
        }
        
        MapBase<T>& operator= (PyObject* rhsp) {
            if(ptr() == rhsp) return *this;
            set (rhsp);
            return *this;
        }
        // Membership
        virtual bool accepts (PyObject *pyob) const {
            return pyob && PyMapping_Check(pyob);
        }
        
        // Clear -- PyMapping Clear is missing
        //
        
        void clear () {
            List k = keys();
            for(List::iterator i = k.begin(); i != k.end(); i++) {
                delItem(*i);
            }
        }      
        
        // Element Access
        T operator[](const std::string& key) const {
            return getItem(key);
        }
        
        mapref<T> operator[](const std::string& key) {
            return mapref<T>(*this, key);
        }
        
        int length () const {
            return PyMapping_Length (ptr());
        }
        
        int hasKey (const std::string& s) const {
            return PyMapping_HasKeyString (ptr(),const_cast<char*>(s.c_str()));
        }
        
        T getItem (const std::string& s) const {
            return Object(
                FromAPI(PyMapping_GetItemString (ptr(),const_cast<char*>(s.c_str())))
                );
        }
        
        virtual void setItem (const std::string& s, const Object& ob) {
            if (PyMapping_SetItemString (ptr(), const_cast<char*>(s.c_str()), *ob)
                == -1)
            {
                throw Exception();
            }
        }
        
        void delItem (const std::string& s) {
            if (PyMapping_DelItemString (ptr(), const_cast<char*>(s.c_str())) == -1){
                throw Exception();
            }
        }
        
        void delItem (const Object& s) {
            if (PyMapping_DelItem (ptr(), *s) == -1){
                throw Exception();
            }
        }
        // Queries
        List keys () const {
            return List(FromAPI(PyMapping_Keys(ptr())));
        }
        
        List values () const { // each returned item is a (key, value) pair
            return List(FromAPI(PyMapping_Values(ptr())));
        }
        
        List items () const {
            return List(FromAPI(PyMapping_Items(ptr())));
        }	   
        };
        
        typedef MapBase<Object> Mapping;
        // ==================================================
        // class Dict
        class Dict: public Mapping {
        public:
            // Constructor
            explicit Dict (PyObject *pyob): Mapping (pyob) {
                validate();
            }
            Dict (const Dict& ob): Mapping(ob) {
                validate();
            }
            // Creation
            Dict () {
                set(FromAPI(PyDict_New ()));
                validate();
            }
            // Assignment acquires new ownership of pointer
            
            Dict& operator= (const Object& rhs) {
                return (*this = *rhs);
            }
            
            Dict& operator= (PyObject* rhsp) {
                if(ptr() == rhsp) return *this;
                set(rhsp);
                return *this;
            }
            // Membership
            virtual bool accepts (PyObject *pyob) const {
                return pyob && PyDict_Check (pyob);
            }
        };
        
        class Callable: public Object {
        protected:
            explicit Callable (): Object() {}
        public:
            // Constructor
            explicit Callable (PyObject *pyob): Object (pyob) {
                validate();
            }
            
            Callable (const Object& ob): Object(ob) {
                validate();
            }
            
            // Assignment acquires new ownership of pointer
            
            Callable& operator= (const Object& rhs) {
                return (*this = *rhs);
            }
            
            Callable& operator= (PyObject* rhsp) {
                if(ptr() == rhsp) return *this;
                set (rhsp);
                return *this;
            }
            
            // Membership
            virtual bool accepts (PyObject *pyob) const {
                return pyob && PyCallable_Check (pyob);
            }
            
            // Call
            Object apply(const Tuple& args) const {
                return Object(FromAPI(PyObject_CallObject(ptr(), *args)));
            }
            
            Object apply(PyObject* args = 0) const {
                return Object(FromAPI(PyObject_CallObject(ptr(), args)));
            }
        }; 
        
        class Module: public Object {
        public:
            // Constructors acquire new ownership of pointer
            explicit Module (PyObject* pyob): Object (pyob) {
                validate();
            }
            
            // Copy constructor acquires new ownership of pointer
            Module (const Module& ob): Object(*ob) {
                validate();
            }
            
            
            Module& operator= (const Object& rhs) {
                return (*this = *rhs);
            }
            
            Module& operator= (PyObject* rhsp) {
                if(ptr() == rhsp) return *this;
                set(rhsp);
                return *this;
            }       
        };
        
        // Numeric interface
        inline Object operator+ (const Object& a) {
            return Object(FromAPI(PyNumber_Positive(*a)));
        }
        inline Object operator- (const Object& a) {
            return Object(FromAPI(PyNumber_Negative(*a)));
        }
        
        inline Object abs(const Object& a) {
            return Object(FromAPI(PyNumber_Absolute(*a)));
        }
        
        inline std::pair<Object,Object> coerce(const Object& a, const Object& b) {
            PyObject *p1, *p2;
            p1 = *a;
            p2 = *b;
            if(PyNumber_Coerce(&p1,&p2) == -1) {
                throw Exception();
            }
            return std::pair<Object,Object>(Object(FromAPI(p1)), Object(FromAPI(p2)));
        }
        
        inline Object operator+ (const Object& a, const Object& b) {
            return Object(FromAPI(PyNumber_Add(*a, *b)));
        }
        inline Object operator+ (const Object& a, int j) {
            return Object(FromAPI(PyNumber_Add(*a, *Int(j))));
        }
        inline Object operator+ (const Object& a, double v) {
            return Object(FromAPI(PyNumber_Add(*a, *Float(v))));
        }
        inline Object operator+ (int j, const Object& b) {
            return Object(FromAPI(PyNumber_Add(*Int(j), *b)));
        }
        inline Object operator+ (double v, const Object& b) {
            return Object(FromAPI(PyNumber_Add(*Float(v), *b)));
        }
        
        inline Object operator- (const Object& a, const Object& b) {
            return Object(FromAPI(PyNumber_Subtract(*a, *b)));
        }
        inline Object operator- (const Object& a, int j) {
            return Object(FromAPI(PyNumber_Subtract(*a, *Int(j))));
        }
        inline Object operator- (const Object& a, double v) {
            return Object(FromAPI(PyNumber_Subtract(*a, *Float(v))));
        }
        inline Object operator- (int j, const Object& b) {
            return Object(FromAPI(PyNumber_Subtract(*Int(j), *b)));
        }
        inline Object operator- (double v, const Object& b) {
            return Object(FromAPI(PyNumber_Subtract(*Float(v), *b)));
        }
        
        inline Object operator* (const Object& a, const Object& b) {
            return Object(FromAPI(PyNumber_Multiply(*a, *b)));
        }
        inline Object operator* (const Object& a, int j) {
            return Object(FromAPI(PyNumber_Multiply(*a, *Int(j))));
        }
        inline Object operator* (const Object& a, double v) {
            return Object(FromAPI(PyNumber_Multiply(*a, *Float(v))));
        }
        inline Object operator* (int j, const Object& b) {
            return Object(FromAPI(PyNumber_Multiply(*Int(j), *b)));
        }
        inline Object operator* (double v, const Object& b) {
            return Object(FromAPI(PyNumber_Multiply(*Float(v), *b)));
        }
        
        inline Object operator/ (const Object& a, const Object& b) {
            return Object(FromAPI(PyNumber_Divide(*a, *b)));
        }
        inline Object operator/ (const Object& a, int j) {
            return Object(FromAPI(PyNumber_Divide(*a, *Int(j))));
        }
        inline Object operator/ (const Object& a, double v) {
            return Object(FromAPI(PyNumber_Divide(*a, *Float(v))));
        }
        inline Object operator/ (int j, const Object& b) {
            return Object(FromAPI(PyNumber_Divide(*Int(j), *b)));
        }
        inline Object operator/ (double v, const Object& b) {
            return Object(FromAPI(PyNumber_Divide(*Float(v), *b)));
        }
        
        inline Object operator% (const Object& a, const Object& b) {
            return Object(FromAPI(PyNumber_Remainder(*a, *b)));
        }
        inline Object operator% (const Object& a, int j) {
            return Object(FromAPI(PyNumber_Remainder(*a, *Int(j))));
        }
        inline Object operator% (const Object& a, double v) {
            return Object(FromAPI(PyNumber_Remainder(*a, *Float(v))));
        }
        inline Object operator% (int j, const Object& b) {
            return Object(FromAPI(PyNumber_Remainder(*Int(j), *b)));
        }
        inline Object operator% (double v, const Object& b) {
            return Object(FromAPI(PyNumber_Remainder(*Float(v), *b)));
        }
        
} // end namespace Py
#endif	// __CXX_Objects__h	
