//----------------------------------*-C++-*----------------------------------//
// Copyright 1998 The Regents of the University of California. 
// All rights reserved. See LEGAL.LLNL for full text and disclaimer.
//---------------------------------------------------------------------------//

#ifndef __CXX_Extensions__h
#define __CXX_Extensions__h
#include "CXX_Objects.h"

extern "C" {
    extern PyObject py_object_initializer;
}

#include <vector>

namespace Py {

    class MethodTable {
    private:
        MethodTable(const MethodTable& m); //unimplemented
        void operator=(const MethodTable& m); //unimplemented

    protected:
        std::vector<PyMethodDef> t; // accumulator of PyMethodDef's
        PyMethodDef *mt; // Actual method table produced when full
        
        static PyMethodDef method (
            const char* method_name, 
            PyCFunction f, 
            int flags = 1,
            const char* doc="") 
        {
            PyMethodDef m;
            m.ml_name = const_cast<char*>(method_name);
            m.ml_meth = f;
            m.ml_flags = flags;
            m.ml_doc = const_cast<char*>(doc);
            return m;
        }
        
        
    public:
        MethodTable() {
            t.push_back(method(0,0));
            mt = 0;
        }
        
        virtual ~MethodTable() {
            delete [] mt;
        }
        
        void add(const char* method_name, PyCFunction f, const char* doc="", int flag=1) {
            if (!mt) {
                t.insert (t.end()-1, method(method_name, f, flag, doc));
            }
            else {
                throw RuntimeError("Too late to add a module method!");
            }
        }
        
        PyMethodDef* table() {    
            if (!mt) {
                int t1size = t.size();
                mt = new PyMethodDef[t1size];
                int j = 0;
                for(std::vector<PyMethodDef>::iterator i = t.begin(); i != t.end(); i++) {
                    mt[j++] = *i;
                }
            }
            return mt;
        }
    };
    // Class to build a module
    class ExtensionModule {
    private:
        char* module_name;
        MethodTable* method_table;
        ExtensionModule(const ExtensionModule&); //unimplemented
        void operator=(const ExtensionModule&); //unimplemented
        
    public:
        ExtensionModule (char* name) {
            method_table = new MethodTable;
            module_name = name;
        }
        
        virtual ~ExtensionModule () {
            delete method_table;
        }
        
        void add(const char* method_name, PyCFunction f, const char* doc="") {
            method_table->add (method_name, f, doc);
        }
        
        // Initialize returns the new module dictionary so you can add to it.
        Dict initialize () {
            // Both Py_InitModule and PyModule_GetDict return borrowed refs
            PyObject* pm = Py_InitModule(module_name, method_table->table());
            return Dict(PyModule_GetDict(pm));
        }
        
    };

    
    template<class R>
        class PythonType {
    private:
        PythonType (const PythonType& tb) {};
        void operator=(const PythonType& t) {};
    protected:
        PyTypeObject* table;
        PySequenceMethods* sequence_table;
        PyMappingMethods* mapping_table;
        PyNumberMethods* number_table;
        PyBufferProcs* buffer_table;
        // every object needs a dealloc
        static void standard_dealloc(PyObject* p) {
	        PyMem_DEL(p);
        }

        // use to initialize methods that must be present
        static void missing_method_inquiry(PyObject*) {
            throw AttributeError("Extension object incompletely defined.");
        }
       
        static void missing_method_unary(PyObject*) {
            throw AttributeError("Extension object incompletely defined.");
        }

        static void missing_method_binary(PyObject*, PyObject*) {
            throw AttributeError("Extension object incompletely defined.");
        }

        static void missing_method_ternary(PyObject*,PyObject*,PyObject*) {
            throw AttributeError("Extension object incompletely defined.");
        }

        static void missing_method_intargfunc(PyObject*,int) {
            throw AttributeError("Extension object incompletely defined.");
        }

        static void missing_method_intintargfunc(PyObject*,int,int) {
            throw AttributeError("Extension object incompletely defined.");
        }

        static void missing_method_coercion(PyObject**,PyObject**) {
            throw AttributeError("Extension object incompletely defined.");
        }

        static void missing_method_getreadbufferproc(PyObject*, int, void**) {
            throw AttributeError("Extension object incompletely defined.");
        }

        static void missing_method_getwritebufferproc(PyObject*,int, void**) {
            throw AttributeError("Extension object incompletely defined.");
        }

        static void missing_method_getsegcountproc(PyObject*,int*) {
            throw AttributeError("Extension object incompletely defined.");
        }

       void init_sequence() {
            if(!sequence_table) {
                sequence_table = new PySequenceMethods;
                table->tp_as_sequence = sequence_table;
                sequence_table->sq_length = (inquiry) missing_method_inquiry;
                sequence_table->sq_concat = (binaryfunc) missing_method_binary;
                sequence_table->sq_repeat = (intargfunc) missing_method_intargfunc;
                sequence_table->sq_item = (intargfunc) missing_method_intargfunc;
                sequence_table->sq_slice = (intintargfunc) missing_method_intintargfunc;
                sequence_table->sq_ass_item = 0;
                sequence_table->sq_ass_slice = 0;
            }
        }
        
        void init_mapping() {
            if(!mapping_table) {
                mapping_table = new PyMappingMethods;
                table->tp_as_mapping = mapping_table;
                mapping_table->mp_length = (inquiry) missing_method_inquiry;
                mapping_table->mp_subscript = (binaryfunc) missing_method_binary;
                mapping_table->mp_ass_subscript = 0;
            }
        }

        void init_number() {
            if(!number_table) {
                number_table = new PyNumberMethods;
                table->tp_as_number = number_table;
                number_table->nb_add = (binaryfunc) missing_method_binary;
                number_table->nb_subtract = (binaryfunc) missing_method_binary;
                number_table->nb_multiply = (binaryfunc) missing_method_binary;
                number_table->nb_divide = (binaryfunc) missing_method_binary;
                number_table->nb_remainder = (binaryfunc) missing_method_binary;
                number_table->nb_divmod = (binaryfunc) missing_method_binary;
                number_table->nb_power = (ternaryfunc) missing_method_ternary;
                number_table->nb_negative = (unaryfunc) missing_method_unary;
                number_table->nb_positive = (unaryfunc) missing_method_unary;
                number_table->nb_absolute = (unaryfunc) missing_method_unary;
                number_table->nb_nonzero = (inquiry) missing_method_inquiry;
                number_table->nb_invert = (unaryfunc) missing_method_unary;
                number_table->nb_lshift = 0;
                number_table->nb_rshift = 0;
                number_table->nb_and = (binaryfunc) missing_method_binary;
                number_table->nb_xor = (binaryfunc) missing_method_binary;
                number_table->nb_or = (binaryfunc) missing_method_binary;
                number_table->nb_coerce = (coercion) missing_method_coercion;
                number_table->nb_int = (unaryfunc) missing_method_unary;
                number_table->nb_long = (unaryfunc) missing_method_unary;
                number_table->nb_float = (unaryfunc) missing_method_unary;
                number_table->nb_oct = (unaryfunc) missing_method_unary;
                number_table->nb_hex = (unaryfunc) missing_method_unary;
            }
        }

        void init_buffer() {
            if(!buffer_table) {
                buffer_table = new PyBufferProcs;
                table->tp_as_buffer = buffer_table;
                buffer_table->bf_getreadbuffer = (getreadbufferproc) missing_method_getreadbufferproc;
                buffer_table->bf_getwritebuffer = (getwritebufferproc) missing_method_getwritebufferproc;
                buffer_table->bf_getsegcount = (getsegcountproc) missing_method_getsegcountproc;
           }
        }
    public:
        // if you define one sequence method you must define 
        // all of them except the assigns

        PythonType (int itemsize = 0) {
            number_table = 0;
            sequence_table = 0;
            mapping_table = 0;
            buffer_table = 0;

            table = new PyTypeObject;
            *reinterpret_cast<PyObject*>(table) = py_object_initializer;
            table->ob_type = &PyType_Type;
            table->ob_size = 0;
            table->tp_name = "unknown";
            table->tp_basicsize = sizeof(R);
            table->tp_itemsize = itemsize;
            table->tp_dealloc = (destructor) standard_dealloc;
            table->tp_print = 0;
            table->tp_getattr = 0;
            table->tp_setattr = 0;
            table->tp_compare = 0;
            table->tp_repr = 0;
            table->tp_as_number = 0;
            table->tp_as_sequence = 0;
            table->tp_as_mapping =  0;
            table->tp_hash = 0;
            table->tp_call = 0;
            table->tp_str = 0;
            table->tp_getattro = 0;
            table->tp_setattro = 0;
            table->tp_as_buffer = 0;
            table->tp_xxx4 = 0L;
            table->tp_doc = 0;
#ifdef COUNT_ALLOCS
            table->tp_alloc = 0;
            table->tp_free = 0;
            table->tp_maxalloc = 0;
            table->tp_next = 0;
#endif
        }
        virtual ~PythonType (){
            delete table;
            delete sequence_table;
            delete mapping_table;
            delete number_table;
            delete buffer_table;
        };
        
        PyTypeObject* type_object () const {return table;}
        
        void name (const char* nam) {
            table->tp_name = const_cast<char *>(nam);
        }

        void doc (const char* d) {
            table->tp_doc = const_cast<char *>(d);
        }
        
        void dealloc(void (*f)(R*)) {
            table->tp_dealloc = (destructor) f;
        }
        
        void print (int (*f)(R*, FILE *, int)) {
            table->tp_print = (printfunc) f;
        }
        
        void getattr (PyObject* (*f)(R*, char*)) {
            table->tp_getattr = (getattrfunc) f;
        }
        
        void setattr (int (*f)(R*, char*, PyObject*)) {
            table->tp_setattr = (setattrfunc) f;
        }
        
        void getattro (PyObject* (*f)(R*, PyObject*)) {
            table->tp_getattro = (getattrofunc) f;
        }
        
        void setattro (int (*f)(R*, PyObject*, PyObject*)) {
            table->tp_setattro = (setattrofunc) f;
        }
        
        void compare (int (*f)(R*, PyObject*)) {
            table->tp_compare = (cmpfunc) f;
        }
        
        void repr (PyObject* (*f)(R*)) {
            table->tp_repr = (reprfunc) f;
        }
        
        void str (PyObject* (*f)(R*)) {
            table->tp_str = (reprfunc) f;
        }
        
        void hash (long (*f)(R*)) {
            table->tp_hash = (hashfunc) f;
        }
        
        void call (PyObject* (*f)(PyObject*, PyObject*, PyObject*)) {
            table->tp_call = (ternaryfunc) f;
        }
        
        // Sequence methods
        void sequence_length(int (*f)(R*)) {
            init_sequence();
            sequence_table->sq_length = (inquiry) f;
        }
        
        void sequence_concat(PyObject* (*f)(R*,PyObject*)) {
            init_sequence();
            sequence_table->sq_concat = (binaryfunc) f;
        }
        
        void sequence_repeat(PyObject* (*f)(R*, int)) {
            init_sequence();
            sequence_table->sq_repeat = (intargfunc) f;
        }
        
        void sequence_item(PyObject* (*f)(R*, int)) {
            init_sequence();
            sequence_table->sq_item = (intargfunc) f;
        }
        
        void sequence_slice(PyObject* (*f)(R*, int, int)) {
            init_sequence();
            sequence_table->sq_slice = (intintargfunc) f;
        }
        
        void sequence_ass_item(int (*f)(R*, int, PyObject*)) {
            init_sequence();
            sequence_table->sq_ass_item = (intobjargproc) f;
        }
        
        void sequence_ass_slice(int (*f)(R*, int, int, PyObject*)) {
            init_sequence();
            sequence_table->sq_ass_slice = (intintobjargproc) f;
        }
        // Mapping
        void mapping_length(int (*f)(R*)) {
            init_mapping();
            mapping_table->mp_length = (inquiry) f;
        }
        
        void mapping_subscript(PyObject* (*f)(R*, PyObject*)) {
            init_mapping();
            mapping_table->mp_subscript = (binaryfunc) f;
        }
        
        void mapping_ass_subscript(int (*f)(R*, PyObject*, PyObject*)) {
            init_mapping();
            mapping_table->mp_ass_subscript = (objobjargproc) f;
        }

        // Number
        void number_nonzero (int (*f)(R*)) {
            init_number();
            number_table->nb_nonzero = (inquiry) f;
        }

        void number_coerce (int (*f)(R**, PyObject**)) {
            init_number();
            number_table->nb_coerce = (coercion) f;
        }

        void number_negative (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_negative = (unaryfunc) f;
        }
        void number_positive (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_positive = (unaryfunc) f;
        }
        void number_absolute (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_absolute = (unaryfunc) f;
        }
        void number_invert (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_invert = (unaryfunc) f;
        }
        void number_int (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_int = (unaryfunc) f;
        }
        void number_float (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_float = (unaryfunc) f;
        }
        void number_long (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_long = (unaryfunc) f;
        }
        void number_oct (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_oct = (unaryfunc) f;
        }
        void number_hex (PyObject* (*f)(R*)) {
            init_number();
            number_table->nb_hex = (unaryfunc) f;
        }

        void number_add (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_add = (binaryfunc) f;
        }
        void number_subtract (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_subtract = (binaryfunc) f;
        }
        void number_multiply (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_multiply = (binaryfunc) f;
        }
        void number_divide (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_divide = (binaryfunc) f;
        }
        void number_remainder (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_remainder = (binaryfunc) f;
        }
        void number_divmod (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_divmod = (binaryfunc) f;
        }
        void number_lshift (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_lshift = (binaryfunc) f;
        }
        void number_rshift (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_rshift = (binaryfunc) f;
        }
        void number_and (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_and = (binaryfunc) f;
        }
        void number_xor (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_xor = (binaryfunc) f;
        }
        void number_or (PyObject* (*f)(R*, PyObject*)) {
            init_number();
            number_table->nb_or = (binaryfunc) f;
        }

        void number_power(PyObject* (*f)(R*, PyObject*, PyObject*)) {
            init_number();
            number_table->nb_power = (ternaryfunc) f;
        }
        // Buffer
        void buffer_getreadbuffer (int (*f)(R*, int, void**)) {
            init_buffer();
            buffer_table->bf_getreadbuffer = (getreadbufferproc) f;
        }

        void buffer_getwritebuffer (int (*f)(R*, int, void**)) {
            init_buffer();
            buffer_table->bf_getwritebuffer = (getwritebufferproc) f;
        }

        void buffer_getsegcount (int (*f)(R*, int*)) {
            init_buffer();
            buffer_table->bf_getsegcount = (getsegcountproc) f;
        }
        
    }; // end of PythonType<R>

    template<class T>
        class PythonExtension: public PyObject {
    public:
        static PythonType<T>& behaviors() {
            static PythonType<T>* p;
            if(!p) p = new PythonType<T>();
            return *p;
        }

        static MethodTable& methods() {
            static MethodTable* p;
            if(!p) p = new MethodTable;
            return *p;
        }

        static PyTypeObject* type_object() {
            return behaviors().type_object();
        }

        static void* operator new(size_t t) {
            return static_cast<void*> (PyObject_NEW(T, type_object()));
        }

        static int check (PyObject *p) {
            // is p like me?
            return p->ob_type == type_object();
        }
        
        static int check (const Object& ob) {
            return check(ob.ptr());
        }

        PyObject* getattr(char* name) {
            return Py_FindMethod(methods().table(), static_cast<PyObject*>(this), name);
        }

    };

    // ExtensionObject<r> is an Object that will accept only r's.
    template<class T>
        class ExtensionObject: public Object {
    public:
        
        explicit ExtensionObject (PyObject *pyob): Object(pyob) {
            validate();
        }
        
        ExtensionObject(const ExtensionObject<T>& other): Object(*other) {
            validate();
        }
        
        ExtensionObject& operator= (const Object& rhs) {
            return (*this = *rhs);
        }
        
        ExtensionObject& operator= (PyObject* rhsp) {
            if(ptr() == rhsp) return *this;
            set(rhsp);
            return *this;
        }
        
        virtual bool accepts (PyObject *pyob) const {
            return (pyob && T::check(pyob));
        }       
    };
 
}// namespace Py
// End of CXX_Extensions.h
#endif

