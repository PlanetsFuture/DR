#ifndef _DR_H
#define _DR_H

// system #include
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string>
#include <ctype.h>

// macro-like inline functions
template<class T>
inline T SQR(const T a)
{
    return a * a;
}

template<class T>
inline const T &MAX(const T &a, const T &b)
{
    return b > a ? (b) : (a);
}

inline float MAX(const double &a, const float &b)
{
    return b > a ? (b) : float (a);
}

inline float MAX(const float &a, const double &b)
{
    return b > a ? float (b) : (a);
}

template<class T>
inline const T &MIN(const T &a, const T &b)
{
    return b < a ? (b) : (a);
}

inline float MIN(const double &a, const float &b)
{
    return b < a ? (b) : float (a);
}

inline float MIN(const float &a, const double &b)
{
    return b < a ? float (b) : (a);
}

template<class T>
inline T SIGN(const T &a, const T &b)
{
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float SIGN(const float &a, const double &b)
{
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float SIGN(const double &a, const float &b)
{
    return (float) (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

template<class T>
inline void SWAP(T &a, T &b)
{
    T dum = a;
    a = b;
    b = dum;
}

// exception handling
#ifndef _USEDRERRORCLASS_
#define throw(message)                              \
    {printf("ERROR: %s\n in file %s at line %d\n",  \
            message,__FILE__,__LINE__); \
        throw(1);}
#else
struct DRerror {
	char *message;
	char *file;
	int line;
    DRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};

#define throw(message) throw(DRerror(message,__FILE__,__LINE__));
void DRcatch(DRerror err) {
	printf("ERROR: %s\n     in file %s at line %d\n",
           err.message, err.file, err.line);
	exit(1);
}
#endif

// Vector und Matrix Classes

#ifdef _USESTDVECTOR_
#define DRvector vector
#else

template<class T>
class DRvector {
private:
    int nn; // size of array. upper index is nn - 1
    T *v;
public:
    DRvector();
    explicit DRvector(int n); // Zero-based array
    DRvector(int n, const T &a); // initialize to constant value
    DRvector(int n, const T *a); // initialize to array
    DRvector(const DRvector &rhs); // copy constructor
    DRvector & operator = (const DRvector &rhs); // assignment
    typedef T value_type; // make T available externally
    inline T & operator[] (const int i); // i'th element
    inline const T & operator[] (const int i) const;
    inline int size() const;
    void resize(int newn); // resize (contents not preserved)
    void assign(int newn, const T &a);
    ~DRvector();
};

// DRvector definitions

template<class T>
DRvector<T>::DRvector() : nn(0), v(NULL) {}

template<class T>
DRvector<T>::DRvector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}

template<class T>
DRvector<T>::DRvector(int n, const T& a) : nn(n), v(n > 0 ? new T[n] : NULL )
{
    for (int i = 0; i < n; i++) v[i] = a;
}

template<class T>
DRvector<T>::DRvector(int n, const T *a) : nn(n), v(n > 0 ? new T[n] : NULL)
{
    for (int i = 0; i < n; i++) v[i] = *a++;
}

template <class T>
DRvector<T>::DRvector(const DRvector<T> &rhs) : nn(rhs.nn),
    v(nn > 0 ? new T[nn] : NULL)
{
    for (int i = 0; i < nn; i++) v[i] = rhs[i];
}

template<class T>
DRvector<T> & DRvector<T>::operator = (const DRvector<T> &rhs)
    /**
     * @postcondition: normal assignment via copying has been performe;
     * if vector and rhs were different sizes, vector
     * has been resized to match the size of rhs.
     */
{
    if (this != rhs)
    {
        if (nn != rhs.nn) {
            if (v != NULL) delete [] (v);
            nn = rhs.nn;
            v = nn > 0 ? new T[nn] : NULL;
        }
        for (int i = 0; i < nn; i++)
            v[i] = rhs[i];
    }
}

template<class T>
inline T & DRvector<T>::operator[] (const int i) // subscripting
{
#ifdef _CHECKBOUNDS_
    if (i < 0 || i >= nn) {
        throw("DRvector subscript out of bounds.");
    }
#endif
    return v[i];
}

template<class T>
inline const T & DRvector<T>::operator[] (const int i) const // subscripting
{
#ifdef _CHECKBOUNDS_
    if (i < 0 || i >= nn) {
        throw("DRvector subscript out of bounds.");
    }
#endif
    return v[i];
}

template<class T>
inline int DRvector<T>::size() const
{
    return nn;
}

template<class T>
void DRvector<T>::resize(int newn)
{
    if (newn != nn) {
        if (v != NULL) delete[] (v);
        nn = newn;
        v = nn > 0 ? new T[nn] : NULL;
    }
}

template<class T>
void DRvector<T>::assign(int newn, const T& a)
{
    if (newn != nn) {
        if (v != NULL) delete[] (v);
        nn = newn;
        v = nn > 0 ? new T[nn] : NULL;
    }
    for (int i = 0; i < nn; i++) v[i] = a;
}

template <class T>
DRvector<T>::~DRvector()
{
    if (v != NULL) delete[] (v);
}

// end of DRvector definitions
#endif // ifdef _USESTDVECTOR_

template<class T>
class DRmatrix {
private:
    int nn;
    int mm;
    T **v;
public:
    DRmatrix();
    DRmatrix(int n, int m); // Zero-based array
    DRmatrix(int n, int m, const T &a); // Initialized to constant
    DRmatrix(int n, int m, const T *a); // Initialize to array
    DRmatrix(const DRmatrix &rhs); // Copy constructor
    DRmatrix & operator = (const DRmatrix &rhs);
    typedef T value_type; // Make T available externally
    inline T* operator[] (const int i); // pointer to row i
    inline const T* operator[] (const int i) const;
    inline int nrows() const;
    inline int ncols() const;
    void resize(int newn, int newm); // resize (contents not preserved)
    void assign(int newn, int newm, const T &a); // resize and assign value
    ~DRmatrix();
};

template<class T>
DRmatrix<T>::DRmatrix() : nn(0), mm(0), v(NULL) {}

template<class T>
DRmatrix<T>::DRmatrix(int n, int m) : nn(n), mm(m), v(n > 0 ? new T*[n] : NULL)
{
    int i, nel = m * n;
    if (v) v[0] = nel > 0 ? new T[nel] : NULL;
    for (i = 1; i < n; i++) v[i] = v[i-1] + m;
}

template<class T>
DRmatrix<T>::DRmatrix(int n, int m, const T &a) : nn(n), mm(m),
    v(n > 0 ? new T*[n] : NULL)
{
    int i, j, nel = m * n;
    if (v) v[0] = nel > 0 ? new T[nel] : NULL;
    for (i = 0; i < n; i++) v[i] = v[i-1] + m;
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++) v[i][j] = a;
}

template<class T>
DRmatrix<T>::DRmatrix(int n, int m, const T *a) : nn(n), mm(m),
    v(n > 0 ? new T*[n] : NULL)
{
    int i, j, nel = m * n;
    if (v) v[0] = nel > 0 ? new T[nel] : NULL;
    for (i = 1; i < n; i++) v[i] = v[i-1] + m;
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++) v[i][j] = *a++;
}

template<class T>
DRmatrix<T>::DRmatrix(const DRmatrix &rhs) : nn(rhs.nn), mm(rhs.mm),
    v(nn > 0 ? new T*[nn] : NULL)
{
    int i, j, nel = mm * nn;
    if (v) v[0] = nel > 0 ? new T[nel] : NULL;
    for (i = 1; i < nn; i++) v[i] =v[i-1] + mm;
    for (i = 0; i < nn; i++)
        for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
}

template<class T>
DRmatrix<T> & DRmatrix<T>::operator = (const DRmatrix<T> &rhs)
    /**
     * @postcondition: normal assignment via copying has been performed;
     * if matrix and rhs were different sizes, matrix
     * has been resized to match the size of rhs
     */
{
    if (this != &rhs) {
        int i, j, nel;
        if (nn != rhs.nn || mm != rhs.mm) {
            if (v != NULL) {
                delete[] (v[0]);
                delete[] (v);
            }
            nn = rhs.nn;
            mm = rhs.mm;
            v = nn > 0 ? new T*[nn] : NULL;
            nel = mm * nn;
            if (v) v[0] = nel > 0 ? new T[nel] : NULL;
            for (i = 1; i < nn; i++) v[i] = v[i-1] + mm;
        }
        for (i = 0; i < nn; i++)
            for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
    }
    return *this;
}

template<class T>
inline T* DRmatrix<T>::operator[] (const int i) // subscripting
{
#ifdef _CHECKBOUNDS_
    if (i < 0 || i >= nn) {
        throw("DRmatrix subscript out of bounds");
    }
#endif
    return v[i];
}

template<class T>
inline const T* DRmatrix<T>::operator[] (const int i) const
{
#ifdef _CHECKBOUNDS_
    if(i < 0 || i >= nn) {
        throw("DRmatrix subscript out of bounds");
    }
#endif
    return v[i];
}

template<class T>
inline int DRmatrix<T>::nrows() const
{
    return nn;
}

template<class T>
inline int DRmatrix<T>::ncols() const
{
    return mm;
}

template<class T>
void DRmatrix<T>::resize(int newn, int newm)
{
    int i, nel;
    if (newn != nn || newm != mm) {
        if (v != NULL) {
            delete[] (v[0]);
            delete[] (v);
        }
        nn = newn;
        mm = newm;
        v = nn > 0 ? new T[nel] : NULL;
        for (i = 1; i < nn; i++) v[i] =v[i-1] + mm;
    }
}

template<class T>
void DRmatrix<T>::assign(int newn, int newm, const T& a)
{
    int i, j, nel;
    if (newn != nn || newm != mm) {
        if (v != NULL) {
            delete[] (v[0]);
            delete[] (v);
        }
        nn = newn;
        mm = newm;
        v = nn > 0 ? new T*[nn] : NULL;
        nel = mm * nn;
        if (v) v[0] = nel > 0 ? new T[nel] : NULL;
        for (i = 1; i < nn; i++) v[i] = v[i-1] + mm;
    }
    for (i = 0; i < nn; i++)
        for (j = 0; j < mm; j++) v[i][j] = a;
}

template<class T>
DRmatrix<T>::~DRmatrix()
{
    if (v != NULL) {
        delete[] (v[0]);
        delete[] (v);
    }
}

template<class T>
class DRMatrix3D {
private:
    int nn;
    int mm;
    int kk;
    T ***v;
public:
    DRMatrix3D();
    DRMatrix3D(int n, int m, int k);
    inline T** operator[] (const int i); // subscripting: pointer to row i
    inline const T* const * operator[] (const int i) const;
    inline int dim1() const;
    inline int dim2() const;
    inline int dim3() const;
    ~DRMatrix3D();
};

template<class T>
DRMatrix3D<T>::DRMatrix3D(): nn(0), mm(0), kk(0), v(NULL) {}

template<class T>
DRMatrix3D<T>::DRMatrix3D(int n, int m, int k) : nn(n), mm(m), kk(k),
    v(new T**[n])
{
    int i, j;
    v[0] = new T*[n * m];
    v[0][0] = new T[n * m * k];
    for (j = 1; j < m; j++) v[0][j] = v[0][j-1] + k;
    for (i = 1; i < n; i++) {
        v[i] = v[i-1] + m;
        v[i][0] = v[i-1][0] + m * k;
        for (j = 1; j < m; j++) v[i][j] = v[i][j-1] + k;
    }
}

template<class T>
inline T** DRMatrix3D<T>::operator[] (const int i) // subscripting:
{
    return v[i];
}

template<class T>
inline const T* const * DRMatrix3D<T>::operator[] (const int i) const
{
    return v[i];
}

template<class T>
inline int DRMatrix3D<T>::dim1() const
{
    return nn;
}

template<class T>
inline int DRMatrix3D<T>::dim2() const
{
    return mm;
}

template<class T>
inline int DRMatrix3D<T>::dim3() const
{
    return kk;
}

template<class T>
DRMatrix3D<T>::~DRMatrix3D()
{
    if (v != NULL) {
        delete[] (v[0][0]);
        delete[] (v[0]);
        delete[] (v);
    }
}

// Vector types
namespace DR {
    typedef const DRvector<int> vector_int_I;
    typedef DRvector<int> vector_int, vector_int_O,
        vector_int_IO;

    typedef const DRvector<unsigned int> vector_uint_I;
    typedef DRvector<unsigned int> vector_uint, vector_uint_O,
        vector_uint_IO;

    typedef const DRvector<long long int> vector_llong_int_I;
    typedef DRvector<long long int> vector_llong_int,
        vector_llong_int_O, vector_llong_int_IO;

    typedef const DRvector<unsigned long long int> vector_ullong_int_I;
    typedef DRvector<unsigned long long int> vector_ullong_int,
        vector_ullong_int_O, vector_ullong_int_IO;

    typedef const DRvector<char> vector_char_I;
    typedef DRvector<char> vector_char, vector_char_O,
        vector_char_IO;

    typedef const DRvector<char*> vector_char_p_I;
    typedef DRvector<char*> vector_char_p,
        vector_char_p_O, vector_char_p_IO;

    typedef const DRvector<unsigned char> vector_uchar_I;
    typedef DRvector<unsigned char> vector_uchar,
        vector_uchar_O, vector_uchar_IO;

    typedef const DRvector<double> vector_double_I;
    typedef DRvector<double> vector_double, vector_double_O,
        vector_double_IO;

    typedef const DRvector<double*> vector_p_I;
    typedef DRvector<double*> vector_double_p, vector_double_p_O,
        vector_double_p_IO;

    typedef const DRvector<std::complex<double> > vector_complex_I;
    typedef DRvector<std::complex<double> > vector_complex,
        vector_complex_O, vector_complex_IO;

    typedef const DRvector<bool> vector_bool_I;
    typedef DRvector<bool> vector_bool, vector_bool_O,
        vector_bool_IO;


// matrix types
    typedef const DRmatrix<int> matrix_int_I;
    typedef DRmatrix<int> matrix_int, matrix_int_O, matrix_int_IO;

    typedef const DRmatrix<unsigned int> matrix_uint_I;
    typedef DRmatrix<unsigned int> matrix_uint,
        matrix_uint_O, matrix_uint_IO;

    typedef const DRmatrix<long long int> matrix_llong_int_I;
    typedef DRmatrix<long long int> matrix_llong_int,
        matrix_llong_int_O, matrix_llong_int_IO;

    typedef const DRmatrix<unsigned long long int> matrix_ullong_int_I;
    typedef DRmatrix<unsigned long long int> matrix_ullong_int,
        matrix_ullong_int_O, matrix_ullong_int_IO;

    typedef const DRmatrix<char> matrix_char_I;
    typedef DRmatrix<char> matrix_char, matrix_char_O,
        matrix_char_IO;

    typedef const DRmatrix<unsigned char> matrix_uchar_I;
    typedef DRmatrix<unsigned char> matrix_uchar,
        matrix_uchar_O, matrix_uchar_IO;

    typedef const DRmatrix<double> matrix_double_I;
    typedef DRmatrix<double> matrix_double, matrix_double_O,
        matrix_double_IO;

    typedef const DRmatrix<bool> matrix_bool_I;
    typedef DRmatrix<bool> matrix_bool, matrix_bool_O,
        matrix_bool_IO;

// 3D matrix types

    typedef const DRMatrix3D<double> matrix_3d_double_I;
    typedef DRMatrix3D<double> matrix_3d_double, matrix_3d_double_O,
        matrix_3d_double_IO;
}

#endif /* _H */
