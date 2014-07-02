#pragma once

#include <fstream>

template< typename A, typename T >
class _Array3_get_index_2
{
public:
    _Array3_get_index_2( A * array, int i0, int i1 )
        : array( array ), i0( i0 ), i1( i1 )
    {}

    T &operator[]( int i2 )
    {
        return (*this->array)( this->i0, this->i1, i2 );
    }

private:
    A * array;
    const int i0, i1;
};

template< typename A, typename T >
class _Array3_get_index_1
{
public:
    _Array3_get_index_1( A * array, int i0 )
        : array( array ), i0( i0 )
    {}

    _Array3_get_index_2< A, T > operator[]( int i1 )
    {
        return _Array3_get_index_2< A, T >( this->array, this->i0, i1 );
    }

private:
    A * array;
    const int i0;
};

template< typename T >
class Array3
{
public:
    Array3()
    {}

    Array3( int n0, int n1, int n2, T * data )
        : n0( n0 ), n1( n1 ), n2( n2 ), data( data )
    {}

    void create( int n0, int n1, int n2 )
    {
        this->n0 = n0;
        this->n1 = n1;
        this->n2 = n2;
        this->data = new T [n0*n1*n2];
    }

    void destroy()
    {
        delete[] this->data;
    }

    T &operator()( int i0, int i1, int i2 )
    {
        return data[ i0 * n1 * n2 + i1 * n2 + i2 ];
    }

    const T &operator()( int i0, int i1, int i2 ) const
    {
        return data[ i0 * n1 * n2 + i1 * n2 + i2 ];
    }

    _Array3_get_index_1< Array3, T > operator[]( int i0 )
    {
        return _Array3_get_index_1< Array3< T >, T >( this, i0 );
    }

    _Array3_get_index_1< const Array3, const T > operator[]( int i0 ) const
    {
        return _Array3_get_index_1< const Array3< T >, const T >( this, i0 );
    }

    void write( std::ofstream & output )
    {
        output.write( reinterpret_cast< char * >( this->data ), sizeof( T ) * this->n0 * this->n1 * this->n2 );
    }

    void read( std::ifstream & input )
    {
        input.read( reinterpret_cast< char * >( this->data ), sizeof( T ) * this->n0 * this->n1 * this->n2 );
    }
    T * get_pointer()
    {
        return this->data;
    }
    void put_pointer(T *reference)
    {
        this->data=reference;
    }
private:
    int n0, n1, n2;
    T * data;
};


template< typename A, typename T >
class _Array2_get_index_1
{
public:
    _Array2_get_index_1( A * array, int i0 )
        : array( array ), i0( i0 )
    {}

    T &operator[]( int i1 )
    {
        return (*this->array)( this->i0, i1 );
    }

private:
    A * array;
    const int i0;
};

template< typename T >
class Array2
{
public:
    Array2()
    {}

    Array2( int n0, int n1, T * data )
        : n0( n0 ), n1( n1 ), data( data )
    {}

    void create( int n0, int n1 )
    {
        this->n0 = n0;
        this->n1 = n1;
        this->data = new T [n0*n1];
    }

    void destroy()
    {
        delete[] this->data;
    }

    T &operator()( int i0, int i1 )
    {
        return data[ i0 * n1 + i1 ];
    }

    const T &operator()( int i0, int i1 ) const
    {
        return data[ i0 * n1 + i1 ];
    }

    _Array2_get_index_1< Array2, T > operator[]( int i0 )
    {
        return _Array2_get_index_1< Array2< T >, T >( this, i0 );
    }

    _Array2_get_index_1< const Array2, const T > operator[]( int i0 ) const
    {
        return _Array2_get_index_1< const Array2< T >, const T >( this, i0 );
    }

    void write( std::ofstream & output )
    {
        output.write( reinterpret_cast< char * >( this->data ), sizeof( T ) * this->n0 * this->n1 );
    }

    void read( std::ifstream & input )
    {
        input.read( reinterpret_cast< char * >( this->data ), sizeof( T ) * this->n0 * this->n1 );
    }
    T * get_pointer()
    {
        return this->data;
    }
    void put_pointer(T *reference)
    {
        this->data=reference;
    }
private:
    int n0, n1;
    T * data;
};


template< typename T >
class Array1
{
public:
    Array1()
    {}

    Array1( int n0, T * data )
        : n0( n0 ), data( data )
    {}

    void create( int n0 )
    {
        this->n0 = n0;
        this->data = new T [n0];
    }

    void destroy()
    {
        delete[] this->data;
    }

    T &operator()( int i0 )
    {
        return data[ i0 ];
    }

    const T &operator()( int i0 ) const
    {
        return data[ i0 ];
    }

    T & operator[]( int i0 )
    {
        return data[ i0 ];
    }

    const T & operator[]( int i0 ) const
    {
        return data[ i0 ];
    }

    void write( std::ofstream & output )
    {
        output.write( reinterpret_cast< char * >( this->data ), sizeof( T ) * this->n0 );
    }

    void read( std::ifstream & input )
    {
        input.read( reinterpret_cast< char * >( this->data ), sizeof( T ) * this->n0 );
    }
    T * get_pointer()
    {
        return this->data;
    }
    void put_pointer(T *reference)
    {
        this->data=reference;
    }
private:
    int n0;
    T * data;
};


template< typename T >
class XYZ
{
public:
    XYZ( const T & x, const T & y, const T & z ) : x( x ), y( y ), z( z ) {}

    double norm2_squared() const { return this->x*this->x + this->y*this->y + this->z*this->z; }
    double norm2() const { return sqrt( this->norm2_squared() ); }

    template< typename S >
    XYZ< S > cast() const { return XYZ< S >( this->x, this->y, this->z ); }

    T x, y, z;
};

template< typename T >
XYZ< T > operator-( const XYZ< T > & r )
{
    return XYZ< T >( -r.x, -r.y, -r.z );
}

template< typename T >
XYZ< T > operator+( const XYZ< T > & l, const XYZ< T > & r )
{
    return XYZ< T >( l.x + r.x, l.y + r.y, l.z + r.z );
}

template< typename T >
XYZ< T > operator+( const T & l, const XYZ< T > & r )
{
    return XYZ< T >( l + r.x, l + r.y, l + r.z );
}

template< typename T >
XYZ< T > operator+( const XYZ< T > & l, const T & r )
{
    return XYZ< T >( l.x + r, l.y + r, l.z + r );
}

template< typename T >
XYZ< T > operator-( const XYZ< T > & l, const XYZ< T > & r )
{
    return XYZ< T >( l.x - r.x, l.y - r.y, l.z - r.z );
}

template< typename T >
XYZ< T > operator-( const T & l, const XYZ< T > & r )
{
    return XYZ< T >( l - r.x, l - r.y, l - r.z );
}

template< typename T >
XYZ< T > operator-( const XYZ< T > & l, const T & r )
{
    return XYZ< T >( l.x - r, l.y - r, l.z - r );
}

template< typename T >
XYZ< T > operator*( const XYZ< T > & l, const XYZ< T > & r )
{
    return XYZ< T >( l.x * r.x, l.y * r.y, l.z * r.z );
}

template< typename T >
XYZ< T > operator*( const T & l, const XYZ< T > & r )
{
    return XYZ< T >( l * r.x, l * r.y, l * r.z );
}

template< typename T >
XYZ< T > operator*( const XYZ< T > & l, const T & r )
{
    return XYZ< T >( l.x * r, l.y * r, l.z * r );
}

template< typename T >
XYZ< T > operator/( const XYZ< T > & l, const XYZ< T > & r )
{
    return XYZ< T >( l.x / r.x, l.y / r.y, l.z / r.z );
}

template< typename T >
XYZ< T > operator/( const T & l, const XYZ< T > & r )
{
    return XYZ< T >( l / r.x, l / r.y, l / r.z );
}

template< typename T >
XYZ< T > operator/( const XYZ< T > & l, const T & r )
{
    return XYZ< T >( l.x / r, l.y / r, l.z / r );
}
