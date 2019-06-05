#ifndef INC_COMPLEXPRECISION
#define INC_COMPLEXPRECISION

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002-2016
 *                       Future Team Aps 
 *                       Denmark
 *
 *                       All Rights Reserved
 *
 *   This source file is subject to the terms and conditions of the
 *   Future Team Software License Agreement which restricts the manner
 *   in which it may be used.
 *   Mail: hve@hvks.com
 *
 *******************************************************************************
*/

/*
 *******************************************************************************
 *
 *
 * Module name     :   complexprecision.h
 * Module ID Nbr   :   
 * Description     :   Arbitrary complex precision class
 *                     Actually it a general complex class that works with both
 *                     standard types like int, float, double and float_precision
 *                     or int_precision or for that matter with any other types
 * --------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version	Author/Date		Description of changes
 * -------  -----------		----------------------
 * 01.01	HVE/030331		Initial release
 * 01.02    hve/060203     Minor declaration bug fixed
 * 01.03    HVE/060217     Error in the formula for exp and log corrected
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/


/* define version string */
static char _VC_[] = "@(#)complexprecision.h 01.03 -- Copyright (C) Future Team Aps";


// Complex Precision template class
template<class _Ty> class complex_precision {
   _Ty re, im;
   public:
      typedef _Ty value_type;

      // constructor
      complex_precision( const _Ty& r = (_Ty)0, const _Ty& i = (_Ty)0 ) : re(r), im(i) {}
      // constructor for any other type to _Ty
      template<class _X> complex_precision( const complex_precision<_X>& a ) : re((_Ty)a.real()), im((_Ty)a.imag()) {}
      
      // Coordinate functions
      _Ty real() const { return re; }
      _Ty imag() const { return im; }
      _Ty real( const _Ty& r )   { return ( re = r ); }
	  _Ty real( const complex_precision<_Ty>& r )   { return ( re = r.real() ); }
      _Ty imag( const _Ty& i )   { return ( im = i ); }
	  _Ty imag( const complex_precision<_Ty>& i )   { return ( im = i.imag() ); }
      _Ty norm() const { return re * re + im * im; }
      _Ty abs () const { _Ty a(re), b(im), c(re);  
                       a = (_Ty)re < (_Ty)0 ? -re : re;
                       b = (_Ty)im < (_Ty)0 ? -im : im;
                       if( a >= b )
                          {
                          c = im / re;
                          return a * sqrt( (_Ty)1 + c * c );
                          }
                       else
                          {
                          c = re / im;
                          return b * sqrt( (_Ty)1 + c * c );
                          }
                       }
      _Ty arg() const  { return atan2( im, re ); }
      complex_precision<_Ty> conj() const { complex_precision<_Ty> x(*this); x.real( re ); x.imag( -im ); return x; }
      _Ty *ref_real()   { return &re; }
      _Ty *ref_imag()   { return &im; }

      // Essential operators
      complex_precision<_Ty>& operator= ( const complex_precision<_Ty>& x )   { re = x.real(); im = x.imag(); return *this; }
      complex_precision<_Ty>& operator+=( const complex_precision<_Ty>& x )   { re += x.real(); im += x.imag(); return *this; }
      complex_precision<_Ty>& operator-=( const complex_precision<_Ty>& x )   { re -= x.real(); im -= x.imag(); return *this; }
      complex_precision<_Ty>& operator*=( const complex_precision<_Ty>& x )   { _Ty w(x.real()); w = re * x.real() - im * x.imag(); im = re * x.imag() + im * x.real(); re = w; return *this; }
      complex_precision<_Ty>& operator/=( const complex_precision<_Ty>& x );  // Too big to have here
     
	  class divide_by_zero {};
   };

template<class _Ty> std::ostream& operator<<( std::ostream& strm, const complex_precision<_Ty>& a )
	{ return strm << "(" << a.real() << "," << a.imag() << ")"; }

template<class _Ty> std::istream& operator>>( std::istream& strm, complex_precision<_Ty>& c ) 
   {
   _Ty re, im; char ch;
   if( strm >> ch && ch != '(')
		strm.putback(ch), strm >> re, im = (_Ty)0;
	else
      if( strm >> re >> ch && ch != ',')
		   if( ch == ')')
			   im = (_Ty)0;
		   else 
            strm.putback( ch ); // strm.setstate(std::ios::failbit);
	   else
         if( strm >> im >> ch && ch != ')')
			   strm.putback( ch ); //, strm.setstate(ios_base::failbit);
	if(!strm.fail())
		c = complex_precision<_Ty>( re, im );

   return strm;
   }


// Arithmetic
template<class _Ty> complex_precision<_Ty> operator+( const complex_precision<_Ty>&, const complex_precision<_Ty>& );  // Binary
template<class _Ty> complex_precision<_Ty> operator+( const complex_precision<_Ty>& );                                 // Unary
template<class _Ty> complex_precision<_Ty> operator-( const complex_precision<_Ty>&, const complex_precision<_Ty>& );  // Binary
template<class _Ty> complex_precision<_Ty> operator-( const complex_precision<_Ty>& );                                 // Unary
template<class _Ty> complex_precision<_Ty> operator*( const complex_precision<_Ty>&, const complex_precision<_Ty>& );  // Binary
template<class _Ty> complex_precision<_Ty> operator/( const complex_precision<_Ty>&, const complex_precision<_Ty>& );  // Binary
                                                                                                                       
// Boolean Comparision Operators
template<class _Ty> bool operator==( const complex_precision<_Ty>&, const complex_precision<_Ty>& );
template<class _Ty> bool operator!=( const complex_precision<_Ty>&, const complex_precision<_Ty>& );

template<class _Ty> complex_precision<_Ty> polar( _Ty magn, _Ty theta = (_Ty)0 )
   {
   if( theta == (_Ty)0 )
      return complex_precision<_Ty>( magn, 0 );

   }

// Essential operators
//
template<class _Ty> complex_precision<_Ty>& complex_precision<_Ty>::operator/=( const complex_precision<_Ty>& y )
   {
   if( y.real() == (_Ty)0 && y.imag() == (_Ty)0 )
      { throw complex_precision<_Ty>::divide_by_zero(); return *this; }

   if( ( y.real() < (_Ty)0 ? -y.real() : y.real() ) >= ( y.imag() < (_Ty)0 ? -y.imag() : y.imag() ) )
      {
      _Ty t(y.imag() / y.real() );        // Force same precision as y 
      _Ty t2(y.real() + y.imag() * t );   // Force same precision as y 
      _Ty t3(re + im * t );               // Force same precision as y 
      im -= re * t;
      im /= t2;
      re = t3 / t2;
      }
   else
      {
      _Ty t(y.real() / y.imag() );        // Force same precision as y 
      _Ty t2(y.real() * t + y.imag() );   // Force same precision as y 
      _Ty t3(re * t + im );               // Force same precision as y 
      im *= t;
      im -= re;
      im /= t2;
      re = t3 / t2;
      }

   return *this;
   }

template<class _Ty> complex_precision<_Ty> operator+( const complex_precision<_Ty>& x, const complex_precision<_Ty>& y )
   {
   complex_precision<_Ty> c(x);

   c += y;
   return c;
   }


template<class _Ty> complex_precision<_Ty> operator+( const complex_precision<_Ty>& x )
   {
   return x;
   }

template<class _Ty> complex_precision<_Ty> operator-( const complex_precision<_Ty>& x, const complex_precision<_Ty>& y )
   {
   complex_precision<_Ty> c(x);

   c -= y;
   return c;
   }

template<class _Ty> complex_precision<_Ty> operator-( const complex_precision<_Ty>& x )
   {
   complex_precision<_Ty> c(x);

   c.real( -c.real() );
   c.imag( -c.imag() );
   return c;
   }

template<class _Ty> complex_precision<_Ty> operator*( const complex_precision<_Ty>& x, const complex_precision<_Ty>& y )
   {
   complex_precision<_Ty> c(x);

   c *= y;
   return c;
   }

template<class _Ty> complex_precision<_Ty> operator/( const complex_precision<_Ty>& x, const complex_precision<_Ty>& y )
   {
   complex_precision<_Ty> c(x);

   c /= y;
   return c;
   }

template<class _Ty> bool operator==( const complex_precision<_Ty>& a, const complex_precision<_Ty>& b )
   {
   return a.real() == b.real() && a.imag() == b.imag();
   }

template<class _Ty> bool operator!=( const complex_precision<_Ty>& a, const complex_precision<_Ty>& b )
   {
   return a.real() != b.real() || a.imag() != b.imag();
   }


template<class _Ty> _Ty abs(const complex_precision<_Ty>& x )
   {
   return x.abs();
   }

template<class _Ty> complex_precision<_Ty> sqrt( const complex_precision<_Ty> x )
   {
   _Ty w(x.real()), c(x.real()), d(x.real());  // Force the local variables to the same precisions as x (Only for float_precision)

   if( x.real() == (_Ty)0 && x.imag() == (_Ty)0 )
      w = (_Ty)0;
   else
      {
      c = x.real() < (_Ty)0 ? -x.real() : x.real();
      d = x.imag() < (_Ty)0 ? -x.imag() : x.imag();
      if( c < d )
         {
         _Ty t(c / d);
         if( t < (_Ty)0 ) t = -t;
         w = sqrt( d ) * sqrt( ( t + sqrt( (_Ty)1 + t * t ) ) / (_Ty)2 );
         }
      else
         {
         _Ty t(d / c);
         w = sqrt( c ) * sqrt( ( (_Ty)1 + sqrt( (_Ty)1 + t * t ) ) / (_Ty)2 );
         }
      }

   if(  w == (_Ty)0 )
      return complex_precision<_Ty>( (_Ty)w );
   
   if( x.real() >= (_Ty)0 )
      return complex_precision<_Ty>( w, x.imag() / ( (_Ty)2 * w ) );

   if( x.imag() >= (_Ty)0 )
      return complex_precision<_Ty>( d / ( (_Ty)2 * w ), w );

   return complex_precision<_Ty>( d / ( (_Ty)2 * w ), -w );
   }

template<class _Ty> complex_precision<_Ty> log( const complex_precision<_Ty> x )
   {
   return complex_precision<_Ty>( log( x.abs() ), atan( x.imag() / x.real() ) );
   }

template<class _Ty> complex_precision<_Ty> log10( const complex_precision<_Ty> x )
   {
   return log( x ) / complex_precision<_Ty>( log( (_Ty)10 ), 0 );
   }

template<class _Ty> complex_precision<_Ty> exp( const complex_precision<_Ty> x )
   {
   return complex_precision<_Ty>( exp( x.real() ) * cos( x.imag() ), exp( x.real() ) * sin( x.imag() ) );
   }

template<class _Ty> complex_precision<_Ty> pow( const complex_precision<_Ty> x, const complex_precision<_Ty> y )
   {
   complex_precision<_Ty> z(x);     // Force same precision as x

   z = log( x );
   z *= y;
   z = exp( z );

   return z;
   }

#endif
