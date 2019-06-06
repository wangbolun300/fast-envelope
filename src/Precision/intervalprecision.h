#ifndef INC_INTERVAL
#define INC_INTERVAL

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002-2014
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
 * Module name     :   intervalprecision.h
 * Module ID Nbr   :   
 * Description     :   Interval Precision arithmetic template class
 *                     Works with the int_precision and float_precision classes
 * --------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version	Author/Date		Description of changes
 * -------  -----------		----------------------
 * 01.01	HVE/020209		Initial release
 * 01.02    HVE/030421		Optimized the * operator to reduce the 8 multiplications to 2.
 * 01.03	HVE/JUN-26-2014	Added is_empty(), contains_zero() method to the class
 * 01.04	HVE/JUN-27-2014	Corrected several errors in in cin >> template function
 * 01.05	HVE/JUN-28-2014	Added is_class() method for getting the interval classification
 *							and width() method for the interval width
 * 01.06	HVE/JUN-30-2014	An error was corrected for interval subtraction of float_preicsion numbers
 *							Also added the method bool contain() for test if a float or interval is included in the interval
 * 01.07	HVE/JUL-6-2014	Corrected an error in /= for the software emulation of of float & double
 * 01.08	HVE/JUL-13-2014	Added Hardware support for interval arithmetic when applicable. Also fix several errors in the 
 *							implementation of sqrt, log, log10, exp and pow functions. Also added new method is_class(), is_empty()
 * 01.09	HVE/JUL-15-2014	Added support for Sin(), Cos() and Tan() interval trigonometric functions.
 * 01.10	HVE/JUL-17-2014	Added support for atan() interval trigonometric function
 * 01.11	HVE/JUL-22-2014 Found a bug that floating point was not reset to near (default by IEEE754) after a hardware supported multiplication
 * 01.12	HVE/JUL-22-2014	Added support for asin() interval trigonometric function
 * 01.13	HVE/JUL-29-2014	Added support for interval versions of LN2, LN10 and PI
 * 01.14	HVE/AUG-10-2014	Added support for mixed mode arithmetic for interval +,- classes
 * 01.15	HVE/JUN-20-2015	Fixed and undeclare variable x when comopiling with no interval hardware support
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/


#include <float.h>
#include <algorithm>

namespace precision {
/* define version string */
static char _VinterP_[] = "@(#)intervalprecision.h 01.15 -- Copyright (C) Future Team Aps";


// HARDWARE_SUPPORT controlled if IEEE754 floating point control can be used for interval arithmetic.
// if not used interval arithmetic will be handle in software 
//#define HARDWARE_SUPPORT

/// The four different interval classification
/// # ZERO			a=0 && b=0
/// # POSITIVE		a>=0 && b>0
/// # NEGATIVE		a<0 && b<=0
/// # MIXED			a<0 && b>0
enum int_class { NO_CLASS, ZERO, POSITIVE, NEGATIVE, MIXED };

// 
// Interval class
// Realisticly the class Type can be float, double, int_precision, & float_precision
// Since float and double are done unsing the Intel cpu (H/W) and using "specilization" 
// the int_precision and float_precision is done using the arbitrary precision packages
// Since their is no way to specific portable ROUND_DOWN and ROUND_UP mode the template class 
// heavily use specilization. For any other type besides float, double, int_precision and float_precision
// the operations is not defined
//
template<class _IT> class interval {
   _IT low, high;
   public:
      typedef _IT value_type;

      // constructor. zero, one or two arguments for type _IT
      interval()							{ low = _IT(0); high = _IT(0); }
	  interval( const _IT& d )				{ low = _IT(d); high = _IT(d); }
	  interval( const _IT& l, const _IT& h) { if( l < h ) { low =l; high = h; } else { low = h; high = l; } }
	  // Constrcutor for mixed type _IT != _X (base types). Allows auto construction of e.g. interval<float_precision> x(float)
	 template <class _X> interval(const _X& x) { low = _IT(x); high = _IT(x); }

      // constructor for any other type to _IT. Both up and down conversion possible
      template<class X> interval( const interval<X>& a ) /*: low(_IT(a.lower())), high( _IT(a.upper())) */ 
	  {
		/*  int is = sizeof(_IT), js = sizeof(X);
		  cout << "_IT=" << typeid(_IT).name() << " X=" << typeid(X).name() << endl;
		  fpdown(); low = _IT(a.lower());
		  fpup(); high = _IT(a.upper());
		  fpnear(); 
		  if (typeid(_IT) == typeid(double) && typeid(X)==typeid(float))  // upscaling float to double
			  {
			  double u = high;
			  u += u * 0.5f * FLT_EPSILON;
			 // high = _IT(u);
			  }
			  */
		  if (a.lower() < a.upper()) { fpdown();  low = _IT(a.lower()); fpup();  high = _IT(a.upper()); fpnear();  }
		  else { fpdown();  low = _IT(a.upper()); fpup();  high = _IT(a.lower()); fpnear(); }
	  }
	 
      // Coordinate functions
      _IT upper() const					{ return high; }
      _IT lower() const					{ return low; }
      _IT upper( const _IT& u )			{ return ( high = u ); }
      _IT lower( const _IT& l )			{ return ( low = l ); }

      _IT center() const				{ return ( high + low ) / _IT(2); }
      _IT radius() const				{ _IT r; r =( high - low ) / _IT(2); if( r < _IT(0) ) r = -r; return r; }
	  _IT width() const					{ _IT r; r = high - low; if (r < _IT(0)) r = -r; return r; }
	
	  bool contain_zero() const			{ return low <= _IT(0) && _IT(0) <= high; }  // Obsolete. use contains() instead.
	  bool contain( const _IT& f=_IT(0)){ return low <= f && f <= high;  }
	  bool contain(const interval<_IT>& i) { return low <= i.lower() && i.upper() <= high; }
	  bool is_empty() const				{ return high < low; }
	
	  enum int_class is_class() const	{ 
										if (low == _IT(0) && high == _IT(0)) return ZERO;
										if (low >= _IT(0) && high > _IT(0)) return POSITIVE;
										if (low < _IT(0) && high <= _IT(0)) return NEGATIVE;
										if (low < _IT(0) && high > _IT(0)) return MIXED;
										return NO_CLASS;
										}
	  // Conversion methods. Safer and less ambiguios than overloading implicit/explivit conversion operators
//	  std::string toString() const					{ return _float_precision_ftoa(this); }
	//  int_precision to_int_precision() const		{ std::string s = _float_precision_ftoainteger(this); return (int_precision)((char *)s.c_str()); }

	  // Operators
	  operator short() const			{ return (short)((high + low) / _IT(2)); }				// Conversion to short 	 
	  operator int() const				{ return (int)( ( high + low ) / _IT(2) ); }			// Conversion to int 
	  operator long() const				{ return (long)((high + low) / _IT(2)); }				// Conversion to long 
	  operator unsigned short() const	{ return (unsigned short)((high + low) / _IT(2)); }		// Conversion to unsigned short 	 
	  operator unsigned int() const		{ return (unsigned int)((high + low) / _IT(2)); }		// Conversion to unsigned int 
	  operator unsigned long() const	{ return (unsigned long)((high + low) / _IT(2)); }		// Conversion to unsigned long 
	  operator double() const			{ return (double)( ( high + low ) / _IT(2) ); }			// Conversion to double 
	  operator float() const			{ return high == low? (float)low : (float)((high + low) / _IT(2)); }				// Conversion to float 
	  operator int_precision() const	{ return (int_precision)((high + low) / _IT(2)); }		// Conversion to int_precision 
	  operator float_precision() const	{ return (float_precision)((high + low) / _IT(2)); }	// Conversion to float_precision 

      _IT *ref_lower()					{ return &low; }
      _IT *ref_upper()					{ return &high; }

      // Essential operators
      interval<_IT>& operator= ( const interval<_IT>& );
      interval<_IT>& operator+=( const interval<_IT>& );
      interval<_IT>& operator-=( const interval<_IT>& );
      interval<_IT>& operator*=( const interval<_IT>& );
      interval<_IT>& operator/=( const interval<_IT>& );
	  interval<_IT>& operator&=( const interval<_IT>& );
	  interval<_IT>& operator|=( const interval<_IT>& );
	  interval<_IT>& operator^=( const interval<_IT>& );

	  // Exception class
	  class bad_int_syntax {};
	  class bad_float_syntax {};
	  class out_of_range   {};
	  class divide_by_zero {};
	  class domain_error   {};
	  class base_error		{};
   };


// Unary and Binary arithmetic
// Arithmetic + Binary and Unary
template <class _IT, class _X> inline interval<_IT> operator+( const interval<_IT>&, const _X&);
template <class _IT, class _X> inline interval<_IT> operator+( const _X&, const interval<_IT>&);
inline interval<float_precision> operator+(const interval<float_precision>&, const float_precision&);	// Specialization for interval<float_precision> and float_precision
inline interval<float_precision> operator+(const float_precision&, const interval<float_precision>&);	// Specialization for ifloat_precision and interval<float_precision>

template<class _IT> interval<_IT> operator+( const interval<_IT>&, const interval<_IT>& ); 
template<class _IT> interval<_IT> operator+( const interval<_IT>& );									// Unary 

// Arithmetic - Binary and Unary
template <class _IT, class _X> inline interval<_IT> operator-(const interval<_IT>&, const _X&);
template <class _IT, class _X> inline interval<_IT> operator-(const _X&, const interval<_IT>&);
inline interval<float_precision> operator-(const interval<float_precision>&, const float_precision&);	// Specialization for interval<float_precision> and float_precision
inline interval<float_precision> operator-(const float_precision&, const interval<float_precision>&);	// Specialization for ifloat_precision and interval<float_precision>

template<class _IT> interval<_IT> operator-( const interval<_IT>&, const interval<_IT>& );
template<class _IT> interval<_IT> operator-( const interval<_IT>& );									// Unary

// Arithmetic * Binary
template <class _IT, class _X> inline interval<_IT> operator*(const interval<_IT>&, const _X&);
template <class _IT, class _X> inline interval<_IT> operator*(const _X&, const interval<_IT>&);
inline interval<float_precision> operator*(const interval<float_precision>&, const float_precision&);	// Specialization for interval<float_precision> and float_precision
inline interval<float_precision> operator*(const float_precision&, const interval<float_precision>&);	// Specialization for ifloat_precision and interval<float_precision>

template<class _IT> interval<_IT> operator*( const interval<_IT>&, const interval<_IT>& );

// Arithmetic / Binary
template <class _IT, class _X> inline interval<_IT> operator/(const interval<_IT>&, const _X&);
template <class _IT, class _X> inline interval<_IT> operator/(const _X&, const interval<_IT>&);
inline interval<float_precision> operator/(const interval<float_precision>&, const float_precision&);	// Specialization for interval<float_precision> and float_precision
inline interval<float_precision> operator/(const float_precision&, const interval<float_precision>&);	// Specialization for ifloat_precision and interval<float_precision>

template<class _IT> interval<_IT> operator/( const interval<_IT>&, const interval<_IT>& );

template<class _IT> interval<_IT> operator&( const interval<_IT>&, const interval<_IT>& );
template<class _IT> interval<_IT> operator|( const interval<_IT>&, const interval<_IT>& );
template<class _IT> interval<_IT> operator^( const interval<_IT>&, const interval<_IT>& );

// Boolean Comparison Operators
template<class _IT> bool operator==(const interval<_IT>&, const interval<_IT>&);
template<class _IT> bool operator!=(const interval<_IT>&, const interval<_IT>&);

// Other functions
template<class _IT> interval<_IT> abs(const interval<_IT>&);

// Manifest Constants like PI, LN2 and LN10
inline interval<float> int_pifloat();
inline interval<double> int_pidouble();
inline interval<float_precision> int_pi(const unsigned int);
inline interval<float> int_ln2float();
inline interval<double> int_ln2double();
inline interval<float_precision> int_ln2(const unsigned int);
inline interval<float> int_ln10float();
inline interval<double> int_ln10double();
inline interval<float_precision> int_ln10(const unsigned int);

// Elementary functions
inline interval<float> sqrt(const interval<float>&);
inline interval<double> sqrt(const interval<double>&);
inline interval<float_precision> sqrt(const interval<float_precision>&);

inline interval<float> log( const interval<float>& );
inline interval<double> log(const interval<double>&);
inline interval<float_precision> log(const interval<float_precision>&);

inline interval<float> log10(const interval<float>&);
inline interval<double> log10(const interval<double>&);
inline interval<float_precision> log10(const interval<float_precision>&);

inline interval<float> exp(const interval<float>&, const float);
inline interval<double> exp(const interval<double>&, const double);
inline interval<float_precision> exp(const interval<float_precision>&);

inline interval<float> pow(const interval<float>&, const float );
inline interval<double> pow( const interval<double>&, const double );
inline interval<float_precision> pow( const interval<float_precision>& );

// Trigonometric functions
inline interval<float> sin(const interval<float>&);
inline interval<double> sin(const interval<double>&);
inline interval<float_precision> sin(const interval<float_precision>&);

inline interval<float> cos(const interval<float>&);
inline interval<double> cos(const interval<double>&);
inline interval<float_precision> cos(const interval<float_precision>&);

inline interval<float> tan(const interval<float>&);
inline interval<double> tan(const interval<double>&);
inline interval<float_precision> tan(const interval<float_precision>&);

inline interval<float> asin(const interval<float>&);
inline interval<double> asin(const interval<double>&);
inline interval<float_precision> asin(const interval<float_precision>&);

inline interval<float> acos(const interval<float>&);
inline interval<double> acos(const interval<double>&);
inline interval<float_precision> acos(const interval<float_precision>&);

inline interval<float> atan(const interval<float>&);
inline interval<double> atan(const interval<double>&);
inline interval<float_precision> atan(const interval<float_precision>&);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//	Floating point control for the IEEE754 hardware. Only fo non managed application
//	Enable by defined #define HARDWARE_SUPPORT
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef HARDWARE_SUPPORT
inline void fpnear()
	{
	unsigned int f87_cw, sse2_cw;
	int cc;
	//_controlfp_s( &currentControl, _RC_NEAR, _MCW_RC);
	cc=__control87_2(_RC_NEAR, _MCW_RC, &f87_cw, &sse2_cw );
	}

inline void fpdown()
	{
	unsigned int f87_cw, sse2_cw;
	int cc;
	//_controlfp_s( &currentControl, _RC_DOWN, _MCW_RC);
	cc=__control87_2(_RC_DOWN, _MCW_RC, &f87_cw, &sse2_cw);
	}

inline void fpup()
	{
	unsigned int f87_cw, sse2_cw;
	int cc;
	//_controlfp_s( &currentControl, _RC_UP, _MCW_RC);
	cc=__control87_2(_RC_UP, _MCW_RC, &f87_cw, &sse2_cw);
	}
#else
inline void fpnear()	{}
inline void fpdown()	{}
inline void fpup()		{}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//   End Floating point control for the IEEE754 hardware
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Output Operator <<
//
template<class _Ty> inline std::ostream& operator<<( std::ostream& strm, interval<_Ty>& a ) { return strm << "[" << a.lower() << "," << a.upper() << "]"; }

// Input operator >>
//
template<class _Ty> inline std::istream& operator>>( std::istream& strm, interval<_Ty>& c ) 
   {
   _Ty l, u; char ch;
   if( strm >> ch && ch != '[')
      strm.putback(ch), strm >> l, u = l;
	else
      if( strm >> l >> ch && ch != ',')
	      if( ch == ']')
	         u = l;
	      else 
            strm.putback( ch ); // strm.setstate(std::ios::failbit);
	   else
         if( strm >> u >> ch && ch != ']')
	         strm.putback( ch ); //, strm.setstate(ios_base::failbit);
	
   if(!strm.fail())
	   c = interval<_Ty>( l, u );

   return strm;
   }
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Essential Operators =,+=,-=,*=,/=
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////


// Assignment operator. Works for all class types
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator=( const interval<_IT>& a )
   {
   low = a.lower();
   high = a.upper();
   return *this;
   }

// += operator. Works all other classes. 
// Please note that this is for all integer classes. interval<int>, interval<long>, interval<int_precision>
// were there os no loss of precision
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator+=( const interval<_IT>& a )
   {
   low += a.lower();
   high += a.upper();
   return *this;
   }

// Specilization for float_precision and +=
//
inline interval<float_precision>& interval<float_precision>::operator+=( const interval<float_precision>& a )
   {
   low.mode( ROUND_DOWN );
   low += a.lower();
   high.mode( ROUND_UP );
   high += a.upper();
   
   return *this;
   }


// Specilization for float and +=
// That can work with both managed an unmanged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
inline interval<float>& interval<float>::operator+=( const interval<float>& a )
	{
#ifdef HARDWARE_SUPPORT
	fpdown();
	low += a.lower();
	fpup();
	high += a.upper();
	fpnear();
#else
	interval<float_precision> rhs(*this), lhs(a);

	rhs += lhs;
	low=rhs.lower();
	high=rhs.upper();
#endif
	return *this;
	}

// Specilization for double and +=
// That can work with both managed an unmanged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
inline interval<double>& interval<double>::operator+=( const interval<double>& a )
	{
#ifdef HARDWARE_SUPPORT
	fpdown();
	low += a.lower();
	fpup();
	high += a.upper();
	fpnear();
#else
	interval<float_precision> rhs(*this), lhs(a);

	rhs += lhs;
	low=rhs.lower();
	high=rhs.upper();
#endif
	return *this;
	}

// -= operator. Works all other classes. 
// Please note that this is for all integer classes. interval<int>, interval<long>, interval<int_precision>
// were there is no loss of precision
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator-=( const interval<_IT>& a )
   {
   low -= a.high;
   high -= a.low;
   return *this;
   }

// Specilization for float_precision and -=
//
inline interval<float_precision>& interval<float_precision>::operator-=( const interval<float_precision>& a )
   {
   low.mode( ROUND_DOWN );
   low -= a.upper();
   high.mode( ROUND_UP );
   high -= a.lower();
   
   return *this;
   }

// Specilization for float and -=
// That can work with both managed an unmanged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
inline interval<float>& interval<float>::operator-=( const interval<float>& a )
	{
#ifdef HARDWARE_SUPPORT
	fpdown();
	low -= a.high;
	fpup();
	high -= a.low;
	fpnear();
#else
	interval<float_precision> rhs(*this), lhs(a);

	rhs -= lhs;
	low=rhs.lower();
	high=rhs.upper();
#endif
	return *this;
	}

// Specilization for double and -=
// That can work with both managed an unmanged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
inline interval<double>& interval<double>::operator-=( const interval<double>& a )
	{
#ifdef HARDWARE_SUPPORT
	fpdown();
	low -= a.high;
	fpup();
	high -= a.low;
	fpnear();
#else
	interval<float_precision> rhs(*this), lhs(a);

	rhs -= lhs;
	low=rhs.lower();
	high=rhs.upper();
#endif
	return *this;
	}

// Works all other classes. 
// Please note that this is for all interger classes. interval<int>, interval<long>, interval<int_precision>
// were there is no loss of precision
// Instead of doing the mindless low = MIN(low*a.high, low*a.low,high*a.low,high*a.high) and
// high = MAX(low*a.high, low*a.low,high*a.low,high*a.high) requiring a total of 8 multiplication
//
//   low, high, a.low, a.high    result
//    +     +     +     +        +  +  [ low*a.low, high*a.high ]
//    +     +     -     +        -  +  [ high*a.low, high*a.high ]
//    +     +     -     -        -  -  [ high*a.low, low*a.high ]
//    -     +     +     +        -  +  [ low*a.high, high*a.high ]  
//    -     +     -     +        -  +  [ MIN(low*a.high,high*a.low), MAX(low*a.low,high*a.high) ]
//    -     +     -     -        -  -  [ high*a.low, low*a.low ]
//    -     -     +     +        -  -  [ low*a.high, high,a.low ]
//    -     -     -     +        -  -  [ low*a.high, low*a.low ]
//    -     -     -     -        +  +  [ high*a.high, low * a.low ]
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator*=( const interval<_IT>& a )
   {
   _IT l, h, t;

   if( low >= _IT(0) ) // 
      { // both low and high >= 0
      if( a.lower() >= _IT(0) )
         { // a.low >=0, a.high >= 0
         l = low * a.lower();
         h = high * a.upper();
         }
      else
         if( a.upper() >= _IT(0) )
            {  //  a.low < 0, a.high >= 0
            l = high * a.lower();
            h = high * a.upper();
            }
         else
            { // a.low and a.high < 0 
            l = high * a.lower();
            h = low * a.upper();
            }
      }
   else
      if( high >= _IT(0) )
         {  // low < 0, high >= 0
         if( a.lower() >= _IT(0) )
            { // a.low >=0, a.high >= 0
            l = low * a.upper();
            h = high * a.upper();
            }
         else
            if( a.upper() >= _IT(0) )
               {  //  a.low < 0, a.high >= 0
               l = low * a.upper(); if ( l > ( t = high * a.lower() ) ) l = t;
               h = high * a.upper(); if ( h < ( t = low * a.lower() ) ) h = t;
               }
            else
               { // a.low and a.high < 0 
               l = high * a.lower();
               h = low * a.lower();
               }
         }
      else
         { // low and high are < 0 
         if( a.lower() >= _IT(0) )
            { // a.low >=0, a.high >= 0
            l = low * a.upper();
            h = high * a.lower();
            }
         else
            if( a.upper() >= _IT(0) )
               {  //  a.low < 0, a.high >= 0
               l = low * a.upper(); 
               h = low * a.lower();
               }
            else
               { // a.low and a.high < 0 
               l = high * a.upper();
               h = low * a.lower();
               }
        }

   low = l;
   high = h;

   return *this;
   }


// Specilization for float_precision and *= operator
//
inline interval<float_precision>& interval<float_precision>::operator*=( const interval<float_precision>& a )
   {
   float_precision l, h, t;

   l.precision( low.precision() );
   h.precision( low.precision() );   
   t.precision( low.precision() );

   l.mode( ROUND_DOWN );
   h.mode( ROUND_UP );

   if( low.sign() > 0 ) // 
      { // both low and high >= 0
      if( a.lower().sign() > 0 )
         { // a.low >=0, a.high >= 0
         l = low;  l *= a.lower();
         h = high; h *= a.upper();
         }
      else
         if( a.upper().sign() > 0 )
            {  //  a.low < 0, a.high >= 0
            l = high;  l *= a.lower();
            h = high; h *= a.upper();
            }
         else
            { // a.low and a.high < 0 
            l = high; l *= a.lower();
            h = low;  h *= a.upper();
            }
      }
   else
      if( high.sign() > 0 )
         {  // low < 0, high >= 0
         if( a.lower().sign() > 0 )
            { // a.low >=0, a.high >= 0
            l = low;  l *= a.upper();
            h = high; h *= a.upper();
            }
         else
            if( a.upper().sign() > 0 )
               {  //  a.low < 0, a.high >= 0
               t.mode( ROUND_DOWN );
               l = low;  l *= a.upper(); if( l > ( t = high, t *= a.lower() ) ) l = t;
               t.mode( ROUND_UP );
               h = high; h *= a.upper(); if( h < ( t = low, t *= a.lower() ) ) h = t;
               }
            else
               { // a.low and a.high < 0 
               l = high; l *= a.lower();
               h = low;  h *= a.lower();
               }
         }
      else
         { // low and high are < 0 
         if( a.lower().sign() > 0 )
            { // a.low >=0, a.high >= 0
            l = low;  l *= a.upper();
            h = high; h *= a.lower();
            }
         else
            if( a.upper().sign() > 0 )
               {  //  a.low < 0, a.high >= 0
               l = low; l *= a.upper(); 
               h = low; h *= a.lower();
               }
            else
               { // a.low and a.high < 0 
               l = high; l *= a.upper();
               h = low; h *= a.lower();
               }
        }

   low = l;
   high = h;

   return *this;
   }

// Specilization for float and *=
// That can work with both managed an unmanged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
inline interval<float>& interval<float>::operator*=( const interval<float>& a )
   {
#ifdef HARDWARE_SUPPORT
	float l, h, t;

	if (low >= 0 ) // 
	{ // both low and high >= 0
		if (a.lower() >= 0 )
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.lower();
			fpup();
			h = high * a.upper();
		}
		else
		if (a.upper() >= 0 )
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = high * a.lower();
			fpup();
			h = high * a.upper();
		}
		else
		{ // a.low and a.high < 0
			fpdown();
			l = high * a.lower();
			fpup();
			h = low * a.upper();
		}
	}
	else
	if (high >= 0 )
	{  // low < 0, high >= 0
		if (a.lower() >= 0 )
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = high * a.upper();
		}
		else
		if (a.upper() >= 0 )
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = low * a.upper(); if (l > (t = high * a.lower())) l = t;
			fpup();
			h = high * a.upper(); if (h < (t = low * a.lower())) h = t;
		}
		else
		{ // a.low and a.high < 0 
			fpdown();
			l = high * a.lower();
			fpup();
			h = low * a.lower();
		}
	}
	else
	{ // low and high are < 0 
		if (a.lower() >= 0 )
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = high * a.lower();
		}
		else
		if (a.upper() >= 0 )
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = low * a.lower();
		}
		else
		{ // a.low and a.high < 0 
			fpdown();
			l = high * a.upper();
			fpup();
			h = low * a.lower();
		}
	}

	low = l;
	high = h;
#else
	interval<float_precision> rhs(*this), lhs(a);

	rhs *= lhs;

	low=rhs.lower();
	high=rhs.upper();
#endif
   return *this;
   }

// Specilization for double and *=
// That can work with both managed an unmanged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
inline interval<double>& interval<double>::operator*=( const interval<double>& a )
	{
#ifdef HARDWARE_SUPPORT
	double l, h, t;

	if (low >= 0) // 
	{ // both low and high >= 0
		if (a.lower() >= 0)
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.lower();
			fpup();
			h = high * a.upper();
		}
		else
		if (a.upper() >= 0)
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = high * a.lower();
			fpup();
			h = high * a.upper();
		}
		else
		{ // a.low and a.high < 0
			fpdown();
			l = high * a.lower();
			fpup();
			h = low * a.upper();
		}
	}
	else
	if (high >= 0)
	{  // low < 0, high >= 0
		if (a.lower() >= 0)
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = high * a.upper();
		}
		else
		if (a.upper() >= 0)
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = low * a.upper(); if (l > (t = high * a.lower())) l = t;
			fpup();
			h = high * a.upper(); if (h < (t = low * a.lower())) h = t;
		}
		else
		{ // a.low and a.high < 0 
			fpdown();
			l = high * a.lower();
			fpup();
			h = low * a.lower();
		}
	}
	else
	{ // low and high are < 0 
		if (a.lower() >= 0)
		{ // a.low >=0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = high * a.lower();
		}
		else
		if (a.upper() >= 0)
		{  //  a.low < 0, a.high >= 0
			fpdown();
			l = low * a.upper();
			fpup();
			h = low * a.lower();
		}
		else
		{ // a.low and a.high < 0 
			fpdown();
			l = high * a.upper();
			fpup();
			h = low * a.lower();
		}
	}

	low = l;
	high = h;
	fpnear();
#else
	interval<float_precision> rhs(*this), lhs(a);

	rhs *= lhs;
	low=rhs.lower();
	high=rhs.upper();
#endif
	return *this;
	}

// Works for all other classes
// Please note that this is for all interger classes. interval<int>, interval<long>
// were there is no loss of precision
// Actually there is specialization for both <int> and <int_precision> further down.
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator/=( const interval<_IT>& b )
   {
   interval<_IT> a, c;

   c.low = (_IT)1 / b.upper();
   c.high = (_IT)1 / b.lower();

   a = interval( low, high );
   c *= a;

   low = c.lower();
   high = c.upper();

   return *this; 
   }

// Specilization for float_precision and /=
//
inline interval<float_precision>& interval<float_precision>::operator/=( const interval<float_precision>& b )
   {
   float_precision l, h;
   interval<float_precision> c(b);

   l.precision(b.upper().precision());
   l = b.upper();
   l.mode( ROUND_DOWN );
   l = _float_precision_inverse( l );

   h.precision(b.lower().precision());
   h = b.lower();
   h.mode( ROUND_UP );
   h = _float_precision_inverse( h );

   c = interval<float_precision>( l , h );
   *this *= c;

   return *this;
   }

// Specilization for float and /=
// That can work with both managed an unmanged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
inline interval<float>& interval<float>::operator/=( const interval<float>& a )
   {
#ifdef HARDWARE_SUPPORT
	interval<float> b, c;

	fpdown();
	c.low = (float)1 / a.upper();
	fpup();
	c.high = (float)1 / a.lower();
	fpnear();

	b = interval(low, high);
	c *= b;

	low = c.lower();
	high = c.upper();
#else
	interval<float_precision> rhs(*this), lhs(a);

	rhs /= lhs;

	low=rhs.lower();
	high=rhs.upper();
#endif
   return *this;
   }

// Specilization for double and /=
// That can work with both managed an unmanged application.
// HARDWARE_SUPPORT indicate if we use the build in ROUND_DOWN or ROUND_UP Mode when performing the operation
// For now Hardware support we covert to float_precision to do the interval arithmetic in software
//
inline interval<double>& interval<double>::operator/=( const interval<double>& a )
	{
#ifdef HARDWARE_SUPPORT
	interval<double> b, c;

	fpdown();
	c.low = (double)1 / a.upper();
	fpup();
	c.high = (double)1 / a.lower();
	fpnear();

	b = interval(low, high);
	c *= b;

	low = c.lower();
	high = c.upper();
#else
	interval<float_precision> rhs(*this), lhs(a);

	rhs /= lhs;
	low=rhs.lower();
	high=rhs.upper();
#endif
	return *this;
	}

// Specilization for int_precision and /=
//
inline interval<int_precision>& interval<int_precision>::operator/=( const interval<int_precision>& b )
   {
   float_precision l = b.upper(), h = b.lower();
   interval<float_precision> c(b), a(*this);

   l.mode( ROUND_DOWN );
   l = _float_precision_inverse( l );

   h.mode( ROUND_UP );
   h = _float_precision_inverse( h );

   c = interval<float_precision>( l , h );
   a *= c;

   low = (int_precision)floor( a.lower() );
   high = (int_precision)ceil( a.upper() );

   return *this;
   }


// Specialization for int and /=
//
inline interval<int>& interval<int>::operator/=( const interval<int>& b )
   {
   double tlow, thigh;
   interval<int> a;
   interval<double> c;

   tlow = 1 / (double)b.upper();
   thigh = 1 / (double)b.lower();

   a = interval( low, high );
   c = interval<double>( tlow, thigh );
   c *= a;

   low = (int)floor( c.lower() );
   high = (int)ceil( c.upper() );

   return *this; 
   }

// Works on all classes. 
// Return the intersection
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator&=(const interval<_IT>& a)
	{
	if (a.lower() > low ) 
		low = a.lower();
	if (a.upper() < high)
		high = a.upper();
	if (low > high)  // Empty set
		{
		low = 0; high = 0;
		}

	return *this;
	}

// Works on all classes. 
// Return the union
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator|=(const interval<_IT>& a)
	{
	if (low > a.upper() || high < a.lower())
		{
		if (a.upper() - a.lower() > high - low)
			{ // return the largest set
			low = a.lower();
			high = a.upper();
			}
		}
	else
		{ // non empty intersection
		if (a.lower() < low)
			low = a.lower();
		if (a.upper() > high)
			high = a.upper();
		}
	}

	
// Works on all classes. 
// Return the set minus
//
template<class _IT> inline interval<_IT>& interval<_IT>::operator^=(const interval<_IT>& a)
	{
	if ( a.lower() < high && a.upper() > low ) // intersection is not empty
		{
		if (a.upper() <= low)
			low = a.upper();
		else
			if (a.high() >= high)
				high = a.lower();
		}

	return *this;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Essential Operators
///
//////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Binary and Unary Operators +,-,*,/
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Binary + operator
// Specialization for float_precision
//
inline interval<float_precision> operator+(const interval<float_precision>& a, const float_precision& b)
	{
	interval<float_precision> c(a);

	c += interval<float_precision>(b);
	return c;
	}

// Binary + operator
// Works for all classes
//
inline interval<float_precision> operator+(float_precision& a, const interval<float_precision>& b )
{
	interval<float_precision> c(b);

	c += interval<float_precision>(a);
	return c;
}

// Binary + operator
// Works for all classes
//
template<class _IT,class _X> inline interval<_IT> operator+(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c += interval<_IT>(_IT(b));
	return c; 
	}

// Binary + operator
// Works for all classes
//
template<class _IT,class _X> inline interval<_IT> operator+( const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(b);

	c += interval<_IT>(_IT(a));
	return c;
	}


// Binary + operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator+(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c += b;
	return c;
	}


// Unary + operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator+( const interval<_IT>& a )
   {
   return a;
   }


// Binary - operator
// Specialization for float_precision
//
inline interval<float_precision> operator-(const interval<float_precision>& a, const float_precision& b)
	{
	interval<float_precision> c(a);

	c -= interval<float_precision>(b);
	return c;
	}

// Binary - operator
// Works for all classes
//
inline interval<float_precision> operator-(float_precision& a, const interval<float_precision>& b)
	{
	interval<float_precision> c(a);

	c -= b;
	return c;
	}

// Binary - operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator-(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c -= interval<_IT>(_IT(b));
	return c;
	}

// Binary - operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator-(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c -= b;
	return c;
	}

// Binary - operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator-( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   c -= b;
   return c;
   }


// Unary - operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator-( const interval<_IT>& a )
   {
   interval<_IT> c(0);

   c -= a;
   return c;
   }

// Binary * operator
// Specialization for float_precision
//
inline interval<float_precision> operator*(const interval<float_precision>& a, const float_precision& b)
	{
	interval<float_precision> c(a);

	c *= interval<float_precision>(b);
	return c;
	}

// Binary * operator
// Works for all classes
//
inline interval<float_precision> operator*(float_precision& a, const interval<float_precision>& b)
	{
	interval<float_precision> c(b);

	c *= interval<float_precision>(a);
	return c;
	}

// Binary * operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator*(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c *= interval<_IT>(_IT(b));
	return c;
	}

// Binary * operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator*(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(b);

	c *= interval<_IT>(_IT(a));
	return c;
	}



// Binary * operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator*( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   c *= b;
   return c;
   }

// Binary / operator
// Specialization for float_precision
//
inline interval<float_precision> operator/(const interval<float_precision>& a, const float_precision& b)
	{
	interval<float_precision> c(a);

	c /= interval<float_precision>(b);
	return c;
	}

// Binary / operator
// Works for all classes
//
inline interval<float_precision> operator/(float_precision& a, const interval<float_precision>& b)
	{
	interval<float_precision> c(a);

	c /= interval<float_precision>(b);
	return c;
	}

// Binary / operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator/(const interval<_IT>& a, const _X& b)
	{
	interval<_IT> c(a);

	c /= interval<_IT>(_IT(b));
	return c;
	}

// Binary / operator
// Works for all classes
//
template<class _IT, class _X> inline interval<_IT> operator/(const _X& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c /= b;
	return c;
	}


// Binary / operator
// Works for all classes
//
template<class _IT> inline interval<_IT> operator/( const interval<_IT>& a, const interval<_IT>& b )
   {
   interval<_IT> c(a);

   if ( c == b && b.is_class() != ZERO )
	  c = interval<_IT>(1,1);
   else
      c /= b;
   
   return c;
   }

// Binary & operator
// Return intersection
// Works for all classes
//
template<class _IT> inline interval<_IT> operator&( const interval<_IT>& a, const interval<_IT>& b )
	{
	interval<_IT> c(a);

	c &= b;
	return c;
	}

// Binary | operator.
// Return union
// Works for all classes
//
template<class _IT> inline interval<_IT> operator|(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c |= b;
	return c;
	}

// Binary ^ operator
// Return set minus
// Works for all classes
//
template<class _IT> inline interval<_IT> operator^(const interval<_IT>& a, const interval<_IT>& b)
	{
	interval<_IT> c(a);

	c ^= b;
	return c;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Binary and Unary Operators
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Boolean Interval for == and !=
///
//////////////////////////////////////////////////////////////////////////////////////

template<class _IT> inline bool operator==(const interval<_IT>& a, const interval<_IT>& b)
	{
	return a.lower() == b.lower() && a.upper() == b.upper();
	}

template<class _IT> inline bool operator!=(const interval<_IT>& a, const interval<_IT>& b)
	{
	return a.lower() != b.lower() || a.upper() != b.upper();
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Boolean operators
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval abs()
///
//////////////////////////////////////////////////////////////////////////////////////

template<class _IT> inline interval<_IT> abs( const interval<_IT>& a )
	{
	if (a.lower() >= _IT(0) )
		return a;
	else
		if (a.upper() <= _IT(0) )
			return -a;

	return interval<_IT>(_IT(0), ( a.upper() > -a.lower() ? a.upper() : -a.lower() ) );
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END interval functions
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval sqrt(), log10(), log(), exp() and pow()
///
//////////////////////////////////////////////////////////////////////////////////////

// Support function for correctly converting and float_precision number back to a float
// By default IEE754 round to nearest and that will create incorrect result for interval arithmetic
//
inline float tofloat(const interval<float_precision>& fp, enum round_mode rm )
{
	float f;
	if (rm == ROUND_DOWN) // Conversion to float introduce an error since it is always round to nearest and it should be round up or round_down
		f = fp.lower();
	else
		f = fp.upper();

	if (f != 0)
		{
		float_precision fp1(f, PDIGIT10(9));
		if (rm == ROUND_DOWN && fp1 > fp.lower())
			f -= f * 0.5f * FLT_EPSILON;
		if (rm == ROUND_UP && fp1 < fp.upper())
			f += f * 0.5f * FLT_EPSILON;
		}
	return f;
}

// Support function for correctly converting and float_precision number back to a float
// By default IEE754 round to nearest and that will create incorrect result for interval arithmetic
//
inline double todouble(const interval<float_precision>& fp, enum round_mode rm)
{
	double  d;
	if (rm == ROUND_DOWN) // Conversion to float introduce an error since it is always round to nearest and it should be round up or round_down
		d = fp.lower();
	else
		d = fp.upper();

	if (d != 0)
		{
		float_precision fp1(d, PDIGIT10(17));
		if (rm == ROUND_DOWN && fp1 > fp.lower())
			d -= d * 0.5 * DBL_EPSILON;
		if (rm == ROUND_UP && fp1 < fp.upper())
			d += d * 0.5 * DBL_EPSILON;
		}
	return d;
}

#ifdef HARDWARE_SUPPORT
// Support function for correctly converting and double number back to a float
// By default IEE754 round to nearest and that will create incorrect result for interval arithmetic
//
inline float tofloat(const interval<double>& di, enum round_mode rm )
	{
	float fres;
	double d;

	switch (rm)
		{
		case ROUND_DOWN: d = di.lower(); fpdown();  break;
		case ROUND_UP: d = di.upper();  fpup();  break;
		}

	fres = (float)d;  
	fpnear();
	return fres;
	}

// Support function for correctly converting and double number back to a float
// By default IEE754 round to nearest and that will create incorrect result for interval arithmetic
//
inline float tofloat(const double& d, enum round_mode rm)
	{
	float fres;

	switch (rm)
		{
		case ROUND_DOWN:  fpdown();  break;
		case ROUND_UP: fpup();  break;
		}

	fres = (float)d;
	fpnear();
	return fres;
	}

// If H/W Support allows us to control the rounding mode then we can do it directly.
// Log(2) for double
inline double ln2double(enum round_mode rm)
	{
	double res;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fldln2;				 Load ln2
		fstp qword ptr[res]; Store result in res
		}
	fpnear();

	return res;
	}

// Log(10) for double
inline double ln10double(enum round_mode rm)
	{
	double res;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		; ln10 = FLDL2T * FLDLN2
		fldl2t;                 Load log2(10)
		fldln2;					Load LN2
		fmulp st(1),st;			Calculate LN2 * lOG2(10)
		fstp qword ptr[res];	Store ln10 in result
		}
	fpnear();

	return res;
	}

// PI for double
inline double pidouble(enum round_mode rm)
	{
	double res;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fldpi;				 Load PI
		fstp qword ptr[res]; Store result in lower
		}
	fpnear();

	return res;
	}

// Sqrt() for double
inline double sqrtdouble( double d, enum round_mode rm )
	{
	double sq = d;

	switch( rm )
	{
	case ROUND_DOWN: fpdown(); break;
	case ROUND_UP: fpup(); break;
	}

	_asm 
		{
		fld qword ptr[sq];  Load lower into floating point stack
		fsqrt;                 Calculate sqrt
		fstp qword ptr[sq]; Store result in lower
		}
	fpnear();

	return sq;
	}

// Sqrt() for float
inline float sqrtfloat(float f, enum round_mode rm)
	{
	double sq = f;
	float fres;

	sq = sqrtdouble(sq, rm);
	fres = tofloat(sq, rm);
	return fres;
	}

// log() for double
inline double logdouble(double d, enum round_mode rm)
	{
	double lg = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[lg];      Load lower into floating point stack
		fldln2;                 Load loge2
		fxch st(1);             Exchange stack top
		fyl2x;                  Calculate y * ln2 x
		fstp qword ptr[lg];     Store result in lower
		}

	fpnear();

	return lg;
	}

// log() for float
inline float logfloat(float f, enum round_mode rm)
	{
	double lg = (double)f;
	float fres;

	lg = logdouble( lg, rm );
	fres = tofloat(lg, rm);
	return fres;
	}

// log10() for double
inline double log10double(double d, enum round_mode rm)
	{
	double lg = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[lg];      Load lower into floating point stack
		fldlg2;                 Load log10(2)
		fxch st(1);             Exchange stack top
		fyl2x;                  Calculate y * ln2 x
		fstp qword ptr[lg];     Store result in lower
		}

	fpnear();

	return lg;
	}

// log10 for float
inline float log10float(float f, enum round_mode rm)
	{
	double lg = (double)f;
	float fres;

	lg = log10double( lg, rm );
	fres = tofloat(lg, rm);
	return fres;
	}

#endif

// Specilization for sqrt(float_precision)
// 
inline interval<float_precision> sqrt( const interval<float_precision>& x )
   {
   float_precision l, u;
   
   l.assign( x.lower() );  // Assign value, precision and mode
   l.mode( ROUND_DOWN );
   l = sqrt( l );

   u.assign( x.upper() );  // Assign value, precision and mode
   u.mode( ROUND_UP );
   u = sqrt( u );

   return interval<float_precision>( l, u );
   }

// sqrt for float using managed code.
//
inline interval<float> sqrt( const interval<float>& x )
   {
   float lower, upper;
#ifdef HARDWARE_SUPPORT
   lower = sqrtfloat( x.lower(), ROUND_DOWN );
   upper = sqrtfloat( x.upper(), ROUND_UP );
#else
   interval<float_precision> fx(x);

   fx=sqrt(fx);
   lower = tofloat( fx, ROUND_DOWN );
   upper = tofloat( fx, ROUND_UP );
#endif
   return interval<float>( lower, upper );
   }

// sqrt for double using managed code.
//
inline interval<double> sqrt( const interval<double>& x )
   { 
   double lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = sqrtdouble( x.lower(), ROUND_DOWN );
	upper = sqrtdouble( x.upper(), ROUND_UP );
#else
	interval<float_precision> fx(x);
   
   fx=sqrt(fx);
   lower = todouble(fx, ROUND_DOWN);
   upper = todouble(fx, ROUND_UP);
 #endif
   return interval<double>( lower, upper );
   }



// Specilization for log float_precision
// 
inline interval<float_precision> log( const interval<float_precision>& x )
   {
   float_precision l, u;
   
   l.assign( x.lower() );  // Assign value, precision and mode
   l.mode( ROUND_DOWN );
   l = log( l );

   u.assign( x.upper() );  // Assign value, precision and mode
   u.mode( ROUND_UP );
   u = log( u );

   return interval<float_precision>( l, u );
   }

// log for float using manged code.
//
inline interval<float> log( const interval<float>& x )
	{
	float lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = logfloat( x.lower(), ROUND_DOWN);
	upper = logfloat( x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx=log(fx);
	lower = tofloat( fx, ROUND_DOWN );
	upper = tofloat( fx, ROUND_UP );
#endif
   return interval<float>( lower, upper );
   }

// log for double using managed code.
//
inline interval<double> log( const interval<double>& x )
	{
	double lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = logdouble( x.lower(), ROUND_DOWN);
	upper = logdouble( x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx=log(fx);
	lower = todouble(fx, ROUND_DOWN);
	upper = todouble(fx, ROUND_UP);
#endif
   return interval<double>( lower, upper );
   }



// Specilization for log float_precision
// 
inline interval<float_precision> log10( const interval<float_precision>& x )
   {
   float_precision l, u;
   
   l.assign( x.lower() );  // Assign value, precision and mode
   l.mode( ROUND_DOWN );
   l = log10( l );

   u.assign( x.upper() );  // Assign value, precision and mode
   u.mode( ROUND_UP );
   u = log10( u );

   return interval<float_precision>( l, u );
   }



// log10 for float using manged code.
//
inline interval<float> log10( const interval<float>& x )
   {
   float lower, upper;
#ifdef HARDWARE_SUPPORT
   lower = log10float( x.lower(), ROUND_DOWN);
   upper = log10float(x.upper(), ROUND_UP);
#else
   interval<float_precision> fx(x);

   fx=log10(fx);
   lower = tofloat( fx, ROUND_DOWN );
   upper = tofloat( fx, ROUND_UP );
#endif
   return interval<float>( lower, upper );
   }

// log10 for double using managed code.
//
inline interval<double> log10( const interval<double>& x )
   {   
	double lower, upper;
#ifdef HARDWARE_SUPPORT
lower = log10double( x.lower(), ROUND_DOWN);
upper = log10double( x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);
 
	fx=log10(fx);
	lower = todouble(fx, ROUND_DOWN);
	upper = todouble(fx, ROUND_UP);
#endif
   return interval<double>( lower, upper );
   }

   
// Specilization for exp float_precision
// 
inline interval<float_precision> exp( const interval<float_precision>& x )
   {
   float_precision l, u;
   
   l.assign( x.lower() );  // Assign value, precision and mode
   l.mode( ROUND_DOWN );
   l = exp( l );

   u.assign( x.upper() );  // Assign value, precision and mode
   u.mode( ROUND_UP );
   u = exp( u );

   return interval<float_precision>( l, u );
   }

#ifdef TEST_ONLY_HVE
// exp for float using manged code.
// ONly for test purposed.
//
inline interval<float> exp2( const interval<float>& x )
   {
	interval<float_precision> fx(x);
   float lower, upper;
   
   fx=exp(fx);
   lower = tofloat(fx, ROUND_DOWN);
   upper = tofloat(fx, ROUND_UP);
   return interval<float>( lower, upper );
   }

// exp for double using managed code.
//
inline interval<double> exp2( const interval<double>& x )
   {
	interval<float_precision> fx(x);
   double lower, upper;
   
   fx=exp(fx);
   lower = todouble(fx, ROUND_DOWN);
   upper = todouble(fx, ROUND_UP);
   return interval<double>( lower, upper );
   }
#endif



// MSC exp() does not allow rounding control
// So we have to do it manually
// Use a taylor series until their is no more change in the result
// exp(x) == 1 + x + x^2/2!+x^3/3!+....
// Equivalent with the same standard C function call
// use argument reduction via exp(x)=(exp(x/2^k)2^k	
// And use Brent enhancement using the double formula:
// expm(x)=exp(x)-1 && expm(2x)=expm(x)(2+expm(x)) on the backend to preseve
// loss of significance digits
//
inline interval<double> exp(const interval<double>& x)
	{
#ifdef HARDWARE_SUPPORT
	int  i, k = 0;
	interval<double> c, res, p0, old;
	const interval<double> c1(1), c2(2);

	c = x; 
	if (x.is_class() == NEGATIVE)
		c = abs(c);

	// Determine Reduction factor
	k = int((log(2) + log(abs(c.center()))) / log(2));
	k = std::min(k, 10);
	if (k > 0)
		{
		i = 2 << (k - 1);
		c /= interval<double>( (double)i );
		}

	p0 = c;
	old = c1;
	res = old + p0;
	for (i = 2; i < 100 && (res.lower() != old.lower() && res.upper() != old.upper()); i++)
		{
		old = res;
		p0 *= ( c / interval<double>((double)i));
		res += p0;
		}
	
	// Brent enhancement avoid loss of significant digits when x is small.
	if (k>0)
		{
		res -= c1;
		for (; k > 0; k--)
			res = (c2 + res)*res;
		res += c1;
		}

	if (x.is_class() == NEGATIVE)
		res = c1 / res;
	
	return res;
#else
	interval<float_precision> fx(x);
	double lower, upper;

	fx = exp(fx);
	lower = todouble( fx, ROUND_DOWN );
	upper = todouble( fx, ROUND_UP );
	return interval<double>(lower, upper);
#endif
	}

// MSC exp() does not allow rounding control for the exp()
// Since we dont normally do it using float arithmetic (as for sqrt(), log() and log10()) we simply just call the interval<double> version
// of exp() and convert back to float preserving as much accuracy as possible
//
inline interval<float> exp(const interval<float>& x)
{
#ifdef HARDWARE_SUPPORT
	interval<double> exp(const interval<double>&);
	interval<double> fx(x);
	float lower, upper;

	fx = exp(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
	return interval<double>(lower, upper);
#else
	interval<float_precision> fx(x);
	float lower, upper;

	fx = exp(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
	return interval<float>(lower, upper);
#endif
}

// Specilization for pow float_precision
// 
inline interval<float_precision> pow( const interval<float_precision>& x, const float_precision& y )
   {
   interval<float_precision> c(x);

   c = log( x );
   c *= interval<float_precision>( y );
   c = exp( c );

   return c;
   }



// MSC pow() does not allow rounding control
// So we have to do it manually
// x^y == exp( y * ln( x ) ) );
// To avoid loss of precision we actually perform the operation using double and then 
// convert the result back to float. This is consistent with the pow() that only takes double as an argument.
// 
inline interval<float> pow(const interval<float>& x, const float y)
	{
	interval<double> c(x);
	float upper, lower;

	c = log(c);
	c *= interval<double>(y);
	c = exp(c);
	lower = (float)c.lower();
	upper = (float)c.upper();
	if (lower > c.lower() )
		lower -= lower * 0.5f * FLT_EPSILON;
	if (upper < c.upper() )
		upper += upper * 0.5f * FLT_EPSILON;

	return interval<float>(lower, upper);
	}

// MSC pow() does not alllow rounding control
// So we have to do it manually
// x^y == exp( y * ln( x ) ) );
// 
inline interval<double> pow(const interval<double>& x, const double y)
	{
	interval<double> c;

	c = log(x);
	c *= interval<double>(y);
	c = exp(c);

	return c;
	}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval sqrt(), log10(), log(), exp(), pow()
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval constants like PI, LN2 and LN10
///
//////////////////////////////////////////////////////////////////////////////////////

// Load manifest constant PI for double
//
inline interval<double> int_pidouble()
	{
	interval<double> pi;
#ifdef HARDWARE_SUPPORT
	pi.lower( pidouble(ROUND_DOWN) );
	pi.upper( pidouble(ROUND_UP) );
#else
	interval<float_precision> fx;

	fx = _float_table(_PI,20);
	pi.lower( todouble(fx,ROUND_DOWN));
	pi.upper( todouble(fx, ROUND_UP ));
#endif
	return pi;
	}

// Load manifest constant PI for float
//
inline interval<float> int_pifloat()
	{
	interval<double> pid;
	interval<float> pif;

	pid = int_pidouble();
	pif.lower(tofloat(pid.lower(), ROUND_DOWN));
	pif.upper(tofloat(pid.upper(), ROUND_UP));
	return pif;;
	}

// Specilization for constant PI float_precision
// 
inline interval<float_precision> int_pi(const unsigned int p = float_precision_ctrl.precision() )
	{
	float_precision fx(0,p+1), l(0,p), u(0,p);

	fx = _float_table(_PI, p+1);
	l.mode(ROUND_DOWN);
	l = fx;
	u.mode(ROUND_UP);
	u = fx;

	return interval<float_precision>(l, u);
	}

// Load manifest constant lN 2 for double
//
inline interval<double> int_ln2double()
	{
	interval<double> ln2;
#ifdef HARDWARE_SUPPORT
	ln2.lower(ln2double(ROUND_DOWN));
	ln2.upper(ln2double(ROUND_UP));
#else
	interval<float_precision> fx;

	fx = _float_table(_LN2, 20);
	ln2.lower(todouble(fx, ROUND_DOWN));
	ln2.upper(todouble(fx, ROUND_UP));
#endif
	return ln2;
	}

// Load manifest constant LN2 for float
//
inline interval<float> int_ln2float()
	{
	interval<double> ln2d;
	interval<float> ln2f;

	ln2d = int_ln2double();
	ln2f.lower(tofloat(ln2d.lower(), ROUND_DOWN));
	ln2f.upper(tofloat(ln2d.upper(), ROUND_UP));
	return ln2f;;
	}

// Specilization for constant LN2 float_precision
// 
inline interval<float_precision> int_ln2(const unsigned int p = float_precision_ctrl.precision())
	{
	float_precision fx(0, p + 1), l(0, p), u(0, p);

	fx = _float_table(_LN2, p + 1);
	l.mode(ROUND_DOWN);
	l = fx;
	u.mode(ROUND_UP);
	u = fx;

	return interval<float_precision>(l, u);
	}

// Load manifest constant ln10 for double
//
inline interval<double> int_ln10double()
	{
	interval<double> ln10;
#ifdef HARDWARE_SUPPORT
	ln10.lower(ln10double(ROUND_DOWN));
	ln10.upper(ln10double(ROUND_UP));
#else
	interval<float_precision> fx;

	fx = _float_table(_LN10, 20);
	ln10.lower(todouble(fx, ROUND_DOWN));
	ln10.upper(todouble(fx, ROUND_UP));
#endif
	return ln10;
	}

// Load manifest constant LN10 for float
//
inline interval<float> int_ln10float()
{
	interval<double> ln10d;
	interval<float> ln10f;

	ln10d = int_ln10double();
	ln10f.lower(tofloat(ln10d.lower(), ROUND_DOWN));
	ln10f.upper(tofloat(ln10d.upper(), ROUND_UP));
	return ln10f;;
	}

// Specilization for constant LN10 float_precision
// 
inline interval<float_precision> int_ln10(const unsigned int p = float_precision_ctrl.precision() )
	{
	float_precision fx(0, p + 2), l(0, p), u(0, p);

	fx = _float_table(_LN10, p + 2);
	l.mode(ROUND_DOWN);
	l = fx;
	u.mode(ROUND_UP);
	u = fx;

	return interval<float_precision>(l, u);
	}


//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval constants
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// Interval sin(), cos(), tan(), asin(), acos(), atan()
///
//////////////////////////////////////////////////////////////////////////////////////

#ifdef HARDWARE_SUPPORT

// Sin() for double
inline double sindouble(double d, enum round_mode rm)
	{
	double res = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[res];  Load lower into floating point stack
		fsin;				 Calculate sin
		fstp qword ptr[res]; Store result in lower
		}
	fpnear();

	return res;
	}

// Sin() for float interval
inline float sinfloat(float f, enum round_mode rm)
	{
	double res = f;
	float fres;

	res = sindouble(res, rm);
	fres = tofloat(res, rm);
	return fres;
	}


// cos() for double
inline double cosdouble(double d, enum round_mode rm)
	{
	double res = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[res];  Load lower into floating point stack
		fcos;				 Calculate cos
		fstp qword ptr[res]; Store result in lower
		}
	fpnear();

	return res;
	}

// Cos() for float interval
inline float cosfloat(float f, enum round_mode rm)
	{
	double res = f;
	float fres;

	res = cosdouble( res, rm );
	fres = tofloat(res, rm);
	return fres;
	}

// tan() for double
inline double tandouble(double d, enum round_mode rm)
	{
	double res = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[res];  Load lower into floating point stack
		fptan;               Calculate tan
		fstp qword ptr[res]; Pop ST(0) and ignore
		fstp qword ptr[res]; Store result
		}
	fpnear();

	return res;
	}

// Tan() for float interval
inline float tanfloat(float f, enum round_mode rm)
	{
	double res = f;
	float fres;

	res = tandouble( res, rm );
	fres = tofloat( res, rm );
	return fres;
	}


// tan() for double
inline double atandouble(double d, enum round_mode rm)
	{
	double res = d;

	switch (rm)
		{
		case ROUND_DOWN: fpdown(); break;
		case ROUND_UP: fpup(); break;
		}

	_asm
		{
		fld qword ptr[res];		Load lower into floating point stack
		fld1;					Load 1.0 on top of stack
		fpatan;					Calculate tan			
		fstp qword ptr[res];	Store result
		}	
	fpnear();

	return res;
	}

// Atan() for float interval
inline float atanfloat(float f, enum round_mode rm)
	{
	double res = f;
	float fres;

	res = atandouble( res, rm );
	fres = tofloat(res, rm);
	return fres;  // or alternatively return tofloat(atandouble((double)f, rm ), rm );
	}

#endif

// Specilization for sin float_precision
// 
inline interval<float_precision> sin(const interval<float_precision>& x)
{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = sin(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = sin(u);

	return interval<float_precision>(l, u);
}

// sin for float using manged code.
//
inline interval<float> sin(const interval<float>& x)
{
	float lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = sinfloat(x.lower(), ROUND_DOWN);
	upper = sinfloat(x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx = sin(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
#endif
	return interval<float>(lower, upper);
}

// Sin for double using managed code.
//
inline interval<double> sin(const interval<double>& x)
{
	double lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = sindouble(x.lower(), ROUND_DOWN);
	upper = sindouble(x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx = sin(fx);
	lower = todouble(fx, ROUND_DOWN);
	upper = todouble(fx, ROUND_UP);
#endif
	return interval<double>(lower, upper);
}

// Specilization for cos float_precision
// 
inline interval<float_precision> cos(const interval<float_precision>& x)
{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = cos(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = cos(u);

	return interval<float_precision>(l, u);
}

// cos for float using manged code.
//
inline interval<float> cos(const interval<float>& x)
{
	float lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = cosfloat(x.lower(), ROUND_DOWN);
	upper = cosfloat(x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx = cos(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
#endif
	return interval<float>(lower, upper);
}

// Cos for double using managed code.
//
inline interval<double> cos(const interval<double>& x)
{
	double lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = cosdouble(x.lower(), ROUND_DOWN);
	upper = cosdouble(x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx = cos(fx);
	lower = todouble(fx, ROUND_DOWN);
	upper = todouble(fx, ROUND_UP);
#endif
	return interval<double>(lower, upper);
}

// Specilization for tan float_precision
// 
inline interval<float_precision> tan(const interval<float_precision>& x)
{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = tan(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = tan(u);

	return interval<float_precision>(l, u);
}

// tan for float using manged code.
//
inline interval<float> tan(const interval<float>& x)
{
	float lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = tanfloat(x.lower(), ROUND_DOWN);
	upper = tanfloat(x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx = tan(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
#endif
	return interval<float>(lower, upper);
}

// Tan for double using managed code.
//
inline interval<double> tan(const interval<double>& x)
{
	double lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = tandouble(x.lower(), ROUND_DOWN);
	upper = tandouble(x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx = tan(fx);
	lower = todouble(fx, ROUND_DOWN);
	upper = todouble(fx, ROUND_UP);
#endif
	return interval<double>(lower, upper);
}

// Specilization for arctan float_precision
// 
inline interval<float_precision> atan(const interval<float_precision>& x)
{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = atan(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = atan(u);

	return interval<float_precision>(l, u);
}

// Arctan for float using manged code.
//
inline interval<float> atan(const interval<float>& x)
{
	float lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = atanfloat(x.lower(), ROUND_DOWN);
	upper = atanfloat(x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx = atan(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
#endif
	return interval<float>(lower, upper);
}



// ArcTan for double using managed code.
//
inline interval<double> atan(const interval<double>& x)
{
	double lower, upper;
#ifdef HARDWARE_SUPPORT
	lower = atandouble(x.lower(), ROUND_DOWN);
	upper = atandouble(x.upper(), ROUND_UP);
#else
	interval<float_precision> fx(x);

	fx = atan(fx);
	lower = todouble(fx, ROUND_DOWN);
	upper = todouble(fx, ROUND_UP);
#endif
	return interval<double>(lower, upper);
}

// MSC asin() does not allow rounding control
// So we have to do it manually
/// Description:
///   Use a taylor series until their is no more change in the result
///   asin(x) == x + x^3/(2*3)+(1*3)x^5/(2*4*5)+(1*3*5)x^7/(2*4*6*7)....
///   Use argument reduction via the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
//	  This function replace the other function using newton iteration. Taylor series is significant
//    faster e.g 50% for 10 digits, 3x for 100 digits and 5x for 1000 digits.
//
inline interval<double> asin(const interval<double>& x)
{
#ifdef HARDWARE_SUPPORT
	int k, sign;
	interval<double> r, u, v, v2, sqrt2, lc, uc;
	const double c1(1), c2(2);

	if (x.lower() >= c1 || x.upper() <= -c1)
		{
		throw interval<double>::domain_error(); return x;
		}

	v = x;
	if (v.lower() < -c1)
		v.lower(-c1);
	if (v.upper() > c1)
		v.upper(c1);

	sign = v.is_class();
	if (sign == NEGATIVE)
		v = -v;

	// Now use the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
	// until argument is less than dlimit
	// Reduce the argument to below 0.5 to make the newton run faster
	sqrt2 = interval<double>(c2);				// Ensure correct number of digits
	sqrt2 = sqrt(sqrt2);	// Now calculate sqrt2 with precision digits
	for (k = 0; v.lower() > 0.5; k++)
		v /= sqrt2 * sqrt(interval<double>(c1) + sqrt(interval<double>(c1) - v * v));

	v2 = v * v;
	r = v;
	u = v;
	// Now iterate using taylor expansion
	for (unsigned int j = 3;; j += 2)
		{
		uc = interval<double>((j - 2) * (j - 2));
		lc = interval<double>(j * j - j);
		v = uc * v2 / lc;
		r *= v;
		if( u.lower() + r.lower() == u.lower() || u.upper() + r.upper() == u.upper())
			break;
		u += r;
		}

	if (k > 0)
		u *= interval<double>(1 << k);

	if (sign == NEGATIVE )
		u = -u;

	return u;
#else
	interval<float_precision> fx(x);
	double lower, upper;

	fx = asin(fx);
	lower = todouble(fx, ROUND_DOWN);
	upper = todouble(fx, ROUND_UP);
	return interval<double>(lower, upper);
#endif
}

// MSC acos() does not allow rounding control
// So we have to do it manually
/// Description:
///   Use a taylor series until their is no more change in the result
///   asin(x) == x + x^3/(2*3)+(1*3)x^5/(2*4*5)+(1*3*5)x^7/(2*4*6*7)....
///   Use argument reduction via the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
//	  This function replace the other function using newton iteration. Taylor series is significant
//    faster e.g 50% for 10 digits, 3x for 100 digits and 5x for 1000 digits.
//
inline interval<double> acos(const interval<double>& x)
{
#ifdef HARDWARE_SUPPORT
	interval<double> pi, res;
	const double c1(1);

	if (x.lower() >= c1 || x.upper() <= -c1)
		{
		throw interval<double>::domain_error(); return x;
		}

	pi.lower( pidouble(ROUND_DOWN) );
	pi.upper( pidouble(ROUND_UP) );
	res = pi * interval<double>(0.5) - asin( x );
	return res;
#else
	interval<float_precision> fx(x);
	double lower, upper;

	fx = acos(fx);
	lower = todouble(fx, ROUND_DOWN);
	upper = todouble(fx, ROUND_UP);
	return interval<double>(lower, upper);
#endif
}

// MSC asin() does not allow rounding control for the asin()
// Since we dont normally do it using float arithmetic (as for sqrt(), log() and log10()) we simply just call the interval<double> version
// of asin() and convert back to float preserving as much accuracy as possible
//
inline interval<float> asin(const interval<float>& x)
{
#ifdef HARDWARE_SUPPORT
	interval<double> asin(const interval<double>&);
	interval<double> fx(x);
	float lower, upper;

	fx = asin(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
	return interval<double>(lower, upper);
#else
	interval<float_precision> fx(x);
	float lower, upper;

	fx = asin(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
	return interval<float>(lower, upper);
#endif
}

// Specilization for arcsin float_precision
// 
inline interval<float_precision> asin(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = asin(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = asin(u);

	return interval<float_precision>(l, u);
	}

// MSC acos() does not allow rounding control for the acos()
// Since we dont normally do it using float arithmetic (as for sqrt(), log() and log10()) we simply just call the interval<double> version
// of acos() and convert back to float preserving as much accuracy as possible
//
inline interval<float> acos(const interval<float>& x)
	{
#ifdef HARDWARE_SUPPORT
	interval<double> acos(const interval<double>&);
	interval<double> fx(x);
	float lower, upper;

	fx = acos(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
	return interval<double>(lower, upper);
#else
	interval<float_precision> fx(x);
	float lower, upper;

	fx = acos(fx);
	lower = tofloat(fx, ROUND_DOWN);
	upper = tofloat(fx, ROUND_UP);
	return interval<float>(lower, upper);
#endif
	}

// Specilization for arcsin float_precision
// 
inline interval<float_precision> acos(const interval<float_precision>& x)
	{
	float_precision l, u;

	l.assign(x.lower());  // Assign value, precision and mode
	l.mode(ROUND_DOWN);
	l = acos(l);

	u.assign(x.upper());  // Assign value, precision and mode
	u.mode(ROUND_UP);
	u = acos(u);

	return interval<float_precision>(l, u);
	}


//////////////////////////////////////////////////////////////////////////////////////
///
/// END Interval sin(), cos(), tan(), asin(), acos(), atan()
///
//////////////////////////////////////////////////////////////////////////////////////
}
#endif
