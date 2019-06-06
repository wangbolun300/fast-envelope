#ifndef INC_FPRECISION
#define INC_FPRECISION

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002-2017
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
 * Module name     :   fprecision.h
 * Module ID Nbr   :   
 * Description     :   Arbitrary floating poiint precision class
 * --------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version	Author/Date		Description of changes
 * -------  -----------		----------------------
 * 01.01	HVE/020209		Initial release
 * 01.02    HVE/030421     Use FFT for Fast multiplication. Clean up the code
 *                         Now it can be included in different compilation units
 *                         Without causing linker errors.
 * 01.03    HVE/030512     Clean up _float_precision_ftoa()
 * 01.04    HVE/030701     A bug in internal base 256 when initialize a 
 *                         float_precision with a double value has been fixed
 *                         The accuracy and speed of the trignometric function
 *                         tan(), sin() & cos() has been greatly improved.
 * 01.05    HVE/030715     Use trisection identity for cos instead of bisection
 *                         identity
 * 01.06    HVE/050121     Doxygen comments added
 * 01.07    HVE/050311     Replacing #define with inline functions
 * 01.08    HVE/060127     A precision bug in /= was fixed. The division was erroneous
 *                         done in fixed PRECISION instead of the maximum of both sides
 * 01.09    HVE/060203     Declaration of atan2 was missing
 * 01.10    HVE/060203     _gcvt replaced with sprintf("%.18g"). _gcvt was Microsoft specific
 *                         Automatically detect big indian/little indian binary floating point format
 *                         Only a problem for the double constructor in RADIX 256 mode
 * 01.11    HVE/060205     A bug in the operator int conversion of float_precision to int has been fixed
 * 01.12    hve/060305     Added construction for all integers, long,short,int,char and unsigned version
 * 01.13	HVE/Jul312007	Sometimes a float precision variable constructor call was initialized prior to 
 *									the constructor call of float_precision_ctrl making the float precision variable
 *									equal to 0. which of course make no sense. The bug has been fixed
 * 01.14	HVE/Jun25 2007		== Comparison with zero fixed a bug that treated +0E+3==+0E0 as different
 *									Also fixed a bug that in the constructor that created 0 as 0E-1 instead of just 0E0
 * 01.15	HVE/Jul 10 2010	Fix a bug when using BASE 256 as the internal representation
 * 01.16	HVE/Jul 11 2010	Fix bugs in F_RADIX=BASE_2 
 * 01.17	HVE/Aug 17 2012 Add ftoainteger to convert a float to a integer ascii string
 * 01.18	HVE/JUN-12-2012	Use umul_short() for single digit multiplication instead of umult_fourier() in *= operator
 * 01.19	HVE/MAR-07-2014 Moved <limits.h> to iprecision.h to ensure it is included even when only working with arbitrary integer
 *							precision and consequencly only use iprecision.h
 * 01.20	HVE/MAY-03-2014	Made the header Visual studio C++ 2013 compatible (avoiding error C2535)
 * 01.21	HVE/JUN-17-2014 Polishing the code for use with other compilers
 * 01.22	HVE/JUN-22-2014	Remove the requirements for microsoft stdafx.h precompiled header plus make (int_precision) available as a 
 *							method .to_int_precision() to avoid ambiguity between different compiles that can't properly resolved it.
 * 01.23	HVE/JUL-6-2014	Replaced the old C Macros MIN(), MAX() with std::min() and std::max()
 * 01.24	HVE/JUL-28-2014	Added "override" of <oper>(int_precision,other type) in iprecision.h  to float_precision operator<oper>(int_precision, float_precision)	in fprecision.h
 * 01.25	HVE/JUL-29-2014 Added template <class _Ty> float_precision operator<oper>(float_precision, const _Ty) and
 *							template <class _Ty> float_precision operator<oper>(const _Ty, const float_precision) for oper -,*,/. "+" was already implemented
 * 01.26	HVE/AUG-19-2016	Added initialization of fprecision through a std::string
 * 01.27	HVE/NOV-07-2016	Added initialization of fprecision through a signed/unsinged int64_t 64bit data  type
 * 01.28	HVE/NOV-13-2016	Added nrooth function declaration
 * 01.29	HVE/JAN-29-2016	Added support for table lookup of exp(1) as a transcendntal constant
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/


/* define version string */
static char _VF_[] = "@(#)fprecision.h 01.29 -- Copyright (C) Future Team Aps";

#include <algorithm>
#include "iprecision.h"


/// The four different ronding modes
/// # ROUND_NEAR  Rounded result is the closest to the infinitely precise result.
/// # ROUND_DOWN  Rounded result is close to but no greater than the infinitely precise result.
/// # ROUND_UP    Rounded result is close to but no less than the infinitely precise result.
/// # ROUND ZERO  Rounded result is close to but no greater in absolute value than the infinitely precise result.
enum round_mode { ROUND_NEAR, ROUND_UP, ROUND_DOWN, ROUND_ZERO };

// The build in constant!
enum table_type { _LN2, _LN10, _PI, _EXP1 };

// Default precision of 20 Radix digits if not specified
static const int PRECISION = 20;

// Float_precision radix. Can be either BASE 2, BASE_10, BASE 16 or BASE_256
static const int F_RADIX = BASE_10; 

inline unsigned char FDIGIT( char x )			{ return F_RADIX <= 10 ? (unsigned char)( x - '0') : (unsigned char)x; }
inline unsigned char FCHARACTER( char x )		{ return F_RADIX <= 10 ? (unsigned char)( x + '0') : (unsigned char)x; }
inline unsigned char FCHARACTER10( char x)		{ return (unsigned char)( x + '0'); }

inline int FCARRY( unsigned int x )				{ return x / F_RADIX; }
inline int FSINGLE( unsigned int x )			{ return x % F_RADIX; }

inline double PLOG10(unsigned int x )			{ return log( x * log( (double)F_RADIX ) / log( (double)BASE_10) ) / log((double)BASE_10); }
inline unsigned int PADJUST( unsigned int x )	{ return (int)( x * ( log( (double)BASE_10 ) / log( (double)F_RADIX ) ) + 0.9 ); }
inline unsigned int PDIGIT10(unsigned int digit)	{ return (unsigned int)ceil( ( digit * log( (double)F_RADIX ) ) /log( (double)BASE_10 ) ); }


///
/// @class float_precision_ctrl
/// @author Henrik Vestermark (hve@hvks.com)
/// @date  2/15/2006
/// @version 1.0
/// @brief  Arbitrary float precision control class
///
/// @todo
///
///// Float Precision control class
///   This keep track of the global settings of default precision and round mode.
///   Everytime a new float_precision constructor is invoked it takes the default 
///   precision and round mode from this float_precision_ctrl class. Unless a precision and/or
///   rounding mode has explicit been specified.
///   Default preicsion is the manifest constant PRECISION 
///   Default rounding mode is ROUND_NEAR (round to nearest)
//
class float_precision_ctrl {
   enum round_mode   mRmode;  // Global Rounding mode. Default Round Nearest
   unsigned int      mPrec;   // Global Number of decimals in mantissa. Default PRECISION.
   
   public:
      // Constructor
      float_precision_ctrl( unsigned int p=PRECISION, enum round_mode rm=ROUND_NEAR ): mRmode(rm), mPrec(p) {}
      
      // Coordinate functions
      enum round_mode mode() const                 { return mRmode; }
      enum round_mode mode( enum round_mode m )    { return( mRmode = m ); }
      unsigned int precision() const               { return mPrec>0 ? mPrec : PRECISION; }
      unsigned precision( unsigned int p )         { mPrec = p > 0 ? p : PRECISION; return mPrec; }
   };

extern float_precision_ctrl float_precision_ctrl;

class float_precision;
   
// Arithmetic + Binary and Unary
template <class _Ty> inline float_precision operator+( float_precision&, const _Ty& );
template <class _Ty> inline float_precision operator+( const _Ty&, const float_precision& );
inline float_precision operator+( int_precision&, float_precision& );					// Override int_precision - other type in iprecision.h
//inline float_precision operator+( float_precision&,int_precision&);					// Override int_precision - other type in iprecision.h
inline float_precision operator+( const float_precision& );								// Unary 
//float_precision operator+( float_precision&, float_precision& );			// Binary. Obsolete

// Arithmetic - Binary and Unary
template <class _Ty> inline float_precision operator-(float_precision&, const _Ty&);
template <class _Ty> inline float_precision operator-(const _Ty&, const float_precision&);
inline float_precision operator-(int_precision&, float_precision&);						// Override int_precision - other type in iprecision.h
//inline float_precision operator-( const float_precision&, const float_precision& );	// Binary. Obsolete
inline float_precision operator-( const float_precision& );								// Unary

// Arithmetic * Binary
template <class _Ty> inline float_precision operator*(float_precision&, const _Ty&);
template <class _Ty> inline float_precision operator*(const _Ty&, const float_precision&);
inline float_precision operator*(int_precision&, float_precision&);						// Override int_precision * other type in iprecision.h
//inline static float_precision operator*( const float_precision&, const float_precision& );  // Binary. Obsolete

// Arithmetic / Binary
template <class _Ty> inline float_precision operator/(float_precision&, const _Ty&);
template <class _Ty> inline float_precision operator/(const _Ty&, const float_precision&);
inline float_precision operator/(int_precision&, float_precision&);						// Override int_precision / other type in iprecision.h
//inline float_precision operator/( const float_precision&, const float_precision& );		// Binary. Obsolete

// Boolean Comparision Operators
inline bool operator> ( const float_precision&, const float_precision& );
inline bool operator< ( const float_precision&, const float_precision& );
inline bool operator==( const float_precision&, const float_precision& );
inline bool operator!=( const float_precision&, const float_precision& );
inline bool operator>=( const float_precision&, const float_precision& );
inline bool operator<=( const float_precision&, const float_precision& );

// Precision Floating point functions equivalent with the std C functions
extern float_precision modf( const float_precision&, float_precision * );
extern float_precision fmod( const float_precision&, const float_precision& );
extern float_precision floor( const float_precision& );
extern float_precision ceil( const float_precision& );
extern float_precision fabs( const float_precision& );  // Obsolete. replaced by overloaded abs(). But here for backward compatitbility
extern float_precision sqrt( const float_precision& );
extern float_precision log10( const float_precision& );
extern float_precision log( const float_precision& );
extern float_precision exp( const float_precision& );
extern float_precision pow( const float_precision&, const float_precision& );
extern float_precision frexp( float_precision&, int * );
extern float_precision ldexp( const float_precision&, int );
extern float_precision abs( const float_precision& );

// Trigonometric functions
extern float_precision atan( const float_precision& );
extern float_precision atan2( const float_precision&, const float_precision& );
extern float_precision asin( const float_precision& );
extern float_precision acos( const float_precision& );
extern float_precision sin( const float_precision& );
extern float_precision cos( const float_precision& );
extern float_precision tan( const float_precision& );

// Hyperbolic functions
extern float_precision sinh( const float_precision& );
extern float_precision cosh( const float_precision& );
extern float_precision tanh( const float_precision& );
extern float_precision asinh( const float_precision& );
extern float_precision acosh( const float_precision& );
extern float_precision atanh( const float_precision& ); 

// Miscelanneous support function
extern float_precision nroot(const float_precision&, unsigned int); 

// Support functions. Works on float_precision 
float_precision _float_precision_inverse( const float_precision& );
float_precision _float_table( enum table_type, unsigned int );
std::string _float_precision_ftoa( const float_precision * );
std::string _float_precision_ftoainteger( const float_precision * );
float_precision _float_precision_atof( const char *, unsigned int, enum round_mode );
float_precision _float_precision_dtof( double, unsigned int, enum round_mode );

// Core Supporting functions. Works directly on string class
int _float_precision_normalize( std::string * );
int _float_precision_rounding( std::string *, int, unsigned int, enum round_mode );
void _float_precision_strip_leading_zeros( std::string * );
void _float_precision_strip_trailing_zeros( std::string * );
void _float_precision_right_shift( std::string *, int );
void _float_precision_left_shift( std::string *, int );
int _float_precision_compare( std::string *, std::string * );
std::string _float_precision_uadd_short( std::string *, unsigned int );
std::string _float_precision_uadd( std::string *, std::string * );
std::string _float_precision_usub_short( int *, std::string *, unsigned int );
std::string _float_precision_usub( int *, std::string *, std::string * );
std::string _float_precision_umul_short( std::string *, unsigned int );
std::string _float_precision_umul( std::string *, std::string * );
std::string _float_precision_umul_fourier( std::string *, std::string * );
std::string _float_precision_udiv_short( unsigned int *, std::string *, unsigned int );
std::string _float_precision_udiv( std::string *, std::string * );
std::string _float_precision_urem( std::string *, std::string * );

///
/// @class float_precision
/// @author Henrik Vestermark (hve@hvks.com)
/// @date  1/24/2005
/// @version 1.0
/// @brief  Arbitrary float precision class
///
/// @todo
///
///// Float Precision class
///   An Arbitrary float always has the format [sign][digit][.[digit]*][E[sign][digits]+] where sign is either '+' or '-'
///   And is always stored in normalized mode after an operation or conversion
///   The length or the representation is always >= 2
///   A null string is considered as an error and an exception is thrown
///   Floating Point Numbers is stored in BASE 10 or BASE 256 depends of the F_RADIX setting
///   Also number is always strip for leading nosignificant zeros
//
class float_precision {
   enum round_mode   mRmode;  // Rounding mode. Default Round Nearest
   unsigned int      mPrec;   // Number of decimals in mantissa. Default 20, We make a shot cut by assuming the number of digits can't exceed 2^32-1
   int               mExpo;   // Exponent as a power of RADIX (not 2 as in IEEE 754). We make a short cut here and use a standard int to hold 
                              // the exponent. This will allow us exponent in the range from -RADIX^2^31 to  RADIX^2^31. Which should be enough
   std::string       mNumber; // The mantissa any length however the fraction point is always after the first digit and is implied

   public:
      // Constructors
	  float_precision()							{ mRmode = float_precision_ctrl.mode();
												mPrec = float_precision_ctrl.precision();
												mExpo = 0;
												mNumber = "+" + FCHARACTER(0);            // Build number
												}
      float_precision( char, unsigned int, enum round_mode );				// When initialized through a char
      float_precision( unsigned char, unsigned int, enum round_mode );		// When initialized through a unsigned char
      float_precision( short, unsigned int, enum round_mode );				// When initialized through a short
      float_precision( unsigned short, unsigned int, enum round_mode );		// When initialized through a unsigned short
      float_precision( int, unsigned int, enum round_mode );				// When initialized through a int
      float_precision( unsigned int, unsigned int, enum round_mode );		// When initialized through a unsigned int
      float_precision( long, unsigned int, enum round_mode );				// When initialized through a long
      float_precision( unsigned long, unsigned int, enum round_mode );		// When initialized through a unsigned long
	  float_precision( int64_t, unsigned int, enum round_mode);				// When initialized through a int64_t
	  float_precision( uint64_t, unsigned int, enum round_mode);			// When initialized through a uint64_t
      float_precision( double, unsigned int, enum round_mode );				// When initialized through a double
      float_precision( const char *, unsigned int, enum round_mode );		// When initialized through a char string
	  float_precision( const std::string&, unsigned int, enum round_mode);	// When initialized through a std::string
      float_precision( const float_precision& s /*= float_precision(0, float_precision_ctrl.precision(), float_precision_ctrl.mode() )*/ ): mNumber(s.mNumber), mRmode(s.mRmode), mPrec(s.mPrec), mExpo(s.mExpo) {}  // When initialized through another float_precision
      float_precision( const int_precision&, unsigned int, enum round_mode ); 

      // Coordinate functions
      std::string get_mantissa() const             { return mNumber; };    // Copy of mantissa
      std::string *ref_mantissa()                  { return &mNumber; }    // Reference of Mantissa
      enum round_mode mode() const                 { return mRmode; }
      enum round_mode mode( enum round_mode m )    { return( mRmode = m ); }
      int exponent() const                         { return mExpo; };
      int exponent( int e )                        { return( mExpo = e ); }  
      int sign() const                             { return CHAR_SIGN( mNumber[0] ); }
      unsigned int precision() const               { return mPrec; }
      unsigned precision( unsigned int p )         { int sign; std::string m;
                                                   sign = CHAR_SIGN( mNumber[0] );
                                                   mPrec = p > 0 ? p : float_precision_ctrl.precision();
                                                   m = mNumber.substr(1); // Bypass sign
                                                   mExpo += _float_precision_rounding( &m, sign, mPrec, mRmode );
                                                   mNumber = SIGN_STRING( sign ) + m;
                                                   return mPrec;
                                                   }
       float_precision epsilon();							// Return Beta^(1-t)
      
      void set_n( std::string mantissa )			{ mNumber = mantissa; }    // Secret function
      float_precision assign( const float_precision& a )  { mRmode = a.mRmode; mPrec = a.mPrec; mExpo = a.mExpo; mNumber = a.mNumber; return *this; }


      int change_sign()   { // Change and return sign   
                          if( mNumber.length() != 2 || FDIGIT( mNumber[1] ) != 0 ) // Don't change sign for +0!
                             if( mNumber[0] == '+' ) mNumber[0] = '-'; else mNumber[0] = '+';
                          return CHAR_SIGN( mNumber[0] );
                          }

	  // Conversion methods. Safer and less ambiguios than overloading implicit/explivit conversion operators
	  std::string toString() const					{ return _float_precision_ftoa(this); }
	  int_precision to_int_precision() const		{ std::string s = _float_precision_ftoainteger(this); return (int_precision)((char *)s.c_str()); }

	  // Implict/explicit conversion operators
      operator char() const;
      operator short() const;
      operator int() const;
      operator long() const;
      operator unsigned char() const;
      operator unsigned short() const;
      operator unsigned int() const;
      operator unsigned long() const;
      operator double() const;
      operator float() const;
      operator int_precision() const;

      // Essential operators
      float_precision& operator= ( const float_precision& );
      float_precision& operator+=( const float_precision& );
      float_precision& operator-=( const float_precision& );
      float_precision& operator*=( const float_precision& );
      float_precision& operator/=( const float_precision& );

      // Specialization
      friend std::ostream& operator<<( std::ostream& strm, const float_precision& d );
      friend std::istream& operator>>( std::istream& strm, float_precision& d );

      // Exception class
      class bad_int_syntax		{}; 
      class bad_float_syntax	{};
      class out_of_range		{};
      class divide_by_zero		{};
      class domain_error		{};
      class base_error			{};
   };




//////////////////////////////////////////////////////////////////////////////////////
///
/// FLOAT PRECISION CONSTRUCTORS
///
//////////////////////////////////////////////////////////////////////////////////////


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 		Char constructor for float_precision
///	@return 		nothing
///	@param      "c"	-	Integer char digit
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initilize with a character
///   Input Always in BASE_10
//
inline float_precision::float_precision( const char c, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   int sign;
   std::string number;

   if( c < '0' || c > '9' )
      throw bad_int_syntax(); 
   else
      {
      mExpo = 0;
      mRmode = m;
      mPrec = p;

      number = ito_precision_string( IDIGIT10( c ), true, F_RADIX );  // Convert to integer
      sign = CHAR_SIGN( number[0] );;                    // Get Sign
      number = number.substr( 1 );                       // Remove sign
      _float_precision_strip_leading_zeros( &number );   // First strip for leading zeros
      mExpo = number.length() -1;                        // Always one digit before the dot
      _float_precision_strip_trailing_zeros( &number );  // Get rid of trailing non-significant zeros
      mExpo += _float_precision_rounding( &number, sign, mPrec, mRmode );    // Perform rounding
      mNumber = SIGN_STRING( sign ) + number;            // Build number
      }
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/05/2006
///	@brief 		Unsigned Char constructor for float_precision
///	@return 		nothing
///	@param      "c"	-	Integer unsigned char digit
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initilize with a character
///   Input Always in BASE_10
//
inline float_precision::float_precision( const unsigned char c, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   int sign;
   std::string number;

   if( c < '0' || c > '9' )
      throw bad_int_syntax(); 
   else
      {
      mExpo = 0;
      mRmode = m;
      mPrec = p;

      number = ito_precision_string( IDIGIT10( c ), true, F_RADIX );  // Convert to integer
      sign = CHAR_SIGN( number[0] );;                    // Get Sign
      number = number.substr( 1 );                       // Remove sign
      _float_precision_strip_leading_zeros( &number );   // First strip for leading zeros
      mExpo = number.length() -1;                        // Always one digit before the dot
      _float_precision_strip_trailing_zeros( &number );  // Get rid of trailing non-significant zeros
      mExpo += _float_precision_rounding( &number, sign, mPrec, mRmode );    // Perform rounding
      mNumber = SIGN_STRING( sign ) + number;            // Build number
      }
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  03/05/2006
///	@brief 		Short constructor for float_precision
///	@return 		nothing
///	@param      "i"	-	Integer number
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initialize with integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline float_precision::float_precision( short i, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   int sign;
   std::string number;

   mRmode = m;
   mPrec = p;
   mExpo = 0;

   number = ito_precision_string( i, true, F_RADIX );          // Convert to integer
   sign = CHAR_SIGN( number[0] );;                    // Get Sign
   number = number.substr( 1 );                       // Remove sign
   _float_precision_strip_leading_zeros( &number );   // First strip for leading zeros
   mExpo = number.length() -1;                        // Always one digit before the dot
   _float_precision_strip_trailing_zeros( &number );  // Get rid of trailing non-significant zeros
   mExpo += _float_precision_rounding( &number, sign, mPrec, mRmode );    // Perform rounding
   mNumber = SIGN_STRING( sign ) + number;            // Build number
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/03/2006
///	@brief 		Unsigned Short constructor for float_precision
///	@return 		nothing
///	@param      "i"	-	Integer number
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initialize with integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline float_precision::float_precision( unsigned short i, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   int sign;
   std::string number;

   mRmode = m;
   mPrec = p;
   mExpo = 0;

   number = ito_precision_string( i, false, F_RADIX );         // Convert to integer
   sign = CHAR_SIGN( number[0] );;                    // Get Sign
   number = number.substr( 1 );                       // Remove sign
   _float_precision_strip_leading_zeros( &number );   // First strip for leading zeros
   mExpo = number.length() -1;                        // Always one digit before the dot
   _float_precision_strip_trailing_zeros( &number );  // Get rid of trailing non-significant zeros
   mExpo += _float_precision_rounding( &number, sign, mPrec, mRmode );    // Perform rounding
   mNumber = SIGN_STRING( sign ) + number;            // Build number
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 		Integer constructor for float_precision
///	@return 		nothing
///	@param      "i"	-	Integer number
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initialize with integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline float_precision::float_precision( int i, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   int sign;
   std::string number;

   mRmode = m;
   mPrec = p;
   mExpo = 0;

   number = ito_precision_string( i, true, F_RADIX );          // Convert to integer
   sign = CHAR_SIGN( number[0] );;                    // Get Sign
   number = number.substr( 1 );                       // Remove sign
   _float_precision_strip_leading_zeros( &number );   // First strip for leading zeros
   mExpo = number.length() -1;                        // Always one digit before the dot
   _float_precision_strip_trailing_zeros( &number );  // Get rid of trailing non-significant zeros
   mExpo += _float_precision_rounding( &number, sign, mPrec, mRmode );    // Perform rounding
   mNumber = SIGN_STRING( sign ) + number;            // Build number
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/03/2006
///	@brief 		Unsigned Integer constructor for float_precision
///	@return 		nothing
///	@param      "i"	-	Integer number
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initialize with integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline float_precision::float_precision( unsigned int i, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   int sign;
   std::string number;

   mRmode = m;
   mPrec = p;
   mExpo = 0;

   number = ito_precision_string( i, false, F_RADIX );         // Convert to integer
   sign = CHAR_SIGN( number[0] );;                    // Get Sign
   number = number.substr( 1 );                       // Remove sign
   _float_precision_strip_leading_zeros( &number );   // First strip for leading zeros
   mExpo = number.length() -1;                        // Always one digit before the dot
   _float_precision_strip_trailing_zeros( &number );  // Get rid of trailing non-significant zeros
   mExpo += _float_precision_rounding( &number, sign, mPrec, mRmode );    // Perform rounding
   mNumber = SIGN_STRING( sign ) + number;            // Build number
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/05/2006
///	@brief 		long constructor for float_precision
///	@return 		nothing
///	@param      "i"	-	Integer number
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initialize with integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline float_precision::float_precision( long i, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   int sign;
   std::string number;

   mRmode = m;
   mPrec = p;
   mExpo = 0;

   number = ito_precision_string( i, true, F_RADIX );          // Convert to integer
   sign = CHAR_SIGN( number[0] );;                    // Get Sign
   number = number.substr( 1 );                       // Remove sign
   _float_precision_strip_leading_zeros( &number );   // First strip for leading zeros
   mExpo = number.length() -1;                        // Always one digit before the dot
   _float_precision_strip_trailing_zeros( &number );  // Get rid of trailing non-significant zeros
   mExpo += _float_precision_rounding( &number, sign, mPrec, mRmode );    // Perform rounding
   mNumber = SIGN_STRING( sign ) + number;            // Build number
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/03/2006
///	@brief 		Unsigned Long constructor for float_precision
///	@return 		nothing
///	@param      "i"	-	Integer number
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initialize with integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline float_precision::float_precision( unsigned long i, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   int sign;
   std::string number;

   mRmode = m;
   mPrec = p;
   mExpo = 0;

   number = ito_precision_string( i, false, F_RADIX );         // Convert to integer
   sign = CHAR_SIGN( number[0] );;                    // Get Sign
   number = number.substr( 1 );                       // Remove sign
   _float_precision_strip_leading_zeros( &number );   // First strip for leading zeros
   mExpo = number.length() -1;                        // Always one digit before the dot
   _float_precision_strip_trailing_zeros( &number );  // Get rid of trailing non-significant zeros
   mExpo += _float_precision_rounding( &number, sign, mPrec, mRmode );    // Perform rounding
   mNumber = SIGN_STRING( sign ) + number;            // Build number
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/NOV/2016
///	@brief 		int64_t 64bit constructor for float_precision
///	@return 		nothing
///	@param      "i"	-	Integer number
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initialize with integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline float_precision::float_precision( int64_t i, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode())
	{
	int sign;
	std::string number;

	mRmode = m;
	mPrec = p;
	mExpo = 0;

	number = i64to_precision_string(i, true, F_RADIX);          // Convert to integer
	sign = CHAR_SIGN(number[0]);;                    // Get Sign
	number = number.substr(1);                       // Remove sign
	_float_precision_strip_leading_zeros(&number);   // First strip for leading zeros
	mExpo = number.length() - 1;                        // Always one digit before the dot
	_float_precision_strip_trailing_zeros(&number);  // Get rid of trailing non-significant zeros
	mExpo += _float_precision_rounding(&number, sign, mPrec, mRmode);    // Perform rounding
	mNumber = SIGN_STRING(sign) + number;            // Build number
	}


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/NOV/2016
///	@brief 		Usinged 64bit constructor for float_precision
///	@return 		nothing
///	@param      "i"	-	64bit Integer number
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate and initialize with integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline float_precision::float_precision( uint64_t i, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode())
	{
	int sign;
	std::string number;

	mRmode = m;
	mPrec = p;
	mExpo = 0;

	number = i64to_precision_string(i, false, F_RADIX);         // Convert to integer
	sign = CHAR_SIGN(number[0]);;                    // Get Sign
	number = number.substr(1);                       // Remove sign
	_float_precision_strip_leading_zeros(&number);   // First strip for leading zeros
	mExpo = number.length() - 1;                        // Always one digit before the dot
	_float_precision_strip_trailing_zeros(&number);  // Get rid of trailing non-significant zeros
	mExpo += _float_precision_rounding(&number, sign, mPrec, mRmode);    // Perform rounding
	mNumber = SIGN_STRING(sign) + number;            // Build number
	}

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 		string constructor for float_precision
///	@return 		nothing
///	@param      "str"	-	Floating point number as a string
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate input and convert to internal representation
///   Always add sign if not specified 
///   Only use core base functions to create multi precision numbers
///   The float can be any integer or decimal float representation
//
inline float_precision::float_precision( const char *str, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   if( str == NULL || *str == '\0' )
      { throw bad_int_syntax(); return; }

   mRmode = m;
   mPrec = p;
   *this = _float_precision_atof( str, p, m );
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  Aug-19-2016
///	@brief 		std::string constructor for float_precision
///	@return 		nothing
///	@param      "str"	-	Floating point number as a std::string
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor
///   Validate input and convert to internal representation
///   Always add sign if not specified 
///   Only use core base functions to create multi precision numbers
///   The float can be any integer or decimal float representation
//

inline float_precision::float_precision(const std::string& str, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode())
	{
	if ( str.empty() )
		{
		throw bad_int_syntax(); return;
		}

	mRmode = m;
	mPrec = p;
	*this = _float_precision_atof( str.c_str(), p, m);
	}


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 		double constructor for float_precision
///	@return 		nothing
///	@param      "d"	-	Floating point number in IEE754 format
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor for double
///   Validate input and convert to internal representation
///   Always add sign if not specified 
///   Only use core base functions to create multi precision numbers
///   The float can be any integer or decimal float representation
//
inline float_precision::float_precision( double d, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   mRmode = m;
   mPrec = p;
   *this = _float_precision_dtof( d, p, m );
   }



///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 		double constructor for float_precision
///	@return 		nothing
///	@param      "ip"	-	arbitrary integer precision 
///	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
///	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
///
///	@todo 
///
/// Description:
///   Constructor for int_precision to float_precision
///   Slow impl.
///    1) conver to Ascii decimal string and then 
///   2) convert it back to floating format
//
inline float_precision::float_precision( const int_precision& ip, unsigned int p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   std::string s;

   s = _int_precision_itoa( const_cast<int_precision*>(&ip) );
   mRmode = m;
   mPrec = p;
   *this = float_precision( (char *)s.c_str(), p, m );
   }

//////////////////////////////////////////////////////////////////////////////////////
///
/// END FLOAT PRECISION CONSTRUCTORS
///
//////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Implict/Explicit conversions to base types int, short, long, char
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator double
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///    conversion to double operator
///
inline float_precision::operator double() const
   {// Conversion to double
   std::string s = _float_precision_ftoa( this ); 
   return (double)atof( s.c_str() ); 
   } 


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator float
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///    Conversion to float
///
inline float_precision::operator float() const
   {// Conversion to long
   return (float)(double)*this;
   }                      

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator long
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to long
///
inline float_precision::operator long() const
   {// Conversion to long
   return (long)(double)*this;
   }                      

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator int
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to int
///
inline float_precision::operator int() const
   {// Conversion to int 
   return (int)(double)*this;
   }                      


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator short
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to short
///
inline float_precision::operator short() const
   {// Conversion to short
   return (short)(double)*this;
   }                      

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator char
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to char
///
inline float_precision::operator char() const
   {// Conversion to char
   return (char)(double)*this;
   }                      

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator unsigned long
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to unsigned long
///
inline float_precision::operator unsigned long() const
   {// Conversion to int 
   return (unsigned long)(double)*this;
   }                      


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator unsigned int
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to unsigned int
///
inline float_precision::operator unsigned int() const
   {// Conversion to int 
   return (unsigned int)(double)*this;
   }                      


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator unsigned short
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to unsigned short
///
inline float_precision::operator unsigned short() const
   {// Conversion to short
   return (unsigned short)(double)*this;
   }                      

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/12/2006
///	@brief 			float_precision::operator unsigned char
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to unsigned char
///
inline float_precision::operator unsigned char() const
   {// Conversion to char
   return (unsigned char)(double)*this;
   }                      


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/30/2007
///	@brief 			float_precision::operator int_precision
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion to int_precision
///
inline float_precision::operator int_precision() const
    {// Conversion to int_precision
    std::string s = _float_precision_ftoainteger( this );
    return (int_precision)( (char *)s.c_str() );
    } 

//////////////////////////////////////////////////////////////////////////////////////
///
/// FLOAT PRECISION OPERATORS
///
//////////////////////////////////////////////////////////////////////////////////////


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	Assign float precision numbers
///	@return 	float_precision&	-	
///	@param   "a"	-	float precsion number to assign
///
///	@todo
///
/// Description:
///   Assign operator
///   Round it to precision and mode of the left hand side
///   Only the exponent and mantissa is assigned
///   Mode and precision is not affected by the assignment.
//
inline float_precision& float_precision::operator=( const float_precision& a )
   {
   int sign;

   mExpo = a.mExpo;   
   sign = a.sign();
   mNumber = a.mNumber.substr(1);
   if( _float_precision_rounding( &mNumber, sign, mPrec, mRmode ) != 0 )  // Round back to left hand side precision
      mExpo++;

   mNumber.insert( (std::string::size_type)0, SIGN_STRING( sign ) );

   return *this;
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	+= float precision numbers
///	@return 	the resulting float_precision number
///	@param   "a"	-	float precsion number to assign
///
///	@todo    Still missing code for x += a where add make sense. fx. if a is so small it does 
///            not affect the result within the given precision is should ignored. same is true
///            if x is insignififcant comapre to a the just assign a to x
///
/// Description:
///   The essential += operator
///   1) Align to same exponent
///   2) Align to same precision
///   3) Add Mantissa
///   4) Add carry to exponent
///   4) Normalize 
///   5) Rounding to precission 
///   Early out algorithm. i.e.
///      - x+=0 return x
///      - x+=a wher x is 0 return a
//
inline float_precision& float_precision::operator+=( const float_precision& a )
   {
   int sign, sign1, sign2, wrap;
   int expo_max, digits_max;
   unsigned int precision_max;
   std::string s, s1, s2;
   
   if( a.mNumber.length() == 2 && FDIGIT( a.mNumber[1] ) == 0 )  // Add zero
      return *this;
   if( mNumber.length() == 2 && FDIGIT( mNumber[1] ) == 0 )      // Add a (not zero) to *this (is zero) Same as *this = a;
      return *this = a;

   // extract sign and unsigned portion of number
   sign1 = a.sign();
   s1 = a.mNumber.substr( 1 ); // Extract Mantissa
   sign2 = CHAR_SIGN( mNumber[0] );
   s2 = mNumber.substr( 1 );   // Extract Mantissa
   expo_max = std::max( mExpo, a.mExpo );
   precision_max = std::max( mPrec, a.precision() );

   // Check if add makes sense. Still missing
   
   // Right shift (padd leading zeros) to the smallest number
   if( a.mExpo != expo_max )
      _float_precision_right_shift( &s1, expo_max - a.mExpo );
   if( mExpo != expo_max ) 
      _float_precision_right_shift( &s2, expo_max - mExpo );

   // Round to same precision 
   if( _float_precision_rounding( &s1, sign1, precision_max, a.mode() ) != 0 ) // If carry when rounding up then one right shift
      _float_precision_right_shift( &s1, 1 );
   if( _float_precision_rounding( &s2, sign2, precision_max, mRmode ) != 0 ) // If carry when rounding up then one right shift 
      _float_precision_right_shift( &s2, 1 );

   // Alignment to same number of digits, so add can be perfomed as integer add
   digits_max = std::max( s1.length(), s2.length() );
   if( s1.length() != digits_max )
      _float_precision_left_shift( &s1, digits_max - s1.length() );
   if( s2.length() != digits_max ) 
      _float_precision_left_shift( &s2, digits_max - s2.length() );
   
   // Now s1 and s2 is aligned to the same exponent. The biggest of the two
   if( sign1 == sign2 )
      {
      s = _float_precision_uadd( &s1, &s2 );
      if( s.length() > s1.length() ) // One more digit
         expo_max++;
      sign = sign1;
      }
   else
      {
      int cmp = _float_precision_compare( &s1, &s2 );
      if( cmp > 0 ) // Since we subctract less the wrap indicater need not to be checked
         {
         s = _float_precision_usub( &wrap, &s1, &s2 );
         sign = sign1;
         }
      else
         if( cmp < 0 )
            {
            s = _float_precision_usub( &wrap, &s2, &s1 );
            sign = sign2;
            }
         else
            {  // Result zero
            sign = 1;
            s = "0"; s[0] = FCHARACTER( 0 );
            expo_max = 0;
            }
      }

   expo_max += _float_precision_normalize( &s );            // Normalize the number
   if( _float_precision_rounding( &s, sign, mPrec, mRmode ) != 0 )  // Round back left hand side precision
      expo_max++;

   mNumber = SIGN_STRING( sign ) + s;
   mExpo = expo_max;

   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	-= float precision numbers
///	@return 	the resulting float_precision number
///	@param   "a"	-	float precsion number to assign
///
///	@todo    
///
/// Description:
///   The essential -= operator
///   n = n - a is the same as n = n + (-a). so change sign and use the += operator instead
//
inline float_precision& float_precision::operator-=( const float_precision& a )
   {
   float_precision b;

   b.precision( a.precision() ); 
   b.mode( a.mode() );
   b = a;
   b.change_sign();
   *this += b;

   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	*= float precision numbers
///	@return 	the resulting float_precision number
///	@param   "a"	-	float precsion number to assign
///
///	@todo    
///
/// Description:
///   The essential *= operator
///   1) Multiply mantissa
///   2) Add exponent
///   3) Normalize
///   4) Rounding to precision
//
inline float_precision& float_precision::operator*=( const float_precision& a )
   {
   int expo_res;
   int sign, sign1, sign2;
   std::string s, s1, s2;

   // extract sign and unsigned portion of number
   sign1 = a.sign();
   s1 = a.mNumber.substr( 1 );
   sign2 = CHAR_SIGN( mNumber[0] );
   s2 = mNumber.substr( 1 );

   sign = sign1 * sign2;
   // Check for multiplication of 1 digit and use umul_short().
   if(s1.length()==1 )
      s=_float_precision_umul_short( &s2, FDIGIT(s1[0]));
   else
     if( s2.length()==1)
         s=_float_precision_umul_short( &s1, FDIGIT(s2[0]));
      else
         s = _float_precision_umul_fourier( &s1, &s2 );  
   expo_res = mExpo + a.mExpo;
   if( s.length() - 1 > s1.length() + s2.length() - 2 ) // A carry 
      expo_res++;
   expo_res += _float_precision_normalize( &s );            // Normalize the number
   if( _float_precision_rounding( &s, sign, mPrec, mRmode ) != 0 )  // Round back left hand side precision
      expo_res++;
   
   if( sign == -1 && s.length() == 1 && FDIGIT( s[0] ) == 0 )  // Avoid -0 as result +0 is right
      sign = 1; // Change sign
   if( s.length() == 1 && FDIGIT( s[0] ) == 0 )  // Result 0 clear exponent
      expo_res = 0;

   mExpo = expo_res;
   mNumber = SIGN_STRING( sign ) + s;

   return *this;
   }



///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/27/2006
///	@brief 	/= float precision numbers
///	@return 	the resulting float_precision number
///	@param   "a"	-	float precsion number to assign
///
///	@todo    
///
/// Description:
///   The essential /= operator
///   We do a /= b as a *= (1/b)
/// Bug 
///   1/27/2006 Inverse was always done with the precision of a instead of the Max precision of both this & a
//
inline float_precision& float_precision::operator/=( const float_precision& a )
   {
   if( mNumber.length() == 2 && FDIGIT( mNumber[1] ) == 0 ) // If divisor is zero the result is zero
      return *this;

   float_precision c;

   c.precision( a.precision() );
   if( a.precision() < mPrec )
      c.precision( mPrec );
   c=a;
   *this *= _float_precision_inverse(c); 
   return *this;
   }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Mixed Mode arithmetic Unary and Binary operators +, -
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  10/Aug/2014
///	@brief 			operator+
///	@return 	float_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for float_precision + float_precision
///		Specialization for float_precision to avoid ambigous overload
///
/*
inline float_precision operator+( float_precision& lhs, float_precision& rhs)
	{
	float_precision c(rhs);

	if (lhs.precision() > c.precision())
		c.precision(lhs.precision());

	return c += lhs;
	}
*/

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/19/2006
///	@brief 			operator+
///	@return 	float_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for float_precision + <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline float_precision operator+( float_precision& lhs, const _Ty& rhs )
   {
   float_precision c(rhs);

   if( lhs.precision() > c.precision() )
      c.precision( lhs.precision() );

   return c += lhs;
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  3/19/2006
///   @version 1.0
///	@brief 			operator+
///	@return 	float_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> + float_precision 
///
template <class _Ty> inline float_precision operator+( const _Ty& lhs, const float_precision& rhs )
   {
   float_precision c(lhs);

   if( rhs.precision() > c.precision() )
      c.precision( rhs.precision() );

   return c += rhs;
   }

///   @author Henrik Vestermark (hve@hvks.com)
///   @date  3/19/2006
///   @version 1.0
///	@brief 			operator+
///	@return 	float_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for int_precision + float_precision 
///
inline float_precision operator+( int_precision& lhs, float_precision& rhs )
   {
   float_precision c(lhs);

   if( rhs.precision() > c.precision() )
      c.precision( rhs.precision() );

   return c += rhs;
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  Aug/9/2014
///   @version 1.0
///	@brief 			operator+
///	@return 	float_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for int_precision + float_precision 
///
inline float_precision operator+( float_precision& lhs, int_precision& rhs )
	{
	float_precision c(rhs);

	if (lhs.precision() > c.precision())
		c.precision(lhs.precision());

	return c += lhs;
	}

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	Unary + float precision number
///	@return 	the resulting float_precision number
///	@param   "a"	-	float precsion number 
///
///	@todo    
///
/// Description:
///   Unary add. Do nothing and return a
//
inline float_precision operator+( const float_precision& a )
   {
   // Otherwise do nothing
   return a;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/29/2014
///	@brief 			operator-
///	@return 	float_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Sub operator for float_precision - <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline float_precision operator-( float_precision& lhs, const _Ty& rhs)
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d -= c;
	}

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/29/2014
///	@brief 			operator-
///	@return 	float_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Sub operator for  <any other type> - float_precision 
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline float_precision operator-( const _Ty& lhs, const float_precision& rhs )
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d -= c;
	}

///   @author Henrik Vestermark (hve@hvks.com)
///   @date  7/29/2014
///   @version 1.0
///	@brief 			operator-
///	@return 	float_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for int_precision - float_precision 
///
inline float_precision operator-(int_precision& lhs, float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c -= rhs;
	}


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	- float precision numbers
///	@return 	the resulting float_precision number
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number to subtract
///
///	@todo    
///
/// Description:
///   Binary subtract two float_precision numbers
///   Implenting using the essential -= operator
//
/*
inline float_precision operator-( const float_precision& a, const float_precision& b )
   {
   unsigned int precision;
   float_precision c;

   precision = a.precision();
   if( precision < b.precision() )
      precision = b.precision();

   c.precision( precision );
   c = a;
   c -= b; 
   
   return c;
   }
*/

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	Unary - float precision number
///	@return 	the resulting float_precision number
///	@param   "a"	-	float precsion number 
///
///	@todo    
///
/// Description:
///   Unary hypen Just change sign
//
inline float_precision operator-( const float_precision& a )
   {
   float_precision b;

   b.precision( a.precision() );
   b = a;
   b.change_sign();
   
   return b;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/29/2014
///	@brief 			operator-
///	@return 	float_precision	*	return addition of lhs * rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Mul operator for float_precision * <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline float_precision operator*(float_precision& lhs, const _Ty& rhs)
	{
	float_precision c(rhs);

	if (lhs.precision() > c.precision())
		c.precision(lhs.precision());

	return c *= lhs;
	}

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/29/2014
///	@brief 			operator-
///	@return 	float_precision	*	return addition of lhs * rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Mul operator for float_precision * <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline float_precision operator*( const _Ty&lhs, const float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c *= rhs;
	}

///   @author Henrik Vestermark (hve@hvks.com)
///   @date  7/29/2014
///   @version 1.0
///	@brief 			operator*
///	@return 	float_precision	-	return addition of lhs * rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for int_precision * float_precision 
///
inline float_precision operator*(int_precision& lhs, float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c *= rhs;
	}


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	* float precision numbers
///	@return 	the resulting float_precision number
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number to multiply
///
///	@todo    
///
/// Description:
///   Binary multiplying two float_precision numbers
///   Implenting using the essential *= operator
//
/*
inline float_precision operator*( const float_precision& a, const float_precision& b )
   {
   unsigned int precision;
   float_precision c;

   precision = a.precision();
   if( precision < b.precision() )
      precision = b.precision();

   c.precision( precision );

   c = a;
   c *= b; 
   
   return c;
   }
*/

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/29/2014
///	@brief 			operator/
///	@return 	float_precision	/	return addition of lhs / rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Div operator for float_precision / <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline float_precision operator/(float_precision& lhs, const _Ty& rhs)
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d /= c;
	}

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  7/29/2014
///	@brief 			operator/
///	@return 	float_precision	-	return addition of lhs / rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Div operator for  <any other type> / float_precision 
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline float_precision operator/(const _Ty& lhs, const float_precision& rhs)
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d /= c;
	}

///   @author Henrik Vestermark (hve@hvks.com)
///   @date  7/29/2014
///   @version 1.0
///	@brief 			operator/
///	@return 	float_precision	-	return addition of lhs / rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for int_precision / float_precision 
///
inline float_precision operator/(int_precision& lhs, float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c /= rhs;
	}


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	/ float precision numbers
///	@return 	the resulting float_precision number
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number to divide
///
///	@todo    
///
/// Description:
///   Binary divide two float_precision numbers
///   Implenting using the essential /= operator
//
/*
inline float_precision operator/( const float_precision& a, const float_precision& b )
   {
   unsigned int precision;
   float_precision c;

   precision = a.precision();
   if( precision < b.precision() )
      precision = b.precision();

   c.precision( precision + 1 );

   c = a;
   c /= b;
   
   return c;
   }
*/


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	== float precision numberss
///	@return 	the boolean result
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number 
///
///	@todo    
///
/// Description:
///   If both operands has the same mantissa length and same exponent
///   and if the mantissa is identical then it's the same. 
///   However a special test of +0 == -0 is done
///   Precsion and rounding mode does not affect the comparison
//
inline bool operator==( const float_precision& a, const float_precision& b )
   {
   if( const_cast<float_precision&>(a).ref_mantissa()->length() == 2 && (*const_cast<float_precision&>(a).ref_mantissa())[1] == FDIGIT(0) && (*const_cast<float_precision&>(b).ref_mantissa())[1] == FDIGIT(0) )
            return true;  // This conditions is only true if +-0 is compare with =-0 and therefore true. Since the mantissa is zero we dont' need to check the exponent
   if( const_cast<float_precision&>(a).ref_mantissa()->length() != const_cast<float_precision&>(b).ref_mantissa()->length() ||
       a.exponent() != b.exponent() ) // Different therefore false
      return false;
   else
      if( (const_cast<float_precision&>(a).ref_mantissa())->compare( *const_cast<float_precision&>(b).ref_mantissa() ) == 0 )   // Same return true
         return true;
   
   return false;
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	compare < float precision numbers
///	@return 	the boolean result
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number 
///
///	@todo    
///
/// Description:
///   1) Test for both operand is zero and return false if condition is meet
///   2) If signs differs then return the boolean result based on that
///   3) Now if same sign and one operand is zero then return the boolean result
///   4) If same sign and not zero check the exponent
///   5) If same sign and same exponent then check the mantissa for boolean result
///   Precsion and rounding mode does not affect the comparison
//
inline bool operator<( const float_precision& a, const float_precision& b )
   {
   int sign1, sign2, cmp;
   bool zero1, zero2;
   
   sign1 = a.sign(); 
   sign2 = b.sign(); 

   zero1 = const_cast<float_precision&>(a).ref_mantissa()->length() == 2 && FDIGIT( ( *const_cast<float_precision&>(a).ref_mantissa())[1] ) == 0 ? true : false;
   zero2 = const_cast<float_precision&>(b).ref_mantissa()->length() == 2 && FDIGIT( ( *const_cast<float_precision&>(b).ref_mantissa())[1] ) == 0 ? true : false;

   if( zero1 == true && zero2 == true )  // Both zero
      return false;

   // Different signs
   if( sign1 < sign2 )
      return true;
   if( sign1 > sign2 )
      return false;

   // Now a &  b has the same sign
   if( zero1 == true )   // If a is zero and a & b has the same sign and b is not zero then a < b
      return true;
   if( zero2 == true )   // If b is zero and a & b has the same sign and a is not zero then a > b
      return false;

   // Same sign and not zero . Check exponent
   if( a.exponent() < b.exponent() )
      return sign1 > 0 ? true : false;
   if( a.exponent() > b.exponent() )
      return sign1 > 0 ? false: true;

   // Same sign & same exponent. Check mantissa
   cmp = (const_cast<float_precision&>(a).ref_mantissa())->compare( *const_cast<float_precision&>(b).ref_mantissa() );
   if( cmp < 0 && sign1 == 1 )
      return true;
   else
      if( cmp > 0 && sign1 == -1 )
         return true;
   
   return false;
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	!= float precision numberss
///	@return 	the boolean result
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number 
///
///	@todo    
///
/// Description:
///   implemented negating the == comparison
//
inline bool operator!=( const float_precision& a, const float_precision& b )
   {
   return b == a ? false : true;
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	> float precision numberss
///	@return 	the boolean result
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number 
///
///	@todo    
///
/// Description:
///   Implemented using the equality a>b => b<a
//
inline bool operator>( const float_precision& a, const float_precision& b )
   {
   return b < a ? true : false;
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	<= float precision numberss
///	@return 	the boolean result
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number 
///
///	@todo    
///
/// Description:
///   Implemented using the equality a<=b => not b<a
//
inline bool operator<=( const float_precision& a, const float_precision& b )
   {
   return b < a ? false : true;
   }


// Boolean >=
//
///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/21/2005
///	@brief 	>= float precision numberss
///	@return 	the boolean result
///	@param   "a"	-	first float precsion number 
///	@param   "b"	-	second float precsion number 
///
///	@todo    
///
/// Description:
///   Implemented using the equality a>=b => not a<b
//
inline bool operator>=( const float_precision& a, const float_precision& b )
   {
   return a < b ? false: true;
   }


//////////////////////////////////////////////////////////////////////////////////////
///
/// END FLOAT PRECISION OPERATORS
///
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
///
/// BEGIN FLOAT PRECISION FUNCTIONS
///
//////////////////////////////////////////////////////////////////////////////////////

// absolute()
//
///	@author Henrik Vestermark (hve@hvks.com)
///	@date  6/25/2014
///	@brief 	abs(a) float precision numbers
///	@return 	the absoluet value
///	@param   "a"	-	first float precsion number 
///
///	@todo    
///
/// Description:
///   return the absolute value of the float_precision number
//
inline float_precision fabs( const float_precision& a )
{
	return abs(a);
}

//////////////////////////////////////////////////////////////////////////////////////
///
/// END FLOAT PRECISION FUNCTIONS
///
//////////////////////////////////////////////////////////////////////////////////////
#endif