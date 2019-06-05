#ifndef INC_PRECISION
#define INC_PRECISION

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
 * Module name     :   iprecision.h
 * Module ID Nbr   :   
 * Description     :   Arbitrary integer precision class
 * --------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version	Author/Date		Description of changes
 * -------  -----------		----------------------
 * 01.01	HVE/020209		Initial release
 * 01.01    HVE/030421     Use FFT for Fast multiplications. Clean up the code
 *                         Now it can be included in different compilation units
 *                         Without causing linker errors
 * 01.02    HVE/030514     % operator was speed optimized to the same speed as the
 *                         the / operator. Furthere more now both /, % operator
 *                         advantages of short_division or short_remainding 
 *                         speeding up the / and % operator
 * 01.03    HVE/030602     A bug for RADIX==BASE_256 in the udiv() and urem()
 *                         function. Quotient was initialized wrongly to "1" instead
 *                         of ICHARACTER(1)
 * 01.04    HVE/030722     A bug in ++ and -- operator did not increment/decrement
 *                         the variable
 * 01.05    HVE/050114     Subtraction of two eaul size number result in -0. this
 *                         has been corrected to +0.
 *                         BASE [2..9] is now also supported. All Base number <= 10
 *                         will be stored as character numbers.
 *                         Binary string can be read by using the notation "0bxxxx"
 *                         where x is either 0 or 1
 * 01.06    HVE/011905     Doxygen comments added
 * 01.07    HVE/050311     Replacing #define with inline functions *
 * 01.08    HVE/050314     Postfix ++ and -- implemented
 * 01.09    HVE/060203     Added const declaration for constant function parameters
 * 01.10    HVE/060217     Added implict conversion between int_precision and standard
 *                         C types double, float, long, int, short and char and the
 *                         corresponding unsigned versions
 * 01.11    HVE/060221     Adding templates for mixed mode arithmetic
 * 01.12    HVE/060307     itostring firs paramater change from unsigned to int to correspond
 *                         to MicroSoft _itoa().
 * 01.13	HVE/20100527   Port to Visual Express studio 10
 * 01.14	HVE/20100609   (unsigned long) operator now used strtoul() instead of wrongly atoi()
 * 01.15	HVE/20100701   Shifting a zero BASE 2 representation left was done wrongly.	
 * 01.16	HVE/20100701   Initializing a precision string to internal BASE 2 resulted in a 
 *							Exception when digit was higher than BASE_2
 * 01.17	HVE/20100702	Added short cuts for shifting left and right with zero or do left/rigth zero shifts
 * 01.18	HVE/20100703	Fix a bug when an addition resultet in zero and RADIX > 10
 * 01.19	HVE/20100703	Fix a bug in >>= & <<= shifting with RADIX > 10
 * 01.20	HVE/20120815	Improve / and % operator for very large numbers
 * 01.21	HVE/SEP 02-2012	Used floating point arithmetic for large integers for / and % operator to significantly speed up operations
 * 01.22	HVE/JUN-12-2012	Use umul_short() for single digit multiplication instead of umult_fourier() in *= operator
 * 01.23	HVE/MAR-07-2014 Moved <limits.h> from fprecision to iprecision.h to ensure it is included even when only working with arbitrary integer
 *							precision and consequencly only use iprecision.h
 * 01.24	HVE/MAY-03-2014	Made the header Visual studio C++ 2013 compatible (avoiding error C2535)
 * 01.25	HVE/JUN-22-2014	Remove the requirements for microsoft stdafx.h precompiled header
 * 01.26	HVE/JUN-28-2014	Added abs(int_precision()
 * 01.27	HVE/NOV-20-2014 improve the % and / operator when rhs is only a single digit
 * 01.28	HVE/AUG-19-2016	Added initialization of iprecision through a std::string
 * 01.29	HVE/NOV-07-2016	Added support for 64bit signerd/unsigned initialization of a int_precision number
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

/* define version string */
static char _VI_[] = "@(#)iprecision.h 01.29 -- Copyright (C) Future Team Aps";

// If _INT_PRECESION_FAST_DIV_REM is defined it will use a magnitude faster div and rem integer operation.
#define _INT_PRECISSION_FAST_DIV_REM

#include <limits.h>
#include <string>
#include <complex>   // Need <complex> to support FFT functions for fast multiplications

// For ANSI please remove comments from the next 3 line
#include <iostream>
#include <ostream>
#include <istream>
#include <cstdlib>
using std::atoi; using std::strtoul;
// End ANSI addition 

#include <stdlib.h> 
#include <string.h>

// RADIX can either be 2, 8, 10, 16 or 256!
static const int BASE_2	  = 2;
static const int BASE_8   = 8;
static const int BASE_10  = 10;
static const int BASE_16  = 16;
static const int BASE_256 = 256;

///
/// @class precision_ctrl
/// @author Henrik Vestermark (hve@hvks.com)
/// @date  2/15/2006
/// @version 1.0
/// @brief  Arbitrary precision control class
///
/// @todo
///
///// Precision control class
///   This keep track of the internal Base for storing int_precision and Float_precision numbers.
///   Default int_precision radix is BASE_10
///   Default float_precision radix is BASE_10
//
class precision_ctrl {
   int mIRadix;			// Internal base of int_precision
   int mFRadix;			// Internal base of float_precision

   public:
      // Constructor
      precision_ctrl( unsigned int ir=BASE_10, unsigned int fr=BASE_10): mIRadix(ir), mFRadix(fr) {}
      
      // Coordinate functions
      inline int I_RADIX() const						{ return mIRadix; }
      inline int I_RADIX( unsigned int ir )		{ return( mIRadix = ir ); }
	  inline int F_RADIX() const						{ return mFRadix; }
	  inline int F_RADIX( unsigned int fr )		{ return( mFRadix = fr ); }
      };

extern precision_ctrl precision_ctrl;

static const int RADIX = BASE_10;			// Set internal base for the arbitrary precision

inline std::string SIGN_STRING( int x )   { return x >=0 ? "+" : "-" ; }
inline int CHAR_SIGN( char x )            { return x == '-' ? -1 : 1; }
inline unsigned char IDIGIT( char x )     { return RADIX <= 10 ? (unsigned char)( x - '0') : (unsigned char)x; }
inline unsigned char IDIGIT10( char x )   { return (unsigned char)( x - '0'); }
inline unsigned char ICHARACTER( char x ) { return RADIX <= 10 ? (unsigned char)( x + '0') : (unsigned char)x; }
inline unsigned char ICHARACTER10( char x){ return (unsigned char)( x + '0'); }
inline int ICARRY( unsigned int x )       { return x / RADIX; }
inline int ISINGLE( unsigned int x )      { return x % RADIX; }


class int_precision;

// Arithmetic
template <class _Ty> inline int_precision operator+( int_precision&, const _Ty& );
template <class _Ty> inline int_precision operator+( const _Ty&, const int_precision& );
inline int_precision operator+( const int_precision& );        // Unary 
inline int_precision operator++( int_precision& );       // Postfix Increment 
inline int_precision operator++( int_precision&, int );  // Postfix Increment 

template <class _Ty> inline int_precision operator-( int_precision&, const _Ty&);
template <class _Ty> inline int_precision operator-( const _Ty&, const int_precision& );
inline int_precision operator-( const int_precision& );        // Unary
inline int_precision operator--( int_precision& );       // Postfix Decrement
inline int_precision operator--( int_precision&, int );  // Postfix Decrement

template <class _Ty> inline int_precision operator*( int_precision&, const _Ty& );
template <class _Ty> inline int_precision operator*( const _Ty&, const int_precision& );
template <class _Ty> inline int_precision operator/( int_precision&, const _Ty& );
template <class _Ty> inline int_precision operator/( const _Ty&, const int_precision& );
template <class _Ty> inline int_precision operator%( int_precision&, const _Ty& );
template <class _Ty> inline int_precision operator%( const _Ty&, const int_precision& );
template <class _Ty> inline int_precision operator<<( int_precision&, const _Ty& );
template <class _Ty> inline int_precision operator<<( const _Ty&, const int_precision& );
template <class _Ty> inline int_precision operator>>( int_precision&, const _Ty& );
template <class _Ty> inline int_precision operator>>( const _Ty&, const int_precision& );
template <class _Ty> inline int_precision operator&( int_precision&, const _Ty& );
template <class _Ty> inline int_precision operator&( const _Ty&, const int_precision& );

// Boolean Comparision Operators
template <class _Ty> inline bool operator==( int_precision&, const _Ty& );
template <class _Ty> inline bool operator==( const _Ty&, const int_precision& );
template <class _Ty> inline bool operator!=( int_precision&, const _Ty& );
template <class _Ty> inline bool operator!=( const _Ty&, const int_precision& );
template <class _Ty> inline bool operator>( int_precision&, const _Ty& );
template <class _Ty> inline bool operator>( const _Ty&, const int_precision& );
template <class _Ty> inline bool operator>=( int_precision&, const _Ty& );
template <class _Ty> inline bool operator>=( const _Ty&, const int_precision& );
template <class _Ty> inline bool operator<=( int_precision&, const _Ty& );
template <class _Ty> inline bool operator<=( const _Ty&, const int_precision& );
template <class _Ty> inline bool operator<( int_precision&, const _Ty& );
template <class _Ty> inline bool operator<( const _Ty&, const int_precision& );

// Integer Precision functions 
extern int_precision abs(const int_precision&);
extern int_precision ipow( const int_precision&, const int_precision& );    // a^b 
extern int_precision ipow_modular( const int_precision&, const int_precision&, const int_precision& ); // a^b%c
extern bool iprime( const int_precision& );

// Core functions that works directly on String class and unsigned arithmetic
void _int_real_fourier( double [], unsigned int, int );
std::string _int_precision_uadd( std::string *, std::string *);
std::string _int_precision_uadd_short( std::string *, unsigned int );
std::string _int_precision_usub( int *, std::string *, std::string *);
std::string _int_precision_usub_short( int *, std::string *, unsigned int );
std::string _int_precision_umul( std::string *, std::string *);
std::string _int_precision_umul_short( std::string *, unsigned int );
std::string _int_precision_umul_fourier( std::string *, std::string *);
std::string _int_precision_udiv( std::string *, std::string *);
std::string _int_precision_udiv_short( unsigned int *, std::string *, unsigned int );
std::string _int_precision_urem( std::string *, std::string *);
std::string _int_precision_uneg( std::string *);
std::string _int_precision_uand( std::string *, std::string *);
int _int_precision_compare( std::string *, std::string * );
void _int_precision_strip_leading_zeros( std::string * );
std::string _int_precision_itoa( const std::string * );
std::string _int_precision_itoa( int_precision * );
std::string _int_precision_atoi( const char *str );
std::string itostring( int, const unsigned );
std::string ito_precision_string( unsigned long, const bool, const int base = RADIX );
std::string i64to_precision_string( uint64_t, const bool, const int base = RADIX);

///
/// @class int_precision
/// @author Henrik Vestermark (hve@hvks.com)
/// @date  9/7/2004
/// @version 1.0
/// @brief  This is an arbitrary integer class
///
/// @todo  	
///
/// Precision class
/// An Arbitrary integer always has the format [sign][digit]+ where sign is either '+' or '-'
/// the length or the representation is always >= 2
/// A null string is considered as an error and an exception is thrown
/// Also number is always strip for leading zeros
/// Since we always initiate to a valid int_precision number, it will never be a empty string
///
class int_precision
	{
   std::string mNumber;

   public:
      // Constructor
	  int_precision()  { mNumber = ICHARACTER(0); }
 	  int_precision( char );				// When initialized through a char
      int_precision( unsigned char );		// When initialized through a unsigned char
      int_precision( short );				// When initialized through an short
      int_precision( unsigned short );		// When initialized through an unsigned short
      int_precision( int );					// When initialized through an int
      int_precision( unsigned int );		// When initialized through an unsigned int
      int_precision( long );				// When initialized through an long
      int_precision( unsigned long );		// When initialized through an unsigned long
      int_precision( const char * );		// When initialized through a char string
	  int_precision( const std::string&  );	// When initialized through a std::string
	  int_precision( const int64_t );		// When initialized through a 64 bit int
	  int_precision( const uint64_t );		// When initialized through a 64 bit unsinged int
	  int_precision( const int_precision& s) : mNumber(s.mNumber) {}  // When initialized through another int_precision

         
      // Coordinate functions
      std::string copy() const   { return mNumber; }
      std::string *pointer()     { return &mNumber; }
      int sign() const           { return CHAR_SIGN( mNumber[0] ); }                    
      int change_sign()          { // Change and return sign   
                                 if( mNumber.length() != 2 || IDIGIT( mNumber[1] ) != 0 ) // Don't change sign for +0!
                                    if( mNumber[0] == '+' ) mNumber[0] = '-'; else mNumber[0] = '+';
                                 return CHAR_SIGN( mNumber[0] );
                                 }
      int size() const           { return mNumber.length(); } // Return number of digits including the sign

	  // Conversion methods. Safer and less ambiguios than overloading implicit/explivit conversion operators
      std::string toString()	{ return _int_precision_itoa(this); }
 
	  // Implict/explicit conversion operators
	  operator long() const;
      operator int() const;
      operator short() const;
      operator char() const;
      operator unsigned long() const;
      operator unsigned int() const;
      operator unsigned short() const;
      operator unsigned char() const;
      operator double() const; 
      operator float() const; 
	  
      // Essential operators
      int_precision& operator=( const int_precision& );
      int_precision& operator+=( const int_precision& );
      int_precision& operator-=( const int_precision& );
      int_precision& operator*=( const int_precision& );
      int_precision& operator/=( const int_precision& );
      int_precision& operator%=( const int_precision& );
      int_precision& operator>>=( const int_precision& );
      int_precision& operator<<=( const int_precision& );
	 // int_precision& operator&=( const int_precision& );
     // int_precision& operator|=( const int_precision& );
	 // int_precision& operator^=( const int_precision& );
	 //
     // Specialization
		friend std::ostream& operator<<( std::ostream& strm, const int_precision& d );
		friend std::istream& operator>>( std::istream& strm, int_precision& d );

      // Exception class
      class bad_int_syntax {}; 
      class out_of_range   {};
      class divide_by_zero {};
   };



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Constructors
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "str"	-	Convert the character string number into a multi precision number
///
///	@todo  
///
/// Description:
///   Constructor 
///   Validate input and convert to internal representation
///   Always add sign if not specified 
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( const char *str )
   {
   std::string s(str);
   
   if( s.empty() )
      { throw bad_int_syntax(); return; }

   mNumber = _int_precision_atoi( str );
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  Aug-19-2016
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "str"	-	Convert the std::string number into a multi precision number
///
///	@todo  
///
/// Description:
///   Constructor 
///   Validate input and convert to internal representation
///   Always add sign if not specified 
///   Only use core base functions to create multi precision numbers
//
//
inline int_precision::int_precision(const std::string& str)
{
	if (str.empty())
		{ throw bad_int_syntax(); return; }

	mNumber = _int_precision_atoi(str.c_str());
}


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "c"	-	the character integer to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initilize with a character
///   Input Always in BASE_10
///   Convert to the internal RADIX 
//
inline int_precision::int_precision( char c )
   {
   if( c < '0' || c > '9' )
      throw bad_int_syntax(); 
   else
      { 
      mNumber = ito_precision_string( IDIGIT10( c ), true );  // Convert to integer
      }
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "c"	-	the character integer to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initilize with a character
///   Input Always in BASE_10
///   Convert to the internal RADIX 
//
inline int_precision::int_precision( unsigned char c )
   {
   if( c < '0' || c > '9' )
      throw bad_int_syntax(); 
   else
      { 
      mNumber = ito_precision_string( IDIGIT10( c ), true );  // Convert to integer
      }
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "i"	-	the binary integer to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initialize with an integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( short i )
   {
   mNumber = ito_precision_string( i, true );
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "i"	-	the binary integer to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initialize with an integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( unsigned short i )
   {
   mNumber = ito_precision_string( i, false );
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "i"	-	the binary integer to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initialize with an integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( int i )
   {
   mNumber = ito_precision_string( i, true );
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "i"	-	the binary unsigned integer to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initialize with an integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( unsigned int i )
   {
   mNumber = ito_precision_string( i, false );
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "i"	-	the binary integer to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initialize with an integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( long i )
   {
   mNumber = ito_precision_string( i, true );
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "i"	-	the binary unsigned long to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initialize with an integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( unsigned long i )
   {
   mNumber = ito_precision_string( i, false );
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  11/7/2016
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "i"	-	the binary int64_t to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initialize with an integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( int64_t i)
	{
	mNumber = i64to_precision_string( (uint64_t)i, true);
	}


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  11/7/2016
///	@brief 	int_precision::int_precision
///	@return 	nothing
///	@param   "i"	-	the binary uint64_t to convert to multi precision number
///
///	@todo  
///
/// Description:
///   Constructor
///   Validate and initialize with an integer
///   Just convert integer to string representation in BASE RADIX
///   The input integer is always BASE_10
///   Only use core base functions to create multi precision numbers
//
inline int_precision::int_precision( uint64_t i)
	{
	mNumber = i64to_precision_string( i, false);
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Implict conversions to base types int, short, long, char
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator long
///	@return 			inline	-	
///
///	@todo  	
///
/// Description:
///   This is the main operator from int_precision to regular long, int, short & char
///   Any explicit or implicit copnversion first convert to standard c long type and then to any other
///   inbuild type int, short, char. As a type long >= int >= short >= char
///   As with regular C type conversion the conversion truncate to desired type and a possible
///   loss in precision is possible
///
inline int_precision::operator long() const
   {// Conversion to long
   long l;
   if( RADIX == BASE_10 ) 
      l = atoi( mNumber.c_str() ); // Do it directly
   else 
      l = atoi( _int_precision_itoa( &mNumber ).c_str() ); // Need to convert from RADIX to BASE_10 )
   return l;
   }    

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator int
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Any explicit or implicit copnversion first convert to standard c long type and then to int
///   As with regular C type conversion the conversion truncate to desired type and a possible
///   loss in precision is possible
///
inline int_precision::operator int() const
   {// Conversion to int
   return (int)(long)*this;
   }    

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator short
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Any explicit or implicit copnversion first convert to standard c long type and then to short
///   As with regular C type conversion the conversion truncate to desired type and a possible
///   loss in precision is possible
///
inline int_precision::operator short() const
   {// Conversion to short
   return (short)((long)*this);
   }    

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator char
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Any explicit or implicit copnversion first convert to standard c long type and then to char
///   As with regular C type conversion the conversion truncate to desired type and a possible
///   loss in precision is possible
///
inline int_precision::operator char() const
   {// Conversion to char
   return (char)((long)*this);
   }    

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator unsigned long
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Any explicit or implicit copnversion first convert to standard c long type and then to int
///   As with regular C type conversion the conversion truncate to desired type and a possible
///   loss in precision is possible
///
inline int_precision::operator unsigned long() const
   {// Conversion to unsigned long
   unsigned long ul;
   if( RADIX == BASE_10 ) 
      ul = strtoul( mNumber.c_str(), NULL, BASE_10 ); // Do it directly
   else 
      ul = strtoul( _int_precision_itoa( &mNumber ).c_str(), NULL, BASE_10 ); // Need to convert from RADIX to BASE_10 )
   return ul;
   }    

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator unsigned int
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Any explicit or implicit copnversion first convert to standard c long type and then to int
///   As with regular C type conversion the conversion truncate to desired type and a possible
///   loss in precision is possible
///
inline int_precision::operator unsigned int() const
   {// Conversion to int
   return (unsigned int)(unsigned long)*this;
   }    

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator unsigned short
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Any explicit or implicit copnversion first convert to standard c long type and then to short
///   As with regular C type conversion the conversion truncate to desired type and a possible
///   loss in precision is possible
///
inline int_precision::operator unsigned short() const
   {// Conversion to short
   return (unsigned short)((unsigned long)*this);
   }    

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator unsigned char
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Any explicit or implicit copnversion first convert to standard c long type and then to char
///   As with regular C type conversion the conversion truncate to desired type and a possible
///   loss in precision is possible
///
inline int_precision::operator unsigned char() const
   {// Conversion to char
   return (unsigned char)((unsigned long)*this);
   }    


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_preceision::operator double
///	@return 			inline	-	
///
///	@todo  
///
/// Description:
///   Conversion from int_precision to double
///
inline int_precision::operator double() const
   {// Conversion to double
   if( RADIX == BASE_10 ) 
      return (double)atof( mNumber.c_str() ); // Do it directly
   else  
      return (double)atof( _int_precision_itoa( &mNumber ).c_str() ); // Need to convert from RADIX to BASE_10
   }                               

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/17/2006
///	@brief 			int_precision::operator float
///	@return 			inline	-	
///
///	@todo  Add to do things	
///
/// Description:
///   Conversion from int_precision to float 
///   Using the double conversion frist and then trunk to float using standard c conversion
///
inline int_precision::operator float() const
   {// Conversion to float
   return (float)((double)*this);
   }    

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Essentialsoperators =, +=, -=, *=, /=, %=, <<=, >>=, &=, |=, ^=
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator=
///	@return 	static int_precision	-	return a=b
///	@param   "a"	-	Assignment operand
///
///	@todo 
///
/// Description:
///   Assign operator
//
inline int_precision& int_precision::operator=( const int_precision& a )
   {
   mNumber = a.mNumber;

   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator+=
///	@return 	static int_precision	-	return a +=b
///	@param   "a"	-	Adding operand
///
///	@todo 
///
/// Description:
///   += operator
//
inline int_precision& int_precision::operator+=( const int_precision& a )
   {
   int sign1, sign2, wrap;
   std::string s1, s2;
	int cmp;

   // extract sign and unsigned portion of number
   sign1 = a.sign(); 
   s1 = a.mNumber.substr( 1 );
   sign2 = CHAR_SIGN( mNumber[0] );
   s2 = mNumber.substr( 1 );

   if( sign1 == sign2 )
      mNumber = SIGN_STRING( sign1 ) + _int_precision_uadd( &s1, &s2 );
   else
      {
		cmp = _int_precision_compare( &s1, &s2 );
		if( cmp > 0 ) // Since we subctract less the wrap indicater need not to be checked
         mNumber = SIGN_STRING(sign1) + _int_precision_usub( &wrap, &s1, &s2 );
      else
			if( cmp < 0 )
				mNumber = SIGN_STRING(sign2) + _int_precision_usub( &wrap, &s2, &s1 );
			else
				{// result is 0
				mNumber = std::string( "+" );
				mNumber.append( 1, ICHARACTER(0) );
				}
		}
   
   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator-=
///	@return 	static int_precision	-	return a -=b
///	@param   "a"	-	Subtracting operand
///
///	@todo 
///
/// Description:
///   -= operator
///   The essential -= operator
///   n = n - a is the same as n = n + (-a); 
//
inline int_precision& int_precision::operator-=( const int_precision& a )
   {
   int_precision b;

   b = a;
   b.change_sign();
   *this += b;

   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator*=
///	@return 	static int_precision	-	return a *=b
///	@param   "a"	-	Multiplying operand
///
///	@todo 
///
/// Description:
///   *= operator
//
inline int_precision& int_precision::operator*=( const int_precision& a )
   {
   int sign1, sign2;
   std::string s, s1, s2;

   // extract sign and unsigned portion of number
   sign1 = a.sign();
   s1 = a.mNumber.substr( 1 );
   sign2 = CHAR_SIGN( mNumber[0] );
   s2 = mNumber.substr( 1 );

   sign1 *= sign2;  // Check for multiplication of 1 digit and use umul_short().
   if(s1.length()==1 )
   	  s=_int_precision_umul_short( &s2, IDIGIT(s1[0]));
   else
	 if( s2.length()==1)
		 s=_int_precision_umul_short( &s1, IDIGIT(s2[0]));
      else
		  s=_int_precision_umul_fourier( &s1, &s2 );

   mNumber = SIGN_STRING( sign1 ) + s;
   if( sign1 == -1 && mNumber.length() == 2 && IDIGIT( mNumber[1] ) == 0 )  // Avoid -0 as result +0 is right
      mNumber[0] = '-';  // Change sign

   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator/=
///	@return 	static int_precision	-	return a /=b
///	@param   "a"	-	Dividing operand
///
///	@todo 
///
/// Description:
///   /= operator
//
inline int_precision& int_precision::operator/=( const int_precision& a )
   {
   int sign1, sign2;
   std::string s1, s2;

#ifdef _INT_PRECISSION_FAST_DIV_REM
   if (this->size()>a.size() + 8 && a.size() != 2)  // Check that lhs is 8 digit larger and that rhs is not a single digit before do the fastremdiv operation
	   {
	   extern int_precision _int_precision_fastdiv( const int_precision&, const int_precision& );
	   int_precision b=*this;
	   *this = _int_precision_fastdiv( b, a );
	   return *this;
       }
#endif

   // extract sign and unsigned portion of number
   sign1 = CHAR_SIGN(mNumber[0] );
   s1 = mNumber.substr( 1 );
   sign2 = a.sign();
   s2 = a.mNumber.substr( 1 );

   sign1 *= sign2;
   mNumber = SIGN_STRING( sign1 ) + _int_precision_udiv( &s1, &s2 );
   if( sign1 == -1 && mNumber.length() == 2 && IDIGIT( mNumber[1] ) == 0 )  // Avoid -0 as result +0 is right
      mNumber[0] = '+';

   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator%=
///	@return 	static int_precision	-	return a %=b
///	@param   "a"	-	Modulus operand
///
///	@todo 
///
/// Description:
///   %= operator
//
inline int_precision& int_precision::operator%=( const int_precision& a )
   {
   int sign1, sign2;
   std::string s1, s2;

#ifdef _INT_PRECISSION_FAST_DIV_REM
   if(this->size()>a.size()+8 && a.size() != 2 )  // Check that lhs is 8 digit larger and that rhs is not a single digit before do the fastremdiv operation
	   {
	   extern int_precision _int_precision_fastrem( const int_precision&, const int_precision& );
	   int_precision b=*this;
	   *this =_int_precision_fastrem( b, a );
	   return *this;
       }
#endif

   // extract sign and unsigned portion of number
   sign1 = CHAR_SIGN(mNumber[0] );
   s1 = mNumber.substr( 1 );
   sign2 = a.sign();
   s2 = a.mNumber.substr( 1 );

   mNumber = SIGN_STRING( sign1 ) + _int_precision_urem( &s1, &s2 );
   if( sign1 == -1 && mNumber.length() == 2 && IDIGIT( mNumber[1] ) == 0 )  // Avoid -0 as result +0 is right
      mNumber[0] = '+';

   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator<<=
///	@return static int_precision	-	return shifting a<<= b
///	@param  "a"	-	Shifting number
///
///	@todo 
///
/// Description:
///   <<= operator
//
inline int_precision& int_precision::operator<<=( const int_precision& a )
   {
   int sign1, sign2, wrap, max_shifts=0;
   std::string s1, s2;
   std::string cshifts, c0, c3, c4, c8;

   c0.insert( (std::string::size_type)0, 1, ICHARACTER( 0 ) );
   // extract sign and unsigned portion of number
   sign1 = CHAR_SIGN(mNumber[0] );
   s1 = mNumber.substr( 1 );
   if( _int_precision_compare( &s1, &c0 ) == 0 )  // Short cut: zero shifting left is still zero.
		return *this;

   sign2 = a.sign();
   s2 = a.mNumber.substr( 1 );
   if( _int_precision_compare( &s2, &c0 ) == 0 )  // Short cut: shift zero left does not change the number.
		return *this;

   if( sign2 < 0 )
      { throw out_of_range(); return *this; }

   // Speed up the operation by shifting the native if possible. Finx max shifting
	for( int i=RADIX; i>=2; max_shifts++, i>>=1 ) ;
	cshifts.insert( (std::string::size_type)0, 1, ICHARACTER( max_shifts ) );
   if( RADIX > BASE_2 )
      {// BASE max Shift bit a a time if possible
      for( ; _int_precision_compare( &s2, &cshifts ) >= 0; s2 = _int_precision_usub_short( &wrap, &s2, max_shifts ) )
         s1 = _int_precision_umul_short( &s1, 1<<max_shifts);  //1<<max_shifts==256 or 2^8
      }
 	else
	   { // BASE 2
		int shift;  // Don't shift. just add zeros number of shift times!
		s2.insert( 0, "+" );
		shift = atoi( _int_precision_itoa( &s2 ).c_str() ); // Need to convert from RADIX to BASE_10 )
		s1.append( shift, ICHARACTER( 0 ) );
		s2 = ICHARACTER(0);
		}

   // Take the remainds of shifts
   if( _int_precision_compare( &s2, &c0 ) > 0 )
      {
      int shift;

      if( RADIX == BASE_10 ) 
         shift = atoi( s2.c_str() ); // Do it directly
      else 
         {
         s2.insert( 0, "+" );
         shift = atoi( _int_precision_itoa( &s2 ).c_str() ); // Need to convert from RADIX to BASE_10
			s2 = s2.substr( 1 );
         }

		if( RADIX >= BASE_10 && RADIX > (2<<(shift-1)) )
			s1 = _int_precision_umul_short( &s1, 2 << ( shift - 1 ) ); 
		else
			{
			for( ; shift > 0; --shift )
				s1 = _int_precision_umul_short( &s1, 2 ); 
			}
      }

   mNumber = SIGN_STRING( sign1 ) + s1;
   
   return *this;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator>>=
///	@return int_precision	-	return shifting a>>= b
///	@param   "a"	-	Shifting operand
///
///	@todo 
///
/// Description:
///   >>= operator
//
inline int_precision& int_precision::operator>>=( const int_precision& a )
   {
   int sign1, sign2, wrap, max_shifts=0;
   unsigned int rem;
   std::string s1, s2;
   std::string cshifts, c0, c3, c4, c8;
   
   c0.insert((std::string::size_type)0, 1, ICHARACTER(0) );
   // extract sign and unsigned portion of number
   sign1 = CHAR_SIGN(mNumber[0] );
   s1 = mNumber.substr( 1 );
   if( _int_precision_compare( &s1, &c0 ) == 0 )  // Short cut: zero shifting right is still zero.
		return *this;
   
   sign2 = a.sign();
   s2 = a.mNumber.substr( 1 );
   if( _int_precision_compare( &s2, &c0 ) == 0 )  // Short cut: shift zero right does not change the number.
		return *this;

   if( sign2 < 0 )
      { throw out_of_range(); return *this; }

   for( int i=RADIX; i>=2; max_shifts++, i>>=1 ) ;
	cshifts.insert( (std::string::size_type)0, 1, ICHARACTER( max_shifts ) );
   if( RADIX > BASE_2 )  // Speed up by alowing shift with 8 (2^8) at a time instead of single shift
      {
     for( ; _int_precision_compare( &s2, &cshifts ) >= 0; s2 = _int_precision_usub_short( &wrap, &s2, max_shifts ) )
         s1 = _int_precision_udiv_short( &rem, &s1, 1 << max_shifts ); //1<<max_shifts 
      }
   else
		{
		int shift;

		s2.insert( 0, "+" );
		shift = atoi( _int_precision_itoa( &s2 ).c_str() ); // Need to convert from RADIX to BASE_10 )
		if( (int)s1.size() <= shift )
			s1 = ICHARACTER( 0 );
		else
			s1 = s1.substr( 0, s1.size() - shift );
		s2 = ICHARACTER(0);
		}

   // Take the remainds of shifts less after applying the speed up trick
   if( _int_precision_compare( &s2, &c0 ) > 0 )
      {
      int shift;

      if( RADIX == BASE_10 ) 
         shift = atoi( s2.c_str() ); // Do it directly
      else 
         {
         s2.insert( 0, "+" );
         shift = atoi( _int_precision_itoa( &s2 ).c_str() ); // Need to convert from RADIX to BASE_10 )
			s2 = s2.substr( 1 );
         }

		if( RADIX >= BASE_10 && RADIX > (2<<(shift-1)) )
			s1 = _int_precision_udiv_short( &rem, &s1, 2 << ( shift - 1 ) );
		else
			for( ; shift > 0; --shift )
				s1 = _int_precision_udiv_short( &rem, &s1, 2 );
      }

   mNumber = SIGN_STRING( sign1 ) + s1;
   if( sign1 == -1 && mNumber.length() == 2 && IDIGIT( mNumber[1] ) == 0 )  // Avoid -0 as result +0 is right
      mNumber[0] = '+';
   
   return *this;
   }
/*
///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2-Sep/2012
///	@brief 	operator&=
///	@return 	static int_precision	-	return a *=b
///	@param   "a"	-	Anding operand
///
///	@todo 
///
/// Description:
///   &= operator
//
inline int_precision& int_precision::operator&=( const int_precision& a )
   {
   int sign1, sign2;
   std::string s1, s2;

   // extract sign and unsigned portion of number
   sign1 = CHAR_SIGN(mNumber[0] );
   s1 = mNumber.substr( 1 );
   sign2 = a.sign();
   s2 = a.mNumber.substr( 1 );
   if( sign1 == '-' && sign2 == '+' ) sign1 = '+';

//   sign1 *= sign2;
   mNumber = SIGN_STRING( sign1 ) + _int_precision_uand( &s1, &s2 );
   if( sign1 == -1 && mNumber.length() == 2 && IDIGIT( mNumber[1] ) == 0 )  // Avoid -0 as result +0 is right
      mNumber[0] = '+';

   return *this;
   }
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Arithmetic
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/19/2006
///	@brief 			operator+
///	@return 	int_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for int_precision + <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline int_precision operator+( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) += rhs; 
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  2/19/2006
///   @version 1.0
///	@brief 			operator+
///	@return 	int_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> + int_precision 
///
template <class _Ty> inline int_precision operator+( const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) += rhs; 
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator unary +
///	@return 	int_precision	-	a
///	@param   "a"	-	operand 
///
///	@todo 
///
/// Description:
///   Unary + operator
///   Do nothing
//
inline int_precision operator+( const int_precision& a )
   {
   return a;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator++ Prefix
///	@return 	int_precision	-	return the incremented a
///	@param   "a"	-	operand
///
///	@todo 
///
/// Description:
///   Increment operator
//
inline int_precision operator++( int_precision& a )
   {
   a += int_precision( 1 );
   return a;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/14/2005
///	@brief 	operator++ Postfix
///	@return 	int_precision	-	return the a before incrementation
///	@param   "a"	-	operand
///
///	@todo 
///
/// Description:
///   Postfix Increment operator
//
inline int_precision operator++( int_precision& a, int )
   {
   int_precision postfix_a(a);

   a += int_precision( 1 );
   return postfix_a;
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/19/2006
///	@brief 			operator-
///	@return 	int_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for int_precision - <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline int_precision operator-( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) -= rhs; 
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  2/19/2006
///   @version 1.0
///	@brief 			operator-
///	@return 	int_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> - int_precision 
///
template <class _Ty> inline int_precision operator-( const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) -= rhs; 
   }



///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator unary - 
///	@return 	int_precision	-	-a
///	@param   "a"	-	operand for sign change
///
///	@todo 
///
/// Description:
///   Unary - operator
///   Change sign
//
inline int_precision operator-( const int_precision& a )
   {
   int_precision b;

   b = a;
   b.change_sign();
   
   return b;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 	operator-- prefix
///	@return 	int_precision	-	return the decremented a
///	@param   "a"	-	operand
///
///	@todo 
///
/// Description:
///   Decrement operator
//
inline int_precision operator--( int_precision& a )
   {
   a -= int_precision( 1 );
   return a;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  3/14/2005
///	@brief 	operator-- postfix
///	@return 	int_precision	-	return the a before decrementation
///	@param   "a"	-	operand
///
///	@todo 
///
/// Description:
///   Postfix Decrement operator
//
int_precision operator--( int_precision& a, int )
   {
   int_precision postfix_a(a);
   a -= int_precision( 1 );
   return postfix_a;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/19/2006
///	@brief 			operator*
///	@return 	int_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for int_precision * <any other type>
///
template <class _Ty> inline int_precision operator*( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) *= rhs; 
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  2/19/2006
///   @version 1.0
///	@brief 			operator*
///	@return 	int_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> * int_precision 
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline int_precision operator*( const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) *= rhs; 
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/19/2006
///	@brief 			operator/
///	@return 	int_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for int_precision / <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline int_precision operator/( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) /= rhs; 
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  2/19/2006
///   @version 1.0
///	@brief 			operator*
///	@return 	int_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> / int_precision 
///
template <class _Ty> inline int_precision operator/( const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) /= rhs; 
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/19/2006
///	@brief 			operator%
///	@return 	int_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for int_precision % <any other type>
///
template <class _Ty> inline int_precision operator%( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) %= rhs; 
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  2/19/2006
///   @version 1.0
///	@brief 			operator%
///	@return 	int_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> % int_precision 
///
template <class _Ty> inline int_precision operator%( const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) %= rhs; 
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/20/2006
///	@brief 			operator<<
///	@return 	int_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for int_precision << <any other type>
///
template <class _Ty> inline int_precision operator<<( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) <<= rhs; 
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  2/20/2006
///   @version 1.0
///	@brief 			operator<<
///	@return 	int_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> << int_precision 
///
template <class _Ty> inline int_precision operator<<(  const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) <<= rhs; 
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/20/2006
///	@brief 			operator>>
///	@return 	int_precision	-	return addition of lhs + rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   Add operator for int_precision >> <any other type>
///
template <class _Ty> inline int_precision operator>>( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) >>= rhs; 
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  2/20/2006
///   @version 1.0
///	@brief 			operator>>
///	@return 	int_precision	-	return addition of lhs - rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> >> int_precision 
///
template <class _Ty> inline int_precision operator>>( const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) >>= rhs; 
   }

/*
///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2-Sep/2012
///	@brief 			operator&
///	@return 	int_precision	-	return addition of lhs & rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
///	@todo  Add to do things	
///
/// Description:
///   And operator for int_precision & <any other type>
///   no const on the lhs parameter to prevent ambigous overload
///
template <class _Ty> inline int_precision operator&( int_precision& lhs, const _Ty& rhs )
   {
   return int_precision(lhs) &= rhs; 
   }


///   @author Henrik Vestermark (hve@hvks.com)
///   @date  2-Sep/2012
///   @version 1.0
///	@brief 			operator&
///	@return 	int_precision	-	return addition of lhs & rhs
///	@param   "lhs"	-	First operand
///	@param   "rhs"	-	Second operand
///
/// Description:
///   Add operator for <any other type> & int_precision 
///
template <class _Ty> inline int_precision operator&( const _Ty& lhs, const int_precision& rhs )
   {
   return int_precision(lhs) &= rhs; 
   }
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    Comparison
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/20/2006
///	@brief 			operator==
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean equal of two precision numbers. Early out algorithm
///   1) if sign different then result is false. We actual don't do this test because of -0==+0 we should not ocuured but just in case
///   2) if length is different then the result is false
///   3) use core compare to determine boolean value
//
template <class _Ty> inline bool operator==( int_precision& a, const _Ty& b )
   {int_precision c(b);
   if( _int_precision_compare( const_cast<int_precision&>(a).pointer(), const_cast<int_precision&>(c).pointer() ) == 0 )   // Same return true
      return true;
   
   return false;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  2/20/2006
///	@brief 			operator==
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean equal of two precision numbers. Early out algorithm
///   1) if sign different then result is false. We actual don't do this test because of -0==+0 we should not ocuured but just in case
///   2) if length is different then the result is false
///   3) use core compare to determine boolean value
//
template <class _Ty> inline bool operator==( const _Ty& a, const int_precision& b )
   {int_precision c(a);
  // if( _int_precision_compare( const_cast<int_precision&>(int_precision(c)).pointer(), const_cast<int_precision&>(b).pointer() ) == 0 )   // Same return true
   if( _int_precision_compare( const_cast<int_precision&>(c).pointer(), const_cast<int_precision&>(b).pointer() ) == 0 )    return true;
    
   return false;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator<
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean less of two precision numbers. Early out algorithm for higher performance
///   1) If sign different determine boolean result based on sign
///   2) Otherwise determine boolean result based length of number amnd the sign
///   3) Same sign and same length. Do a core comparison and return the result
//
template <class _Ty> inline bool operator<( int_precision& a, const _Ty& c )
   {
   int sign1, sign2, cmp;
   int_precision b(c);

   sign1 = a.sign(); 
   sign2 = b.sign(); 

   // Different signs
   if( sign1 > sign2 )
      return false;
   if( sign1 < sign2 )
      return true;

   // Same sign
   if( sign1 == 1 && a.size() < b.size() ) // Different therefore true
      return true;
   if( sign1 == 1 && a.size() > b.size() ) // Different therefore false
      return false;
   if( sign1 == -1 && a.size() > b.size() )
      return true;
   if( sign1 == -1 && a.size() < b.size() )
      return false;

   // Same sign and same length
   cmp = _int_precision_compare( const_cast<int_precision&>(a).pointer(), const_cast<int_precision&>(b).pointer() );
   if( cmp < 0 && sign1 == 1 )
      return true;
   else
      if( cmp > 0 && sign1 == -1 )
         return true;
   
   return false;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator<
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean less of two precision numbers. Early out algorithm for higher performance
///   1) If sign different determine boolean result based on sign
///   2) Otherwise determine boolean result based length of number amnd the sign
///   3) Same sign and same length. Do a core comparison and return the result
//
template <class _Ty> inline bool operator<( const _Ty& c, const int_precision& b )
   {
   int sign1, sign2, cmp;
   int_precision a(c);

   sign1 = a.sign(); 
   sign2 = b.sign(); 

   // Different signs
   if( sign1 > sign2 )
      return false;
   if( sign1 < sign2 )
      return true;

   // Same sign
   if( sign1 == 1 && a.size() < b.size() ) // Different therefore true
      return true;
   if( sign1 == 1 && a.size() > b.size() ) // Different therefore false
      return false;
   if( sign1 == -1 && a.size() > b.size() )
      return true;
   if( sign1 == -1 && a.size() < b.size() )
      return false;

   // Same sign and same length
   cmp = _int_precision_compare( const_cast<int_precision&>(a).pointer(), const_cast<int_precision&>(b).pointer() );
   if( cmp < 0 && sign1 == 1 )
      return true;
   else
      if( cmp > 0 && sign1 == -1 )
         return true;
   
   return false;
   }


///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator!=
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean not equal of two precision numbers
//
template <class _Ty> inline bool operator!=( int_precision& a, const _Ty& b )
   {
   return a == b ? false : true;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator!=
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean not equal of two precision numbers
//
template <class _Ty> inline bool operator!=( const _Ty& a, const int_precision& b )
   {
   return a == b ? false : true;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator>
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean greater of two precision numbers
//
template <class _Ty> inline bool operator>( int_precision& a, const _Ty& b )
   {
   return b < a ? true : false;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator>
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean greater of two precision numbers
//
template <class _Ty> inline bool operator>( const _Ty& a, const int_precision& b )
   { int_precision c(a);
   return b < c ? true : false;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator<=
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean less or equal of two precision numbers
//
template <class _Ty> inline bool operator<=( int_precision& a, const _Ty& b )
   {
   return b < a ? false : true;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator<=
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean less or equal of two precision numbers
//
template <class _Ty> inline bool operator<=( const _Ty& a, const int_precision& b )
   { int_precision c(a);
   return b < c ? false : true;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator>=
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean greater or equal of two precision numbers
//
template <class _Ty> inline bool operator>=( int_precision& a, const _Ty& b )
   {
   return a < b ? false: true;
   }

///	@author Henrik Vestermark (hve@hvks.com)
///	@date  1/19/2005
///	@brief 			operator>=
///	@return 			bool	- 	the boolean value of the operator
///	@param "a"	-	First operand number to compare
///	@param "b"	-	Second operand number to 
///
///	@todo
///
/// Description:
///   Boolean less or equal of two precision numbers
//
template <class _Ty> inline bool operator>=( const _Ty& a, const int_precision& b )
   {
   return a < b ? false: true;
   }
  
#endif