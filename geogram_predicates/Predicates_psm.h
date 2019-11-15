#define GEOGRAM_PSM
#ifndef GEO_STATIC_LIBS
#define GEO_DYNAMIC_LIBS
#endif
/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */


/*
 *  This file is a PSM (pluggable software module)
 *   generated from the distribution of Geogram.
 *
 *  See Geogram documentation on:
 *   http://alice.loria.fr/software/geogram/doc/html/index.html
 *
 *  See documentation of the functions bundled in this PSM on:
 *   http://alice.loria.fr/software/geogram/doc/html/namespaceGEO_1_1PCK.html
 */



/******* extracted from ../api/defs.h *******/

#ifndef GEOGRAM_API_DEFS
#define GEOGRAM_API_DEFS


/*
 * Deactivate warnings about documentation
 * We do that, because CLANG's doxygen parser does not know
 * some doxygen commands that we use (retval, copydoc) and
 * generates many warnings for them...
 */
#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif


#if defined(GEO_DYNAMIC_LIBS)
   #if defined(_MSC_VER)
      #define GEO_IMPORT __declspec(dllimport)
      #define GEO_EXPORT __declspec(dllexport)
   #elif defined(__GNUC__)
      #define GEO_IMPORT
      #define GEO_EXPORT __attribute__ ((visibility("default")))
   #else
      #define GEO_IMPORT
      #define GEO_EXPORT
   #endif
#else
   #define GEO_IMPORT
   #define GEO_EXPORT
#endif

#ifdef geogram_EXPORTS
#define GEOGRAM_API GEO_EXPORT
#else
#define GEOGRAM_API GEO_IMPORT
#endif


#define NO_GEOGRAM_API

typedef int GeoMesh;

typedef unsigned char geo_coord_index_t;

/*
 * If GARGANTUA is defined, then geogram is compiled
 * with 64 bit indices.
 */
#ifdef GARGANTUA

#include <stdint.h>

typedef uint64_t geo_index_t;

typedef int64_t geo_signed_index_t;

#else

typedef unsigned int geo_index_t;

typedef int geo_signed_index_t;

#endif

typedef double geo_coord_t;

typedef int geo_boolean;

enum {
    GEO_FALSE = 0,
    GEO_TRUE = 1
};

#endif


/******* extracted from ../basic/common.h *******/

#ifndef GEOGRAM_BASIC_COMMON
#define GEOGRAM_BASIC_COMMON


// iostream should be included before anything else,
// otherwise 'cin', 'cout' and 'cerr' will be uninitialized.
#include <iostream>



namespace GEO {

    enum {
	GEOGRAM_NO_HANDLER = 0,
	GEOGRAM_INSTALL_HANDLERS = 1
    };

    void GEOGRAM_API initialize(int flags = GEOGRAM_INSTALL_HANDLERS);

    void GEOGRAM_API terminate();
}


#if (defined(NDEBUG) || defined(GEOGRAM_PSM)) && !defined(GEOGRAM_PSM_DEBUG)
#undef GEO_DEBUG
#undef GEO_PARANOID
#else
#define GEO_DEBUG
#define GEO_PARANOID
#endif

// =============================== LINUX defines ===========================

#if defined(__ANDROID__)
#define GEO_OS_ANDROID
#endif

#if defined(__linux__)

#define GEO_OS_LINUX
#define GEO_OS_UNIX

#ifndef GEO_OS_ANDROID
#define GEO_OS_X11
#endif

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(__INTEL_COMPILER)
#  define GEO_COMPILER_INTEL
#elif defined(__clang__)
#  define GEO_COMPILER_CLANG
#elif defined(__GNUC__)
#  define GEO_COMPILER_GCC
#else
#  error "Unsupported compiler"
#endif

// The following works on GCC and ICC
#if defined(__x86_64)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== WINDOWS defines =========================

#elif defined(_WIN32) || defined(_WIN64)

#define GEO_OS_WINDOWS

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(_MSC_VER)
#  define GEO_COMPILER_MSVC
#elif defined(__MINGW32__) || defined(__MINGW64__)
#  define GEO_COMPILER_MINGW
#endif

#if defined(_WIN64)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== APPLE defines ===========================

#elif defined(__APPLE__)

#define GEO_OS_APPLE
#define GEO_OS_UNIX

#if defined(_OPENMP)
#  define GEO_OPENMP
#endif

#if defined(__clang__)
#  define GEO_COMPILER_CLANG
#elif defined(__GNUC__)
#  define GEO_COMPILER_GCC
#else
#  error "Unsupported compiler"
#endif

#if defined(__x86_64) || defined(__ppc64__)
#  define GEO_ARCH_64
#else
#  define GEO_ARCH_32
#endif

// =============================== Emscripten defines  ======================

#elif defined(__EMSCRIPTEN__)

#define GEO_OS_UNIX
#define GEO_OS_LINUX
#define GEO_OS_EMSCRIPTEN
#define GEO_ARCH_64
#define GEO_COMPILER_EMSCRIPTEN

// =============================== Unsupported =============================
#else
#error "Unsupported operating system"
#endif

#if defined(GEO_COMPILER_GCC)   || \
    defined(GEO_COMPILER_CLANG) || \
    defined(GEO_COMPILER_MINGW) || \
    defined(GEO_COMPILER_EMSCRIPTEN)
#define GEO_COMPILER_GCC_FAMILY
#endif

#ifdef DOXYGEN_ONLY
// Keep doxygen happy
#define GEO_OS_WINDOWS
#define GEO_OS_APPLE
#define GEO_OS_ANDROID
#define GEO_ARCH_32
#define GEO_COMPILER_INTEL
#define GEO_COMPILER_MSVC
#endif

#define CPP_CONCAT_(A, B) A ## B

#define CPP_CONCAT(A, B) CPP_CONCAT_(A, B)

#if defined(GOMGEN)
#define GEO_NORETURN
#elif defined(GEO_COMPILER_GCC_FAMILY) || \
      defined(GEO_COMPILER_INTEL)
#define GEO_NORETURN __attribute__((noreturn))
#else
#define GEO_NORETURN
#endif

#if defined(GOMGEN)
#define GEO_NORETURN_DECL
#elif defined(GEO_COMPILER_MSVC)
#define GEO_NORETURN_DECL __declspec(noreturn)
#else
#define GEO_NORETURN_DECL
#endif

#if defined(GEO_COMPILER_CLANG) || defined(GEO_COMPILER_EMSCRIPTEN)
#if __has_feature(cxx_noexcept)
#define GEO_NOEXCEPT noexcept
#endif
#endif

#ifndef GEO_NOEXCEPT
#define GEO_NOEXCEPT throw()
#endif

#define FOR(I,UPPERBND) for(index_t I = 0; I<index_t(UPPERBND); ++I)

#endif


/******* extracted from ../basic/argused.h *******/

#ifndef GEOGRAM_BASIC_ARGUSED
#define GEOGRAM_BASIC_ARGUSED



namespace GEO {

    template <class T>
    inline void geo_argused(const T&) {
    }
}

#endif


/******* extracted from ../basic/numeric.h *******/

#ifndef GEOGRAM_BASIC_NUMERIC
#define GEOGRAM_BASIC_NUMERIC

#include <cmath>
#include <float.h>
#include <limits.h>
#include <algorithm> // for std::min / std::max

// Visual C++ ver. < 2010 does not have C99 stdint.h,
// using a fallback portable one.
#if defined(GEO_OS_WINDOWS) && (_MSC_VER < 1600)
#else
#include <stdint.h>
#endif

#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace GEO {

    namespace Numeric {


        typedef void* pointer;


        typedef int8_t int8;


        typedef int16_t int16;


        typedef int32_t int32;


        typedef int64_t int64;


        typedef uint8_t uint8;


        typedef uint16_t uint16;


        typedef uint32_t uint32;


        typedef uint64_t uint64;


        typedef float float32;


        typedef double float64;

        inline float32 max_float32() {
            return std::numeric_limits<float32>::max();
        }

        inline float32 min_float32() {
            // Note: numeric_limits<>::min() is not
            // what we want (it returns the smallest
            // positive non-denormal).
            return -max_float32();
        }

        inline float64 max_float64() {
            return std::numeric_limits<float64>::max();
        }

        inline float64 min_float64() {
            // Note: numeric_limits<>::min() is not
            // what we want (it returns the smallest
            // positive non-denormal).
            return -max_float64();
        }

        bool GEOGRAM_API is_nan(float32 x);

        bool GEOGRAM_API is_nan(float64 x);

        void GEOGRAM_API random_reset();

        int32 GEOGRAM_API random_int32();

        float32 GEOGRAM_API random_float32();

        float64 GEOGRAM_API random_float64();

        template <class T, bool is_numeric>
        struct LimitsHelper : std::numeric_limits<T> {
        };

        template <class T>
        struct LimitsHelper<T, true> : std::numeric_limits<T> {

            static const size_t size = sizeof(T);

            static const size_t numbits = 8 * sizeof(T);
        };

        template <class T>
        struct Limits :
            LimitsHelper<T, std::numeric_limits<T>::is_specialized> {
        };
    }



    enum Sign {

        NEGATIVE = -1,

        ZERO = 0,

        POSITIVE = 1
    };

    template <class T>
    inline Sign geo_sgn(const T& x) {
        return (x > 0) ? POSITIVE : (
            (x < 0) ? NEGATIVE : ZERO
        );
    }

    template <class T>
    inline T geo_sqr(T x) {
        return x * x;
    }

    template <class T>
    inline void geo_clamp(T& x, T min, T max) {
        if(x < min) {
            x = min;
        } else if(x > max) {
            x = max;
        }
    }

    typedef geo_index_t index_t;

    inline index_t max_index_t() {
        return std::numeric_limits<index_t>::max();
    }

    typedef geo_signed_index_t signed_index_t;

    inline signed_index_t max_signed_index_t() {
        return std::numeric_limits<signed_index_t>::max();
    }

    inline signed_index_t min_signed_index_t() {
        return std::numeric_limits<signed_index_t>::min();
    }

    typedef geo_coord_index_t coord_index_t;

    inline double round(double x) {
	return ((x - floor(x)) > 0.5 ? ceil(x) : floor(x));
    }
}

#endif


/******* extracted from ../basic/psm.h *******/

#ifndef GEOGRAM_BASIC_PSM
#define GEOGRAM_BASIC_PSM


#include <assert.h>
#include <iostream>
#include <string>

#ifndef GEOGRAM_PSM
#define GEOGRAM_PSM
#endif

#ifndef GEOGRAM_BASIC_ASSERT

#define geo_assert(x) assert(x)
#define geo_range_assert(x, min_val, max_val) \
    assert((x) >= (min_val) && (x) <= (max_val))
#define geo_assert_not_reached assert(0)

#ifdef GEO_DEBUG
#define geo_debug_assert(x) assert(x)
#define geo_debug_range_assert(x, min_val, max_val) \
    assert((x) >= (min_val) && (x) <= (max_val))
#else
#define geo_debug_assert(x)
#define geo_debug_range_assert(x, min_val, max_val)
#endif

#ifdef GEO_PARANOID
#define geo_parano_assert(x) geo_assert(x)
#define geo_parano_range_assert(x, min_val, max_val) \
    geo_range_assert(x, min_val, max_val)
#else
#define geo_parano_assert(x)
#define geo_parano_range_assert(x, min_val, max_val)
#endif

#endif

#ifndef geo_cite
#define geo_cite(x)
#endif

#ifndef geo_cite_with_info
#define geo_cite_with_info(x,y)
#endif

#ifndef GEOGRAM_BASIC_LOGGER

namespace GEO {
    namespace Logger {
        inline std::ostream& out(const std::string& name) {
            return std::cout << " [" << name << "]";
        }

        inline std::ostream& err(const std::string& name) {
            return std::cerr << "E[" << name << "]";
        }

        inline std::ostream& warn(const std::string& name) {
            return std::cerr << "W[" << name << "]";
        }
    }

}

#endif

#ifndef FPG_UNCERTAIN_VALUE
#define FPG_UNCERTAIN_VALUE 0
#endif



#ifndef GEOGRAM_BASIC_THREAD_SYNC
#define GEOGRAM_SPINLOCK_INIT 0

namespace GEO {
    namespace Process {

        typedef int spinlock;

        inline void acquire_spinlock(spinlock& x) {
            // Not implemented yet for PSMs
            geo_argused(x);
            geo_assert_not_reached;
        }

        inline void release_spinlock(spinlock& x) {
            // Not implemented yet for PSMs
            geo_argused(x);
            geo_assert_not_reached;
        }
    }
}
#endif

#define GEOGRAM_WITH_PDEL

#endif

/******* extracted from predicates.h *******/

#ifndef GEOGRAM_NUMERICS_PREDICATES
#define GEOGRAM_NUMERICS_PREDICATES



namespace GEO {

    namespace PCK {

	enum SOSMode { SOS_ADDRESS, SOS_LEXICO };

	void GEOGRAM_API set_SOS_mode(SOSMode m);

	SOSMode GEOGRAM_API get_SOS_mode();

        Sign GEOGRAM_API side1_SOS(
            const double* p0, const double* p1,
            const double* q0,
            coord_index_t DIM
        );

        Sign GEOGRAM_API side2_SOS(
            const double* p0, const double* p1, const double* p2,
            const double* q0, const double* q1,
            coord_index_t DIM
        );

        Sign GEOGRAM_API side3_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3,
            const double* q0, const double* q1, const double* q2,
            coord_index_t DIM
        );

        Sign GEOGRAM_API side3_3dlifted_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3,
            double h0, double h1, double h2, double h3,
            const double* q0, const double* q1, const double* q2,
	    bool SOS=true
        );

        Sign GEOGRAM_API side4_SOS(
            const double* p0,
            const double* p1, const double* p2,
            const double* p3, const double* p4,
            const double* q0, const double* q1,
            const double* q2, const double* q3,
            coord_index_t DIM
        );


        Sign GEOGRAM_API side4_3d(
            const double* p0,
            const double* p1, const double* p2,
            const double* p3, const double* p4
        );

        Sign GEOGRAM_API side4_3d_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4
        );

         Sign GEOGRAM_API in_sphere_3d_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3,
            const double* p4
         );


         Sign GEOGRAM_API in_circle_2d_SOS(
             const double* p0, const double* p1, const double* p2,
             const double* p3
         );


         Sign GEOGRAM_API in_circle_3d_SOS(
             const double* p0, const double* p1, const double* p2,
             const double* p3
         );


        Sign GEOGRAM_API in_circle_3dlifted_SOS(
            const double* p0, const double* p1, const double* p2,
            const double* p3,
            double h0, double h1, double h2, double h3,
	    bool SOS=true
        );

        Sign GEOGRAM_API orient_2d(
            const double* p0, const double* p1, const double* p2
        );


#ifndef GEOGRAM_PSM
        inline Sign orient_2d(
            const vec2& p0, const vec2& p1, const vec2& p2
        ) {
            return orient_2d(p0.data(),p1.data(),p2.data());
        }
#endif


        Sign GEOGRAM_API orient_2dlifted_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3,
            double h0, double h1, double h2, double h3
        );


        Sign GEOGRAM_API orient_3d(
            const double* p0, const double* p1,
            const double* p2, const double* p3
        );


#ifndef GEOGRAM_PSM
        inline Sign orient_3d(
            const vec3& p0, const vec3& p1,
            const vec3& p2, const vec3& p3
        ) {
            return orient_3d(p0.data(),p1.data(),p2.data(),p3.data());
        }
#endif

        Sign GEOGRAM_API orient_3dlifted(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4,
            double h0, double h1, double h2, double h3, double h4
        );


        Sign GEOGRAM_API orient_3dlifted_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4,
            double h0, double h1, double h2, double h3, double h4
        );


	Sign GEOGRAM_API det_3d(
	    const double* p0, const double* p1, const double* p2
	);

	Sign GEOGRAM_API det_4d(
	    const double* p0, const double* p1,
	    const double* p2, const double* p3
	);

	Sign GEOGRAM_API det_compare_4d(
	    const double* p0, const double* p1,
	    const double* p2, const double* p3,
	    const double* p4
	);

	bool GEOGRAM_API aligned_3d(
	    const double* p0, const double* p1, const double* p2
	);

	Sign GEOGRAM_API dot_3d(
	    const double* p0, const double* p1, const double* p2
	);

	Sign GEOGRAM_API dot_compare_3d(
	    const double* v0, const double* v1, const double* v2
	);

	bool points_are_identical_2d(
	    const double* p1,
	    const double* p2
	);

	bool GEOGRAM_API points_are_identical_3d(
	    const double* p1,
	    const double* p2
	);

	bool GEOGRAM_API points_are_colinear_3d(
	    const double* p1,
	    const double* p2,
	    const double* p3
        );

	inline Sign orient_3d_inexact(
	    const double* p0, const double* p1,
	    const double* p2, const double* p3
	);

        void GEOGRAM_API show_stats();

        void GEOGRAM_API initialize();

        void GEOGRAM_API terminate();
    }
}

#endif

