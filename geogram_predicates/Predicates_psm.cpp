#include "Predicates_psm.h"

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



/******* extracted from ../basic/memory.h *******/

#ifndef GEOGRAM_BASIC_MEMORY
#define GEOGRAM_BASIC_MEMORY

#include <vector>
#include <string.h>
#include <stdlib.h>

#ifdef GEO_OS_WINDOWS

#include <windows.h>
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#else

#include <unistd.h>

#endif


namespace GEO {

    namespace Memory {

        typedef unsigned char byte;


        typedef unsigned char word8;


        typedef unsigned short word16;


        typedef unsigned int word32;


        typedef byte* pointer;


	typedef void (*function_pointer)();

        inline void clear(void* addr, size_t size) {
            ::memset(addr, 0, size);
        }

        inline void copy(void* to, const void* from, size_t size) {
            ::memcpy(to, from, size);
        }

	inline pointer function_pointer_to_generic_pointer(function_pointer fptr) {
	    // I know this is ugly, but I did not find a simpler warning-free
	    // way that is portable between all compilers.
	    pointer result = nullptr;
	    ::memcpy(&result, &fptr, sizeof(pointer));
	    return result;
	}

	inline function_pointer generic_pointer_to_function_pointer(pointer ptr) {
	    // I know this is ugly, but I did not find a simpler warning-free
	    // way that is portable between all compilers.
	    function_pointer result = nullptr;
	    ::memcpy(&result, &ptr, sizeof(pointer));
	    return result;
	}

	inline function_pointer generic_pointer_to_function_pointer(void* ptr) {
	    // I know this is ugly, but I did not find a simpler warning-free
	    // way that is portable between all compilers.
	    function_pointer result = nullptr;
	    ::memcpy(&result, &ptr, sizeof(pointer));
	    return result;
	}

#define GEO_MEMORY_ALIGNMENT 64

        template <int DIM>
        struct PointAlignment {
            static const size_t value = 1;
        };

        template <>
        struct PointAlignment<2> {
            static const size_t value = 16;
        };

        template <>
        struct PointAlignment<3> {
            static const size_t value = 8;
        };

        template <>
        struct PointAlignment<4> {
            static const size_t value = 32;
        };

        template <>
        struct PointAlignment<6> {
            static const size_t value = 16;
        };

        template <>
        struct PointAlignment<8> {
            static const size_t value = 64;
        };

#define geo_dim_alignment(dim) GEO::Memory::PointAlignment<dim>::value

        inline void* aligned_malloc(
            size_t size, size_t alignment = GEO_MEMORY_ALIGNMENT
        ) {
#if   defined(GEO_OS_ANDROID)
            // Alignment not supported under Android.
            geo_argused(alignment);
            return malloc(size);
#elif defined(GEO_COMPILER_INTEL)
            return _mm_malloc(size, alignment);
#elif defined(GEO_COMPILER_GCC) || defined(GEO_COMPILER_CLANG)
            void* result;
            return posix_memalign(&result, alignment, size) == 0
                   ? result : nullptr;
#elif defined(GEO_COMPILER_MSVC)
            return _aligned_malloc(size, alignment);
#else
            geo_argused(alignment);
            return malloc(size);
#endif
        }

        inline void aligned_free(void* p) {
#if   defined(GEO_OS_ANDROID)
            // Alignment not supported under Android.
            free(p);
#elif defined(GEO_COMPILER_INTEL)
            _mm_free(p);
#elif defined(GEO_COMPILER_GCC_FAMILY)
            free(p);
#elif defined(GEO_COMPILER_MSVC)
            _aligned_free(p);
#else
            free(p);
#endif
        }

#if   defined(GEO_OS_ANDROID)
#define geo_decl_aligned(var) var
#elif defined(GEO_COMPILER_INTEL)
#define geo_decl_aligned(var) __declspec(aligned(GEO_MEMORY_ALIGNMENT)) var
#elif defined(GEO_COMPILER_GCC_FAMILY)
#define geo_decl_aligned(var) var __attribute__((aligned(GEO_MEMORY_ALIGNMENT)))
#elif defined(GEO_COMPILER_MSVC)
#define geo_decl_aligned(var) __declspec(align(GEO_MEMORY_ALIGNMENT)) var
#elif defined(GEO_COMPILER_EMSCRIPTEN)
#define geo_decl_aligned(var) var
#endif

#if   defined(GEO_OS_ANDROID)
#define geo_assume_aligned(var, alignment)
#elif defined(GEO_COMPILER_INTEL)
#define geo_assume_aligned(var, alignment) \
    __assume_aligned(var, alignment)
#elif defined(GEO_COMPILER_CLANG)
#define geo_assume_aligned(var, alignment)
        // GCC __builtin_assume_aligned is not yet supported by clang-3.3
#elif defined(GEO_COMPILER_GCC)
#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 7
#define geo_assume_aligned(var, alignment) \
        *(void**) (&var) = __builtin_assume_aligned(var, alignment)
        // the GCC way of specifying that a pointer is aligned returns
        // the aligned pointer (I can't figure out why). It needs to be
        // affected otherwise it is not taken into account (verified by
        // looking at the output of gcc -S)
#else
#define geo_assume_aligned(var, alignment)
#endif
#elif defined(GEO_COMPILER_MSVC)
#define geo_assume_aligned(var, alignment)
        // TODO: I do not know how to do that with MSVC
#elif defined(GEO_COMPILER_EMSCRIPTEN)
#define geo_assume_aligned(var, alignment)
#elif defined(GEO_COMPILER_MINGW)
#define geo_assume_aligned(var, alignment)
#endif

#if   defined(GEO_COMPILER_INTEL)
#define geo_restrict __restrict
#elif defined(GEO_COMPILER_GCC_FAMILY)
#define geo_restrict __restrict__
#elif defined(GEO_COMPILER_MSVC)
#define geo_restrict __restrict
#elif defined(GEO_COMPILER_EMSCRIPTEN)
#define geo_restrict
#endif

        inline bool is_aligned(
            void* p, size_t alignment = GEO_MEMORY_ALIGNMENT
        ) {
            return (reinterpret_cast<size_t>(p) & (alignment - 1)) == 0;
        }

        inline void* align(void* p) {
            size_t offset = (
                GEO_MEMORY_ALIGNMENT -
                (reinterpret_cast<size_t>(p) & (GEO_MEMORY_ALIGNMENT - 1))
            ) & (GEO_MEMORY_ALIGNMENT - 1);
            return reinterpret_cast<char*>(p) + offset;
        }

#define geo_aligned_alloca(size) \
    GEO::Memory::align(alloca(size + GEO_MEMORY_ALIGNMENT - 1))

        template <class T, int ALIGN = GEO_MEMORY_ALIGNMENT>
        class aligned_allocator {
        public:

            typedef T value_type;


            typedef T* pointer;


            typedef T& reference;


            typedef const T* const_pointer;


            typedef const T& const_reference;


            typedef ::std::size_t size_type;


            typedef ::std::ptrdiff_t difference_type;

            template <class U>
            struct rebind {

                typedef aligned_allocator<U> other;
            };

            pointer address(reference x) {
                return &x;
            }

            const_pointer address(const_reference x) {
                return &x;
            }

            pointer allocate(
                size_type nb_elt, ::std::allocator<void>::const_pointer hint = nullptr
            ) {
                geo_argused(hint);
                pointer result = static_cast<pointer>(
                    aligned_malloc(sizeof(T) * nb_elt, ALIGN)
                );
                return result;
            }

            void deallocate(pointer p, size_type nb_elt) {
                geo_argused(nb_elt);
                aligned_free(p);
            }

            size_type max_size() const {
                ::std::allocator<char> a;
                return a.max_size() / sizeof(T);
            }

            void construct(pointer p, const_reference val) {
                new (static_cast<void*>(p))value_type(val);
            }

            void destroy(pointer p) {
                p->~value_type();
#ifdef GEO_COMPILER_MSVC
                (void) p; // to avoid a "unreferenced variable" warning
#endif
            }

            template <class T2, int A2> operator aligned_allocator<T2, A2>() {
                return aligned_allocator<T2,A2>();
            }
        };

        template <typename T1, int A1, typename T2, int A2>
        inline bool operator== (
            const aligned_allocator<T1, A1>&, const aligned_allocator<T2, A2>&
        ) {
            return true;
        }

        template <typename T1, int A1, typename T2, int A2>
        inline bool operator!= (
            const aligned_allocator<T1, A1>&, const aligned_allocator<T2, A2>&
        ) {
            return false;
        }
    }



    template <class T>
    class vector : public ::std::vector<T, Memory::aligned_allocator<T> > {
        typedef ::std::vector<T, Memory::aligned_allocator<T> > baseclass;

    public:
        vector() :
            baseclass() {
        }

        explicit vector(index_t size) :
            baseclass(size) {
        }

        explicit vector(index_t size, const T& val) :
            baseclass(size, val) {
        }

        index_t size() const {
            //   casts baseclass::size() from size_t (64 bits)
            //   to index_t (32 bits), because all
            //   indices in Vorpaline are supposed to fit in 32 bits (index_t).
            // TODO: geo_debug_assert(baseclass::size() < max index_t)
            return index_t(baseclass::size());
        }

        T& operator[] (index_t i) {
            geo_debug_assert(i < size());
            return baseclass::operator[] (i);
        }

        const T& operator[] (index_t i) const {
            geo_debug_assert(i < size());
            return baseclass::operator[] (i);
        }

        T& operator[] (signed_index_t i) {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }

        const T& operator[] (signed_index_t i) const {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }


#ifdef GARGANTUA // If compiled with 64 bits index_t

        T& operator[] (int i) {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }

        const T& operator[] (int i) const {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }

        T& operator[] (unsigned int i) {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }

        const T& operator[] (unsigned int i) const {
            geo_debug_assert(i >= 0 && index_t(i) < size());
            return baseclass::operator[] (index_t(i));
        }
#endif

        T* data() {
            return size() == 0 ? nullptr : &(*this)[0];
        }

        const T* data() const {
            return size() == 0 ? nullptr : &(*this)[0];
        }

    };

    template <>
    class vector<bool> : public ::std::vector<bool> {
        typedef ::std::vector<bool> baseclass;

    public:

        vector() :
            baseclass() {
        }


        explicit vector(index_t size) :
            baseclass(size) {
        }


        explicit vector(index_t size, bool val) :
            baseclass(size, val) {
        }


        index_t size() const {
            //   casts baseclass::size() from size_t (64 bits)
            //   to index_t (32 bits), because all
            //   indices in Vorpaline are supposed to fit in 32 bits (index_t).
            // TODO: geo_debug_assert(baseclass::size() < max index_t)
            return index_t(baseclass::size());
        }

        // TODO: operator[] with bounds checking (more complicated
        // than just returning bool&, check implementation in STL).
    };
}

#endif


/******* extracted from ../basic/vecg.h *******/

#ifndef GEOGRAM_BASIC_VECG
#define GEOGRAM_BASIC_VECG


#include <iostream>
#include <cfloat>
#include <cmath>


namespace GEO {

    template <index_t DIM, class T>
    class vecng {
    public:

        static const index_t dim = DIM;


        typedef vecng<DIM, T> vector_type;


        typedef T value_type;

        vecng() {
            for(index_t i = 0; i < DIM; i++) {
                data_[i] = T(0);
            }
        }

        // This one should never be called :
        // a template constructor cannot be a copy constructor

        template <class T2>
        explicit vecng(const vecng<DIM, T2>& v) {
            for(index_t i = 0; i < DIM; i++) {
                data_[i] = T(v[i]);
            }
        }

        // to avoid compilation problems
        template <class T2, index_t DIM2>
        explicit vecng(
            const vecng<DIM2, T2>& v
        ) {
            geo_debug_assert(DIM2 == DIM);
            for(index_t i = 0; i < DIM; i++) {
                data_[i] = T(v[i]);
            }
        }

        template <class T2>
        explicit vecng(const T2* v) {
            for(index_t i = 0; i < DIM; i++) {
                data_[i] = T(v[i]);
            }
        }

        index_t dimension() const {
            return DIM;
        }

        T* data() {
            return data_;
        }

        const T* data() const {
            return data_;
        }

        inline T& operator[] (index_t i) {
            geo_debug_assert(i < DIM);
            return data()[i];
        }

        inline const T& operator[] (index_t i) const {
            geo_debug_assert(i < DIM);
            return data()[i];
        }

        inline T length2() const {
            T result = T(0);
            for(index_t i = 0; i < DIM; i++) {
                result += data_[i] * data_[i];
            }
            return result;
        }

        inline T length() const {
            return sqrt(length2());
        }

        inline T distance2(const vector_type& v) const {
            T result(0);
            for(index_t i = 0; i < DIM; i++) {
                result += geo_sqr(v.data_[i] - data_[i]);
            }
            return result;
        }

        inline T distance(const vector_type& v) const {
            return sqrt(distance2(v));
        }

        // operators

        inline vector_type& operator+= (const vector_type& v) {
            for(index_t i = 0; i < DIM; i++) {
                data_[i] += v.data_[i];
            }
            return *this;
        }

        inline vector_type& operator-= (const vector_type& v) {
            for(index_t i = 0; i < DIM; i++) {
                data_[i] -= v.data_[i];
            }
            return *this;
        }

        template <class T2>
        inline vector_type& operator*= (T2 s) {
            for(index_t i = 0; i < DIM; i++) {
                data_[i] *= T(s);
            }
            return *this;
        }

        template <class T2>
        inline vector_type& operator/= (T2 s) {
            for(index_t i = 0; i < DIM; i++) {
                data_[i] /= T(s);
            }
            return *this;
        }

        inline vector_type operator+ (const vector_type& v) const {
            vector_type result(*this);
            for(index_t i = 0; i < DIM; i++) {
                result.data_[i] += v.data_[i];
            }
            return result;
        }

        inline vector_type operator- (const vector_type& v) const {
            vector_type result(*this);
            for(index_t i = 0; i < DIM; i++) {
                result.data_[i] -= v.data_[i];
            }
            return result;
        }

        template <class T2>
        inline vector_type operator* (T2 s) const {
            vector_type result(*this);
            for(index_t i = 0; i < DIM; i++) {
                result.data_[i] *= T(s);
            }
            return result;
        }

        template <class T2>
        inline vector_type operator/ (T2 s) const {
            vector_type result(*this);
            for(index_t i = 0; i < DIM; i++) {
                result.data_[i] /= T(s);
            }
            return result;
        }

        inline vector_type operator- () const {
            vector_type result;
            for(index_t i = 0; i < DIM; i++) {
                result.data_[i] = -data_[i];
            }
            return result;
        }

    private:
        T data_[DIM];
    };

    template <index_t DIM, class T>
    inline T dot(
        const vecng<DIM, T>& v1, const vecng<DIM, T>& v2
    ) {
        T result = 0;
        for(index_t i = 0; i < DIM; i++) {
            result += v1[i] * v2[i];
        }
        return result;
    }

    template <class T2, index_t DIM, class T>
    inline vecng<DIM, T> operator* (
        T2 s, const vecng<DIM, T>& v
    ) {
        vecng<DIM, T> result;
        for(index_t i = 0; i < DIM; i++) {
            result[i] = T(s) * v[i];
        }
        return result;
    }

    // Compatibility with GLSL

    template <index_t DIM, class T>
    inline T length(const vecng<DIM, T>& v) {
        return v.length();
    }

    template <index_t DIM, class T>
    inline T length2(const vecng<DIM, T>& v) {
        return v.length2();
    }

    template <index_t DIM, class T>
    inline T distance2(
        const vecng<DIM, T>& v1, const vecng<DIM, T>& v2
    ) {
        return v2.distance2(v1);
    }

    template <index_t DIM, class T>
    inline T distance(
        const vecng<DIM, T>& v1, const vecng<DIM, T>& v2
    ) {
        return v2.distance(v1);
    }

    template <index_t DIM, class T>
    inline vecng<DIM, T> normalize(
        const vecng<DIM, T>& v
    ) {
        T s = length(v);
        if(s > 1e-30) {
            s = T(1) / s;
        }
        return s * v;
    }

    template <index_t DIM, class T>
    inline vecng<DIM, T> mix(
        const vecng<DIM, T>& v1, const vecng<DIM, T>& v2, T s
    ) {
        return (T(1) - s) * v1 + s * v2;
    }



    template <class T>
    class vecng<2, T> {
    public:

        static const index_t dim = 2;


        typedef vecng<dim, T> vector_type;


        typedef T value_type;


        vecng() :
            x(0),
            y(0) {
        }

        vecng(T x_in, T y_in) :
            x(x_in),
            y(y_in) {
        }


        template <class T2>
        explicit vecng(const vecng<dim, T2>& v) :
            x(v.x),
            y(v.y) {
        }


        template <class T2>
        explicit vecng(const T2* v) :
            x(v[0]),
            y(v[1]) {
        }


        inline T length2() const {
            return x * x + y * y;
        }


        inline T length() const {
            return sqrt(x * x + y * y);
        }


        inline T distance2(const vector_type& v) const {
            T dx = v.x - x;
            T dy = v.y - y;
            return dx * dx + dy * dy;
        }


        inline T distance(const vector_type& v) const {
            return sqrt(distance2(v));
        }


        inline vector_type& operator+= (const vector_type& v) {
            x += v.x;
            y += v.y;
            return *this;
        }


        inline vector_type& operator-= (const vector_type& v) {
            x -= v.x;
            y -= v.y;
            return *this;
        }


        template <class T2>
        inline vector_type& operator*= (T2 s) {
            x *= T(s);
            y *= T(s);
            return *this;
        }


        template <class T2>
        inline vector_type& operator/= (T2 s) {
            x /= T(s);
            y /= T(s);
            return *this;
        }


        inline vector_type operator+ (const vector_type& v) const {
            return vector_type(x + v.x, y + v.y);
        }


        inline vector_type operator- (const vector_type& v) const {
            return vector_type(x - v.x, y - v.y);
        }


        template <class T2>
        inline vector_type operator* (T2 s) const {
            return vector_type(x * T(s), y * T(s));
        }


        template <class T2>
        inline vector_type operator/ (T2 s) const {
            return vector_type(x / T(s), y / T(s));
        }


        inline vector_type operator- () const {
            return vector_type(-x, -y);
        }


        index_t dimension() const {
            return dim;
        }


        T* data() {
            return &x;
        }


        const T* data() const {
            return &x;
        }


        inline T& operator[] (index_t i) {
            geo_debug_assert(i < dim);
            return data()[i];
        }


        inline const T& operator[] (index_t i) const {
            geo_debug_assert(i < dim);
            return data()[i];
        }


        T x;

        T y;
    };

    template <class T>
    inline T dot(
        const vecng<2, T>& v1, const vecng<2, T>& v2
    ) {
        return v1.x * v2.x + v1.y * v2.y;
    }

    template <class T>
    inline T det(
        const vecng<2, T>& v1, const vecng<2, T>& v2
    ) {
        return v1.x * v2.y - v1.y * v2.x;
    }

    template <class T2, class T>
    inline vecng<2, T> operator* (
        T2 s, const vecng<2, T>& v
    ) {
        return vecng<2, T>(T(s) * v.x, T(s) * v.y);
    }



    template <class T>
    class vecng<3, T> {
    public:

        static const index_t dim = 3;


        typedef vecng<dim, T> vector_type;


        typedef T value_type;


        vecng() :
            x(0),
            y(0),
            z(0) {
        }

        vecng(T x_in, T y_in, T z_in) :
            x(x_in),
            y(y_in),
            z(z_in) {
        }


        template <class T2>
        explicit vecng(const vecng<dim, T2>& v) :
            x(v.x),
            y(v.y),
            z(v.z) {
        }


        template <class T2>
        explicit vecng(const T2* v) :
            x(v[0]),
            y(v[1]),
            z(v[2]) {
        }


        inline T length2() const {
            return x * x + y * y + z * z;
        }


        inline T length() const {
            return sqrt(x * x + y * y + z * z);
        }


        inline T distance2(const vector_type& v) const {
            T dx = v.x - x;
            T dy = v.y - y;
            T dz = v.z - z;
            return dx * dx + dy * dy + dz * dz;
        }


        inline T distance(const vector_type& v) const {
            return sqrt(distance2(v));
        }


        inline vector_type& operator+= (const vector_type& v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }


        inline vector_type& operator-= (const vector_type& v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }


        template <class T2>
        inline vector_type& operator*= (T2 s) {
            x *= T(s);
            y *= T(s);
            z *= T(s);
            return *this;
        }


        template <class T2>
        inline vector_type& operator/= (T2 s) {
            x /= T(s);
            y /= T(s);
            z /= T(s);
            return *this;
        }


        inline vector_type operator+ (const vector_type& v) const {
            return vector_type(x + v.x, y + v.y, z + v.z);
        }


        inline vector_type operator- (const vector_type& v) const {
            return vector_type(x - v.x, y - v.y, z - v.z);
        }


        template <class T2>
        inline vector_type operator* (T2 s) const {
            return vector_type(x * T(s), y * T(s), z * T(s));
        }


        template <class T2>
        inline vector_type operator/ (T2 s) const {
            return vector_type(x / T(s), y / T(s), z / T(s));
        }


        inline vector_type operator- () const {
            return vector_type(-x, -y, -z);
        }


        index_t dimension() const {
            return dim;
        }


        T* data() {
            return &x;
        }


        const T* data() const {
            return &x;
        }


        inline T& operator[] (index_t i) {
            geo_debug_assert(i < dim);
            return data()[i];
        }


        inline const T& operator[] (index_t i) const {
            geo_debug_assert(i < dim);
            return data()[i];
        }


        T x;

        T y;

        T z;
    };

    template <class T>
    inline T dot(
        const vecng<3, T>& v1, const vecng<3, T>& v2
    ) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    template <class T>
    inline vecng<3, T> cross(
        const vecng<3, T>& v1, const vecng<3, T>& v2
    ) {
        return vecng<3, T>(
            v1.y * v2.z - v1.z * v2.y,
            v1.z * v2.x - v1.x * v2.z,
            v1.x * v2.y - v1.y * v2.x
        );
    }

    template <class T2, class T>
    inline vecng<3, T> operator* (
        T2 s, const vecng<3, T>& v
    ) {
        return vecng<3, T>(T(s) * v.x, T(s) * v.y, T(s) * v.z);
    }



    template <class T>
    class vecng<4, T> {
    public:

        static const index_t dim = 4;


        typedef vecng<dim, T> vector_type;


        typedef T value_type;


        vecng() :
            x(0),
            y(0),
            z(0),
            w(0) {
        }

        vecng(T x_in, T y_in, T z_in, T w_in) :
            x(x_in),
            y(y_in),
            z(z_in),
            w(w_in) {
        }


        template <class T2>
        explicit vecng(const vecng<dim, T2>& v) :
            x(v.x),
            y(v.y),
            z(v.z),
            w(v.w) {
        }


        template <class T2>
        explicit vecng(const T2* v) :
            x(v[0]),
            y(v[1]),
            z(v[2]),
            w(v[3]) {
        }


        inline T length2() const {
            return x * x + y * y + z * z + w * w;
        }


        inline T length() const {
            return sqrt(x * x + y * y + z * z + w * w);
        }


        inline T distance2(const vector_type& v) const {
            T dx = v.x - x;
            T dy = v.y - y;
            T dz = v.z - z;
            T dw = v.w - w;
            return dx * dx + dy * dy + dz * dz + dw * dw;
        }


        inline T distance(const vector_type& v) const {
            return sqrt(distance2(v));
        }


        index_t dimension() const {
            return dim;
        }


        inline vector_type& operator+= (const vector_type& v) {
            x += v.x;
            y += v.y;
            z += v.z;
            w += v.w;
            return *this;
        }


        inline vector_type& operator-= (const vector_type& v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            w -= v.w;
            return *this;
        }


        template <class T2>
        inline vector_type& operator*= (T2 s) {
            x *= T(s);
            y *= T(s);
            z *= T(s);
            w *= T(s);
            return *this;
        }


        template <class T2>
        inline vector_type& operator/= (T2 s) {
            x /= T(s);
            y /= T(s);
            z /= T(s);
            w /= T(s);
            return *this;
        }


        inline vector_type operator+ (const vector_type& v) const {
            return vector_type(x + v.x, y + v.y, z + v.z, w + v.w);
        }


        inline vector_type operator- (const vector_type& v) const {
            return vector_type(x - v.x, y - v.y, z - v.z, w - v.w);
        }


        template <class T2>
        inline vector_type operator* (T2 s) const {
            return vector_type(x * T(s), y * T(s), z * T(s), w * T(s));
        }


        template <class T2>
        inline vector_type operator/ (T2 s) const {
            return vector_type(x / T(s), y / T(s), z / T(s), w / T(s));
        }


        inline vector_type operator- () const {
            return vector_type(-x, -y, -z, -w);
        }


        T* data() {
            return &x;
        }


        const T* data() const {
            return &x;
        }


        inline T& operator[] (index_t i) {
            geo_debug_assert(i < dim);
            return data()[i];
        }


        inline const T& operator[] (index_t i) const {
            geo_debug_assert(i < dim);
            return data()[i];
        }


        T x;

        T y;

        T z;

        T w;
    };

    template <class T>
    inline T dot(
        const vecng<4, T>& v1, const vecng<4, T>& v2
    ) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
    }

    template <class T2, class T>
    inline vecng<4, T> operator* (
        T2 s, const vecng<4, T>& v
    ) {
        return vecng<4, T>(T(s) * v.x, T(s) * v.y, T(s) * v.z, T(s) * v.w);
    }

    template <index_t DIM, class T>
    inline std::ostream& operator<< (
        std::ostream& out, const GEO::vecng<DIM, T>& v
    ) {
        const char* sep = "";
        for(index_t i = 0; i < DIM; i++) {
            out << sep << v[i];
            sep = " ";
        }
        return out;
    }

    template <index_t DIM, class T>
    inline std::istream& operator>> (
        std::istream& in, GEO::vecng<DIM, T>& v
    ) {
        for(index_t i = 0; i < DIM; i++) {
            in >> v[i];
        }
        return in;
    }
}

#endif


/******* extracted from ../basic/matrix.h *******/

#ifndef GEOGRAM_BASIC_MATRIX
#define GEOGRAM_BASIC_MATRIX



namespace GEO {



    inline double det2x2(
        double a11, double a12,
        double a21, double a22
    ) {
        return a11*a22-a12*a21 ;
    }

    inline double det3x3(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33
    ) {
    return
         a11*det2x2(a22,a23,a32,a33)
        -a21*det2x2(a12,a13,a32,a33)
        +a31*det2x2(a12,a13,a22,a23);
    }


    inline double det4x4(
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34,
        double a41, double a42, double a43, double a44
    ) {
        double m12 = a21*a12 - a11*a22;
        double m13 = a31*a12 - a11*a32;
        double m14 = a41*a12 - a11*a42;
        double m23 = a31*a22 - a21*a32;
        double m24 = a41*a22 - a21*a42;
        double m34 = a41*a32 - a31*a42;

        double m123 = m23*a13 - m13*a23 + m12*a33;
        double m124 = m24*a13 - m14*a23 + m12*a43;
        double m134 = m34*a13 - m14*a33 + m13*a43;
        double m234 = m34*a23 - m24*a33 + m23*a43;

        return (m234*a14 - m134*a24 + m124*a34 - m123*a44);
    }


    template <index_t DIM, class FT>
    class Matrix {
    public:

        typedef Matrix<DIM, FT> matrix_type;


        typedef FT value_type;


        static const index_t dim = DIM;

        inline Matrix() {
            load_identity();
        }

        explicit Matrix(const FT* vals) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = *vals;
                    ++vals;
                }
            }
        }

        inline index_t dimension() const {
            return DIM;
        }

        inline void load_zero() {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = FT(0);
                }
            }
        }

        inline void load_identity() {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = (i == j) ? FT(1) : FT(0);
                }
            }
        }

        inline bool is_identity() const {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    FT rhs = ((i == j) ? FT(1) : FT(0));
                    if(coeff_[i][j] != rhs) {
                        return false;
                    }
                }
            }
            return true;
        }

        inline FT& operator() (index_t i, index_t j) {
            geo_debug_assert(i < DIM);
            geo_debug_assert(j < DIM);
            return coeff_[i][j];
        }

        inline const FT& operator() (index_t i, index_t j) const {
            geo_debug_assert(i < DIM);
            geo_debug_assert(j < DIM);
            return coeff_[i][j];
        }

        inline matrix_type& operator+= (const matrix_type& m) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] += m.coeff_[i][j];
                }
            }
            return *this;
        }

        inline matrix_type& operator-= (const matrix_type& m) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] -= m.coeff_[i][j];
                }
            }
            return *this;
        }

        inline matrix_type& operator*= (FT val) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] *= val;
                }
            }
            return *this;
        }

        inline matrix_type& operator/= (FT val) {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] /= val;
                }
            }
            return *this;
        }

        inline matrix_type operator+ (const matrix_type& m) const {
            matrix_type result = *this;
            result += m;
            return result;
        }

        inline matrix_type operator- (const matrix_type& m) const {
            matrix_type result = *this;
            result -= m;
            return result;
        }

        inline matrix_type operator* (FT val) const {
            matrix_type result = *this;
            result *= val;
            return result;
        }

        inline matrix_type operator/ (FT val) const {
            matrix_type result = *this;
            result /= val;
            return result;
        }

        matrix_type operator* (const matrix_type& m) const {
            matrix_type result;
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    result.coeff_[i][j] = FT(0);
                    for(index_t k = 0; k < DIM; k++) {
                        result.coeff_[i][j] += coeff_[i][k] * m.coeff_[k][j];
                    }
                }
            }
            return result;
        }

        matrix_type inverse() const {
            matrix_type result;
            bool invertible = compute_inverse(result);
            geo_assert(invertible);
            return result;
        }


        bool compute_inverse(matrix_type& result) const {
            FT val=FT(0.0), val2=FT(0.0);
            matrix_type tmp = (*this);

            result.load_identity();

            for(index_t i = 0; i != DIM; i++) {
                val = tmp(i, i);                     /* find pivot */
                index_t ind = i;
                for(index_t j = i + 1; j != DIM; j++) {
                    if(fabs(tmp(j, i)) > fabs(val)) {
                        ind = j;
                        val = tmp(j, i);
                    }
                }

                if(ind != i) {
                    for(index_t j = 0; j != DIM; j++) {
                        val2 = result(i, j);
                        result(i, j) = result(ind, j);
                        result(ind, j) = val2;           /* swap columns */
                        val2 = tmp(i, j);
                        tmp(i, j) = tmp(ind, j);
                        tmp(ind, j) = val2;
                    }
                }

                if(val == 0.0) {
                    return false;
                }

                for(index_t j = 0; j != DIM; j++) {
                    tmp(i, j) /= val;
                    result(i, j) /= val;
                }

                for(index_t j = 0; j != DIM; j++) {
                    if(j == i) {
                        continue;                       /* eliminate column */
                    }
                    val = tmp(j, i);
                    for(index_t k = 0; k != DIM; k++) {
                        tmp(j, k) -= tmp(i, k) * val;
                        result(j, k) -= result(i, k) * val;
                    }
                }
            }

            return true;
        }

        matrix_type transpose() const {
            matrix_type result;
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j < DIM; j++) {
                    result(i, j) = (* this)(j, i);
                }
            }
            return result;
        }



        inline const FT* data() const {
            return &(coeff_[0][0]);
        }



        inline FT* data() {
            return &(coeff_[0][0]);
        }

        void get_lower_triangle(FT* store) const {
            for(index_t i = 0; i < DIM; i++) {
                for(index_t j = 0; j <= i; j++) {
                    *store++ = coeff_[i][j];
                }
            }
        }

    private:
        FT coeff_[DIM][DIM];
    };



    template <index_t DIM, class FT>
    inline std::ostream& operator<< (
        std::ostream& output, const Matrix<DIM, FT>& m
    ) {
        const char* sep = "";
        for(index_t i = 0; i < DIM; i++) {
            for(index_t j = 0; j < DIM; j++) {
                output << sep << m(i, j);
                sep = " ";
            }
        }
        return output;
    }

    template <index_t DIM, class FT>
    inline std::istream& operator>> (
        std::istream& input, Matrix<DIM, FT>& m
    ) {
        for(index_t i = 0; i < DIM; i++) {
            for(index_t j = 0; j < DIM; j++) {
                input >> m(i, j);
            }
        }
        return input;
    }



    template <index_t DIM, class FT> inline
    void mult(const Matrix<DIM, FT>& M, const FT* x, FT* y) {
        for(index_t i = 0; i < DIM; i++) {
            y[i] = 0;
            for(index_t j = 0; j < DIM; j++) {
                y[i] += M(i, j) * x[j];
            }
        }
    }



    template <index_t DIM, class FT> inline
    vecng<DIM,FT> operator*(
        const Matrix<DIM, FT>& M, const vecng<DIM,FT>& x
    ) {
        vecng<DIM,FT> y;
        for(index_t i = 0; i < DIM; i++) {
            y[i] = 0;
            for(index_t j = 0; j < DIM; j++) {
                y[i] += M(i, j) * x[j];
            }
        }
        return y;
    }

    template <index_t DIM, class FT> inline
    vecng<DIM,FT> mult(
        const Matrix<DIM, FT>& M, const vecng<DIM,FT>& x
    ) {
        vecng<DIM,FT> y;
        for(index_t i = 0; i < DIM; i++) {
            y[i] = 0;
            for(index_t j = 0; j < DIM; j++) {
                y[i] += M(i, j) * x[j];
            }
        }
        return y;
    }



}

#endif


/******* extracted from multi_precision.h *******/

#ifndef GEOGRAM_NUMERICS_MULTI_PRECISION
#define GEOGRAM_NUMERICS_MULTI_PRECISION

#include <iostream>
#include <new>


namespace GEO {

    extern double expansion_splitter_;
    extern double expansion_epsilon_;

    inline void two_sum(double a, double b, double& x, double& y) {
        x = a + b;
        double bvirt = x - a;
        double avirt = x - bvirt;
        double bround = b - bvirt;
        double around = a - avirt;
        y = around + bround;
    }

    inline void two_diff(double a, double b, double& x, double& y) {
        x = a - b;
        double bvirt = a - x;
        double avirt = x + bvirt;
        double bround = bvirt - b;
        double around = a - avirt;
        y = around + bround;
    }

    inline void split(double a, double& ahi, double& alo) {
        double c = expansion_splitter_ * a;
        double abig = c - a;
        ahi = c - abig;
        alo = a - ahi;
    }

    inline void two_product(double a, double b, double& x, double& y) {
#ifdef FP_FAST_FMA
        // If the target processor supports the FMA (Fused Multiply Add)
        // instruction, then the product of two doubles into a length-2
        // expansion can be implemented as follows. Thanks to Marc Glisse
        // for the information.
        // Note: under gcc, automatic generations of fma() for a*b+c needs
        // to be deactivated, using -ffp-contract=off, else it may break
        // other functions such as fast_expansion_sum_zeroelim().
        x = a*b;
        y = fma(a,b,-x);
#else
        x = a * b;
        double ahi, alo;
        split(a, ahi, alo);
        double bhi, blo;
        split(b, bhi, blo);
        double err1 = x - (ahi * bhi);
        double err2 = err1 - (alo * bhi);
        double err3 = err2 - (ahi * blo);
        y = (alo * blo) - err3;
#endif
    }

    inline void square(double a, double& x, double& y) {
#ifdef FP_FAST_FMA
        // If the target processor supports the FMA (Fused Multiply Add)
        // instruction, then the product of two doubles into a length-2
        // expansion can be implemented as follows. Thanks to Marc Glisse
        // for the information.
        // Note: under gcc, automatic generations of fma() for a*b+c needs
        // to be deactivated, using -ffp-contract=off, else it may break
        // other functions such as fast_expansion_sum_zeroelim().
        x = a*a;
        y = fma(a,a,-x);
#else
        x = a * a;
        double ahi, alo;
        split(a, ahi, alo);
        double err1 = x - (ahi * ahi);
        double err3 = err1 - ((ahi + ahi) * alo);
        y = (alo * alo) - err3;
#endif
    }



    class GEOGRAM_API expansion {
    public:
        index_t length() const {
            return length_;
        }

        index_t capacity() const {
            return capacity_;
        }

        void set_length(index_t new_length) {
            geo_debug_assert(new_length <= capacity());
            length_ = new_length;
        }

        const double& operator[] (index_t i) const {
            // Note: we allocate capacity+1 storage
            // systematically, since basic functions
            // may access one additional value (without
            // using it)
            geo_debug_assert(i <= capacity_);
            return x_[i];
        }

        double& operator[] (index_t i) {
            // Note: we allocate capacity+1 storage
            // systematically, since basic functions
            // may access one additional value (without
            // using it)
            geo_debug_assert(i <= capacity_);
            return x_[i];
        }

        double* data() {
            return x_;
        }

        const double* data() const {
            return x_;
        }

        static size_t bytes(index_t capa) {
            // --> 2*sizeof(double) because x_ is declared of size [2]
            // to avoid compiler's warning.
            // --> capa+1 to have an additional 'sentry' at the end
            // because fast_expansion_sum_zeroelim() may access
            // an entry past the end (without using it).
            return
                sizeof(expansion) - 2 * sizeof(double) +
                (capa + 1) * sizeof(double);
        }

        expansion(index_t capa) :
            length_(0),
            capacity_(capa) {
        }

#ifdef CPPCHECK
        // cppcheck does not understand that the result
        // of alloca() is passed to the placement syntax
        // of operator new.
    expansion& new_expansion_on_stack(index_t capa);
#else
#define new_expansion_on_stack(capa)                           \
    (new (alloca(expansion::bytes(capa)))expansion(capa))
#endif

        static expansion* new_expansion_on_heap(index_t capa);

        static void delete_expansion_on_heap(expansion* e);

        // ========================== Initialization from doubles

	expansion& assign(double a) {
	    set_length(1);
	    x_[0] = a;
	    return *this;
	}

        static index_t sum_capacity(double a, double b) {
            geo_argused(a);
            geo_argused(b);
            return 2;
        }

        expansion& assign_sum(double a, double b) {
            set_length(2);
            two_sum(a, b, x_[1], x_[0]);
            return *this;
        }

        static index_t diff_capacity(double a, double b) {
            geo_argused(a);
            geo_argused(b);
            return 2;
        }

        expansion& assign_diff(double a, double b) {
            set_length(2);
            two_diff(a, b, x_[1], x_[0]);
            return *this;
        }

        static index_t product_capacity(double a, double b) {
            geo_argused(a);
            geo_argused(b);
            return 2;
        }

        expansion& assign_product(double a, double b) {
            set_length(2);
            two_product(a, b, x_[1], x_[0]);
            return *this;
        }

        static index_t square_capacity(double a) {
            geo_argused(a);
            return 2;
        }

        expansion& assign_square(double a) {
            set_length(2);
            square(a, x_[1], x_[0]);
            return *this;
        }

        // ====== Initialization from expansion and double

        static index_t sum_capacity(const expansion& a, double b) {
            geo_argused(b);
            return a.length() + 1;
        }

        expansion& assign_sum(const expansion& a, double b);

        static index_t diff_capacity(const expansion& a, double b) {
            geo_argused(b);
            return a.length() + 1;
        }

        expansion& assign_diff(const expansion& a, double b);

        static index_t product_capacity(const expansion& a, double b) {
            geo_argused(b);
            // TODO: implement special case where the double argument
            // is a power of two.
            return a.length() * 2;
        }

        expansion& assign_product(const expansion& a, double b);

        // ========================== Initialization from expansions

        static index_t sum_capacity(const expansion& a, const expansion& b) {
            return a.length() + b.length();
        }

        expansion& assign_sum(const expansion& a, const expansion& b);

        static index_t sum_capacity(
            const expansion& a, const expansion& b, const expansion& c
        ) {
            return a.length() + b.length() + c.length();
        }

        expansion& assign_sum(
            const expansion& a, const expansion& b, const expansion& c
        );

        static index_t sum_capacity(
            const expansion& a, const expansion& b,
            const expansion& c, const expansion& d
        ) {
            return a.length() + b.length() + c.length() + d.length();
        }

        expansion& assign_sum(
            const expansion& a, const expansion& b,
            const expansion& c, const expansion& d
        );

        static index_t diff_capacity(const expansion& a, const expansion& b) {
            return a.length() + b.length();
        }

        expansion& assign_diff(const expansion& a, const expansion& b);

        static index_t product_capacity(
            const expansion& a, const expansion& b
        ) {
            return a.length() * b.length() * 2;
        }

        expansion& assign_product(const expansion& a, const expansion& b);

        static index_t product_capacity(
            const expansion& a, const expansion& b, const expansion& c
        ) {
            return a.length() * b.length() * c.length() * 4;
        }

        expansion& assign_product(
            const expansion& a, const expansion& b, const expansion& c
        );

        static index_t square_capacity(const expansion& a) {
            if(a.length() == 2) {
                return 6;
            }                                  // see two_square()
            return a.length() * a.length() * 2;
        }

        expansion& assign_square(const expansion& a);

        // ====== Determinants =============================

        static index_t det2x2_capacity(
            const expansion& a11, const expansion& a12,
            const expansion& a21, const expansion& a22
        ) {
            return
                product_capacity(a11, a22) +
                product_capacity(a21, a12);
        }

        expansion& assign_det2x2(
            const expansion& a11, const expansion& a12,
            const expansion& a21, const expansion& a22
        );

        static index_t det3x3_capacity(
            const expansion& a11, const expansion& a12, const expansion& a13,
            const expansion& a21, const expansion& a22, const expansion& a23,
            const expansion& a31, const expansion& a32, const expansion& a33
        ) {
            // Development w.r.t. first row
            index_t c11_capa = det2x2_capacity(a22, a23, a32, a33);
            index_t c12_capa = det2x2_capacity(a21, a23, a31, a33);
            index_t c13_capa = det2x2_capacity(a21, a22, a31, a32);
            return 2 * (
                a11.length() * c11_capa +
                a12.length() * c12_capa +
                a13.length() * c13_capa
            );
        }

        expansion& assign_det3x3(
            const expansion& a11, const expansion& a12, const expansion& a13,
            const expansion& a21, const expansion& a22, const expansion& a23,
            const expansion& a31, const expansion& a32, const expansion& a33
        );

        static index_t det_111_2x3_capacity(
            const expansion& a21, const expansion& a22, const expansion& a23,
            const expansion& a31, const expansion& a32, const expansion& a33
        ) {
            return
                det2x2_capacity(a22, a23, a32, a33) +
                det2x2_capacity(a23, a21, a33, a31) +
                det2x2_capacity(a21, a22, a31, a32);
        }

        expansion& assign_det_111_2x3(
            const expansion& a21, const expansion& a22, const expansion& a23,
            const expansion& a31, const expansion& a32, const expansion& a33
        );

        // ======= Geometry-specific initializations =======

        static index_t sq_dist_capacity(coord_index_t dim) {
            return index_t(dim) * 6;
        }

        expansion& assign_sq_dist(
            const double* p1, const double* p2, coord_index_t dim
        );

        static index_t dot_at_capacity(coord_index_t dim) {
            return index_t(dim) * 8;
        }

        expansion& assign_dot_at(
            const double* p1, const double* p2, const double* p0,
            coord_index_t dim
        );


        static index_t length2_capacity(
            const expansion& x, const expansion& y, const expansion& z
        ) {
            return square_capacity(x) + square_capacity(y) + square_capacity(z);
        }

        expansion& assign_length2(
            const expansion& x, const expansion& y, const expansion& z
        );

        // =============== some general purpose functions =========

        static void initialize();

        expansion& negate() {
            for(index_t i = 0; i < length_; ++i) {
                x_[i] = -x_[i];
            }
            return *this;
        }

        expansion& scale_fast(double s) {
            // TODO: debug assert is_power_of_two(s)
            for(index_t i = 0; i < length_; ++i) {
                x_[i] *= s;
            }
            return *this;
        }

        double estimate() const {
            double result = 0.0;
            for(index_t i = 0; i < length(); ++i) {
                result += x_[i];
            }
            return result;
        }

        Sign sign() const {
            if(length() == 0) {
                return ZERO;
            }
            return geo_sgn(x_[length() - 1]);
        }

        std::ostream& show(std::ostream& os) const {
            for(index_t i = 0; i < length(); ++i) {
                os << i << ':' << x_[i] << ' ';
            }
            return os << std::endl;
        }

    protected:
        static index_t sub_product_capacity(
            index_t a_length, index_t b_length
        ) {
            return a_length * b_length * 2;
        }

        expansion& assign_sub_product(
            const double* a, index_t a_length, const expansion& b
        );

#define expansion_sub_product(a, a_length, b)           \
    new_expansion_on_stack(                       \
        sub_product_capacity(a_length, b.length()) \
    )->assign_sub_product(a, a_length, b)

    private:
        expansion(const expansion& rhs);

        expansion& operator= (const expansion& rhs);

    private:
        index_t length_;
        index_t capacity_;
        double x_[2];  // x_ is in fact of size [capacity_]

        friend class expansion_nt;
    };

    // =============== arithmetic operations ===========================

#define expansion_create(a)	      \
    new_expansion_on_stack(1)->assign(a)


#define expansion_sum(a, b)            \
    new_expansion_on_stack(           \
        expansion::sum_capacity(a, b)   \
    )->assign_sum(a, b)

#define expansion_sum3(a, b, c)          \
    new_expansion_on_stack(            \
        expansion::sum_capacity(a, b, c) \
    )->assign_sum(a, b, c)


#define expansion_sum4(a, b, c, d)          \
    new_expansion_on_stack(              \
        expansion::sum_capacity(a, b, c, d) \
    )->assign_sum(a, b, c, d)

#define expansion_diff(a, b)             \
    new_expansion_on_stack(             \
        expansion::diff_capacity(a, b)   \
    )->assign_diff(a, b)

#define expansion_product(a, b)            \
    new_expansion_on_stack(               \
        expansion::product_capacity(a, b)  \
    )->assign_product(a, b)

#define expansion_product3(a, b, c)           \
    new_expansion_on_stack(                 \
        expansion::product_capacity(a, b, c)  \
    )->assign_product(a, b, c)

#define expansion_square(a)             \
    new_expansion_on_stack(             \
        expansion::square_capacity(a)   \
    )->assign_square(a)

    // =============== determinants =====================================

#define expansion_det2x2(a11, a12, a21, a22)          \
    new_expansion_on_stack(                        \
        expansion::det2x2_capacity(a11, a12, a21, a22) \
    )->assign_det2x2(a11, a12, a21, a22)

#define expansion_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33)   \
    new_expansion_on_stack(                                             \
        expansion::det3x3_capacity(a11,a12,a13,a21,a22,a23,a31,a32,a33) \
    )->assign_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33)

#define expansion_det_111_2x3(a21, a22, a23, a31, a32, a33)           \
    new_expansion_on_stack(                                      \
        expansion::det_111_2x3_capacity(a21, a22, a23, a31, a32, a33) \
    )->assign_det_111_2x3(a21, a22, a23, a31, a32, a33)

    // =============== geometric functions ==============================

#define expansion_sq_dist(a, b, dim)           \
    new_expansion_on_stack(                  \
        expansion::sq_dist_capacity(dim)     \
    )->assign_sq_dist(a, b, dim)

#define expansion_dot_at(a, b, c, dim)           \
    new_expansion_on_stack(                   \
        expansion::dot_at_capacity(dim)       \
    )->assign_dot_at(a, b, c, dim)


#define expansion_length2(x,y,z)              \
    new_expansion_on_stack(                   \
       expansion::length2_capacity(x,y,z)     \
    )->assign_length2(x,y,z)



    Sign GEOGRAM_API sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,
        const expansion& a10,const expansion& a11
    );

    Sign GEOGRAM_API sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,const expansion& a02,
        const expansion& a10,const expansion& a11,const expansion& a12,
        const expansion& a20,const expansion& a21,const expansion& a22
    );

    Sign GEOGRAM_API sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,
        const expansion& a02,const expansion& a03,
        const expansion& a10,const expansion& a11,
        const expansion& a12,const expansion& a13,
        const expansion& a20,const expansion& a21,
        const expansion& a22,const expansion& a23,
        const expansion& a30,const expansion& a31,
        const expansion& a32,const expansion& a33
    );


}

#endif


/******* extracted from multi_precision.cpp *******/


// This makes sure the compiler will not optimize y = a*x+b
// with fused multiply-add, this would break the exact
// predicates.
#ifdef GEO_COMPILER_MSVC
#pragma fp_contract(off)
#endif


namespace {

    using namespace GEO;



    bool expansion_length_stat_ = false;
    std::vector<index_t> expansion_length_histo_;

    class ExpansionStatsDisplay {
    public:
        ~ExpansionStatsDisplay() {
            for(index_t i = 0; i < expansion_length_histo_.size(); ++i) {
                std::cerr << "expansion len " << i
                    << " : " << expansion_length_histo_[i] << std::endl;
            }
        }
    };

    ExpansionStatsDisplay expansion_stats_display_;



    class Pools {
    public:

        Pools() : pools_(1024,nullptr) {
            chunks_.reserve(1024);
        }

        ~Pools() {
            for(index_t i=0; i<chunks_.size(); ++i) {
                delete[] chunks_[i];
            }
        }

        void* malloc(size_t size) {
            if(size >= pools_.size()) {
                return ::malloc(size);
            }
            if(pools_[size] == nullptr) {
                new_chunk(size);
            }
            void* result = pools_[size];
            pools_[size] = *static_cast<void**>(pools_[size]);
            return result;
        }

        void free(void* ptr, size_t size) {
            if(size >= pools_.size()) {
                ::free(ptr);
                return;
            }
            *static_cast<void**>(ptr) = pools_[size];
            pools_[size] = ptr;
        }


    protected:
        static const index_t POOL_CHUNK_SIZE = 512;

        void new_chunk(size_t size_in) {
            size_t size = (size_in / 8 + 1)*8; // Align memory.
            Memory::pointer chunk = new Memory::byte[size * POOL_CHUNK_SIZE];
            for(index_t i=0; i<POOL_CHUNK_SIZE-1; ++i) {
                Memory::pointer cur = chunk + size * i;
                Memory::pointer next = cur + size;
                *reinterpret_cast<void**>(cur) = next;
            }
            *reinterpret_cast<void**>(chunk + (size-1)*POOL_CHUNK_SIZE) =
		pools_[size_in];
            pools_[size_in] = chunk;
            chunks_.push_back(chunk);
        }


    private:
        std::vector<void*> pools_;

        std::vector<Memory::pointer> chunks_;

    };

    static Pools pools_;



    inline void fast_two_sum(double a, double b, double& x, double& y) {
        x = a + b;
        double bvirt = x - a;
        y = b - bvirt;
    }

#ifdef REMOVE_ME
    inline void fast_two_diff(double a, double b, double& x, double& y) {
        x = a - b;
        double bvirt = a - x;
        y = bvirt - b;
    }
#endif

    inline void two_one_sum(
        double a1, double a0, double b, double& x2, double& x1, double& x0
    ) {
        double _i;
        two_sum(a0, b, _i, x0);
        two_sum(a1, _i, x2, x1);
    }

    inline void two_two_sum(
        double a1, double a0, double b1, double b0,
        double& x3, double& x2, double& x1, double& x0
    ) {
        double _j, _0;
        two_one_sum(a1, a0, b0, _j, _0, x0);
        two_one_sum(_j, _0, b1, x3, x2, x1);
    }

    inline void two_product_presplit(
        double a, double b, double bhi, double blo, double& x, double& y
    ) {
        x = a * b;
        double ahi;
        double alo;
        split(a, ahi, alo);
        double err1 = x - (ahi * bhi);
        double err2 = err1 - (alo * bhi);
        double err3 = err2 - (ahi * blo);
        y = (alo * blo) - err3;
    }

    inline void two_product_2presplit(
        double a, double ahi, double alo,
        double b, double bhi, double blo,
        double& x, double& y
    ) {
        x = a * b;
        double err1 = x - (ahi * bhi);
        double err2 = err1 - (alo * bhi);
        double err3 = err2 - (ahi * blo);
        y = (alo * blo) - err3;
    }

    inline void two_square(
        double a1, double a0,
        double* x
    ) {
        double _0, _1, _2;
        double _j, _k, _l;
        square(a0, _j, x[0]);
        _0 = a0 + a0;
        two_product(a1, _0, _k, _1);
        two_one_sum(_k, _1, _j, _l, _2, x[1]);
        square(a1, _j, _1);
        two_two_sum(_j, _1, _l, _2, x[5], x[4], x[3], x[2]);
    }

    void two_two_product(
        const double* a,
        const double* b,
        double* x
    ) {
        double _0, _1, _2;
        double _i, _j, _k, _l, _m, _n;

        // If the target processor supports the FMA (Fused Multiply Add)
        // instruction, then the product of two doubles into a length-2
        // expansion can be implemented as follows. Thanks to Marc Glisse
        // for the information.
        // Note: under gcc, automatic generations of fma() for a*b+c needs
        // to be deactivated, using -ffp-contract=off, else it may break
        // other functions such as fast_expansion_sum_zeroelim().
#ifdef FP_FAST_FMA
        two_product(a[0],b[0],_i,x[0]);
        two_product(a[1],b[0],_j,_0);
        two_sum(_i, _0, _k, _1);
        fast_two_sum(_j, _k, _l, _2);
        two_product(a[0], b[1], _i, _0);
        two_sum(_1, _0, _k, x[1]);
        two_sum(_2, _k, _j, _1);
        two_sum(_l, _j, _m, _2);
        two_product(a[1], b[1], _j, _0);
        two_sum(_i, _0, _n, _0);
        two_sum(_1, _0, _i, x[2]);
        two_sum(_2, _i, _k, _1);
        two_sum(_m, _k, _l, _2);
        two_sum(_j, _n, _k, _0);
        two_sum(_1, _0, _j, x[3]);
        two_sum(_2, _j, _i, _1);
        two_sum(_l, _i, _m, _2);
        two_sum(_1, _k, _i, x[4]);
        two_sum(_2, _i, _k, x[5]);
        two_sum(_m, _k, x[7], x[6]);
#else
        double a0hi, a0lo;
        split(a[0], a0hi, a0lo);
        double bhi, blo;
        split(b[0], bhi, blo);
        two_product_2presplit(
            a[0], a0hi, a0lo, b[0], bhi, blo, _i, x[0]
        );
        double a1hi, a1lo;
        split(a[1], a1hi, a1lo);
        two_product_2presplit(
            a[1], a1hi, a1lo, b[0], bhi, blo, _j, _0
        );
        two_sum(_i, _0, _k, _1);
        fast_two_sum(_j, _k, _l, _2);
        split(b[1], bhi, blo);
        two_product_2presplit(
            a[0], a0hi, a0lo, b[1], bhi, blo, _i, _0
        );
        two_sum(_1, _0, _k, x[1]);
        two_sum(_2, _k, _j, _1);
        two_sum(_l, _j, _m, _2);
        two_product_2presplit(
            a[1], a1hi, a1lo, b[1], bhi, blo, _j, _0
        );
        two_sum(_i, _0, _n, _0);
        two_sum(_1, _0, _i, x[2]);
        two_sum(_2, _i, _k, _1);
        two_sum(_m, _k, _l, _2);
        two_sum(_j, _n, _k, _0);
        two_sum(_1, _0, _j, x[3]);
        two_sum(_2, _j, _i, _1);
        two_sum(_l, _i, _m, _2);
        two_sum(_1, _k, _i, x[4]);
        two_sum(_2, _i, _k, x[5]);
        two_sum(_m, _k, x[7], x[6]);
#endif
    }

    void grow_expansion_zeroelim(
        const expansion& e, double b, expansion& h
    ) {
        double Q, hh;
        double Qnew;
        index_t eindex, hindex;
        index_t elen = e.length();

        hindex = 0;
        Q = b;
        for(eindex = 0; eindex < elen; eindex++) {
            double enow = e[eindex];
            two_sum(Q, enow, Qnew, hh);
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    void scale_expansion_zeroelim(
        const expansion& e, double b, expansion& h
    ) {
        double Q, sum;
        double hh;
        double product1;
        double product0;
        index_t eindex, hindex;

        // If the target processor supports the FMA (Fused Multiply Add)
        // instruction, then the product of two doubles into a length-2
        // expansion can be implemented as follows. Thanks to Marc Glisse
        // for the information.
        // Note: under gcc, automatic generations of fma() for a*b+c needs
        // to be deactivated, using -ffp-contract=off, else it may break
        // other functions such as fast_expansion_sum_zeroelim().
#ifndef FP_FAST_FMA
        double bhi, blo;
#endif
        index_t elen = e.length();

        // Sanity check: e and h cannot be the same.
        geo_debug_assert(&e != &h);

#ifdef FP_FAST_FMA
        two_product(e[0], b, Q, hh);
#else
        split(b, bhi, blo);
        two_product_presplit(e[0], b, bhi, blo, Q, hh);
#endif

        hindex = 0;
        if(hh != 0) {
            h[hindex++] = hh;
        }
        for(eindex = 1; eindex < elen; eindex++) {
            double enow = e[eindex];
#ifdef FP_FAST_FMA
            two_product(enow, b,  product1, product0);
#else
            two_product_presplit(enow, b, bhi, blo, product1, product0);
#endif
            two_sum(Q, product0, sum, hh);
            if(hh != 0) {
                h[hindex++] = hh;
            }
            fast_two_sum(product1, sum, Q, hh);
            if(hh != 0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    void fast_expansion_sum_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        double Q;
        double Qnew;
        double hh;
        index_t eindex, findex, hindex;
        double enow, fnow;
        index_t elen = e.length();
        index_t flen = f.length();

        // sanity check: h cannot be e or f
        geo_debug_assert(&h != &e);
        geo_debug_assert(&h != &f);

        enow = e[0];
        fnow = f[0];
        eindex = findex = 0;
        if((fnow > enow) == (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = f[++findex];
        }
        hindex = 0;
        if((eindex < elen) && (findex < flen)) {
            if((fnow > enow) == (fnow > -enow)) {
                fast_two_sum(enow, Q, Qnew, hh);
                enow = e[++eindex];
            } else {
                fast_two_sum(fnow, Q, Qnew, hh);
                fnow = f[++findex];
            }
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
            while((eindex < elen) && (findex < flen)) {
                if((fnow > enow) == (fnow > -enow)) {
                    two_sum(Q, enow, Qnew, hh);
                    enow = e[++eindex];
                } else {
                    two_sum(Q, fnow, Qnew, hh);
                    fnow = f[++findex];
                }
                Q = Qnew;
                if(hh != 0.0) {
                    h[hindex++] = hh;
                }
            }
        }
        while(eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            enow = e[++eindex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        while(findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            fnow = f[++findex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    void fast_expansion_diff_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        double Q;
        double Qnew;
        double hh;
        index_t eindex, findex, hindex;
        double enow, fnow;
        index_t elen = e.length();
        index_t flen = f.length();

        // sanity check: h cannot be e or f
        geo_debug_assert(&h != &e);
        geo_debug_assert(&h != &f);

        enow = e[0];
        fnow = -f[0];
        eindex = findex = 0;
        if((fnow > enow) == (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = -f[++findex];
        }
        hindex = 0;
        if((eindex < elen) && (findex < flen)) {
            if((fnow > enow) == (fnow > -enow)) {
                fast_two_sum(enow, Q, Qnew, hh);
                enow = e[++eindex];
            } else {
                fast_two_sum(fnow, Q, Qnew, hh);
                fnow = -f[++findex];
            }
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
            while((eindex < elen) && (findex < flen)) {
                if((fnow > enow) == (fnow > -enow)) {
                    two_sum(Q, enow, Qnew, hh);
                    enow = e[++eindex];
                } else {
                    two_sum(Q, fnow, Qnew, hh);
                    fnow = -f[++findex];
                }
                Q = Qnew;
                if(hh != 0.0) {
                    h[hindex++] = hh;
                }
            }
        }
        while(eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            enow = e[++eindex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        while(findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            fnow = -f[++findex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }
}



namespace GEO {

    double expansion_splitter_;
    double expansion_epsilon_;

    void expansion::initialize() {
        // Taken from Jonathan Shewchuk's exactinit.
        double half;
        double check, lastcheck;
        int every_other;

        every_other = 1;
        half = 0.5;
        expansion_epsilon_ = 1.0;
        expansion_splitter_ = 1.0;
        check = 1.0;
        // Repeatedly divide `epsilon' by two until it is too small to add to
        // one without causing roundoff.  (Also check if the sum is equal to
        // the previous sum, for machines that round up instead of using exact
        // rounding.  Not that this library will work on such machines anyway.
        do {
            lastcheck = check;
            expansion_epsilon_ *= half;
            if(every_other) {
                expansion_splitter_ *= 2.0;
            }
            every_other = !every_other;
            check = 1.0 + expansion_epsilon_;
        } while((check != 1.0) && (check != lastcheck));
        expansion_splitter_ += 1.0;
    }

    static Process::spinlock expansions_lock = GEOGRAM_SPINLOCK_INIT;

    expansion* expansion::new_expansion_on_heap(index_t capa) {
	Process::acquire_spinlock(expansions_lock);
        if(expansion_length_stat_) {
            if(capa >= expansion_length_histo_.size()) {
                expansion_length_histo_.resize(capa + 1);
            }
            expansion_length_histo_[capa]++;
        }
        Memory::pointer addr = Memory::pointer(
            pools_.malloc(expansion::bytes(capa))
        );
	Process::release_spinlock(expansions_lock);
        expansion* result = new(addr)expansion(capa);
        return result;
    }

    void expansion::delete_expansion_on_heap(expansion* e) {
	Process::acquire_spinlock(expansions_lock);
        pools_.free(e, expansion::bytes(e->capacity()));
	Process::release_spinlock(expansions_lock);
    }

    // ====== Initialization from expansion and double ===============

    expansion& expansion::assign_sum(const expansion& a, double b) {
        geo_debug_assert(capacity() >= sum_capacity(a, b));
        grow_expansion_zeroelim(a, b, *this);
        return *this;
    }

    expansion& expansion::assign_diff(const expansion& a, double b) {
        geo_debug_assert(capacity() >= diff_capacity(a, b));
        grow_expansion_zeroelim(a, -b, *this);
        return *this;
    }

    expansion& expansion::assign_product(const expansion& a, double b) {
        // TODO: implement special case where the double argument
        // is a power of two.
        geo_debug_assert(capacity() >= product_capacity(a, b));
        scale_expansion_zeroelim(a, b, *this);
        return *this;
    }

    // =============  expansion sum and difference =========================

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b));
        fast_expansion_sum_zeroelim(a, b, *this);
        return *this;
    }

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b, c));
        expansion& ab = expansion_sum(a, b);
        this->assign_sum(ab, c);
        return *this;
    }

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b,
        const expansion& c, const expansion& d
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b, c));
        expansion& ab = expansion_sum(a, b);
        expansion& cd = expansion_sum(c, d);
        this->assign_sum(ab, cd);
        return *this;
    }

    expansion& expansion::assign_diff(const expansion& a, const expansion& b) {
        geo_debug_assert(capacity() >= diff_capacity(a, b));
        fast_expansion_diff_zeroelim(a, b, *this);
        return *this;
    }

    // =============  expansion product ==================================

    // Recursive helper function for product implementation
    expansion& expansion::assign_sub_product(
        const double* a, index_t a_length, const expansion& b
    ) {
        geo_debug_assert(
            capacity() >= sub_product_capacity(a_length, b.length())
        );
        if(a_length == 1) {
            scale_expansion_zeroelim(b, a[0], *this);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            const double* a1 = a;
            index_t a1_length = a_length / 2;
            const double* a2 = a1 + a1_length;
            index_t a2_length = a_length - a1_length;
            expansion& a1b = expansion_sub_product(a1, a1_length, b);
            expansion& a2b = expansion_sub_product(a2, a2_length, b);
            this->assign_sum(a1b, a2b);
        }
        return *this;
    }

    expansion& expansion::assign_product(
        const expansion& a, const expansion& b
    ) {
        geo_debug_assert(capacity() >= product_capacity(a, b));
        if(a.length() == 0 || b.length() == 0) {
            x_[0] = 0.0;
            set_length(0);
        } else if(a.length() == 1 && b.length() == 1) {
            two_product(a[0], b[0], x_[1], x_[0]);
            set_length(2);
        } else if(a.length() == 1) {
            scale_expansion_zeroelim(b, a[0], *this);
        } else if(b.length() == 1) {
            scale_expansion_zeroelim(a, b[0], *this);
        } else if(a.length() == 2 && b.length() == 2) {
            two_two_product(a.data(), b.data(), x_);
            set_length(8);
        } else {
            // Recursive distillation: the shortest expansion
            // is split into two parts.
            if(a.length() < b.length()) {
                const double* a1 = a.data();
                index_t a1_length = a.length() / 2;
                const double* a2 = a1 + a1_length;
                index_t a2_length = a.length() - a1_length;
                expansion& a1b = expansion_sub_product(a1, a1_length, b);
                expansion& a2b = expansion_sub_product(a2, a2_length, b);
                this->assign_sum(a1b, a2b);
            } else {
                const double* b1 = b.data();
                index_t b1_length = b.length() / 2;
                const double* b2 = b1 + b1_length;
                index_t b2_length = b.length() - b1_length;
                expansion& ab1 = expansion_sub_product(b1, b1_length, a);
                expansion& ab2 = expansion_sub_product(b2, b2_length, a);
                this->assign_sum(ab1, ab2);
            }
        }
        return *this;
    }

    expansion& expansion::assign_product(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        const expansion& bc = expansion_product(b, c);
        this->assign_product(a, bc);
        return *this;
    }

    expansion& expansion::assign_square(const expansion& a) {
        geo_debug_assert(capacity() >= square_capacity(a));
        if(a.length() == 1) {
            square(a[0], x_[1], x_[0]);
            set_length(2);
        } else if(a.length() == 2) {
            two_square(a[1], a[0], x_);
            set_length(6);
        } else {
            this->assign_product(a, a);
        }
        return *this;
    }

    // =============  determinants ==========================================

    expansion& expansion::assign_det2x2(
        const expansion& a11, const expansion& a12,
        const expansion& a21, const expansion& a22
    ) {
        const expansion& a11a22 = expansion_product(a11, a22);
        const expansion& a12a21 = expansion_product(a12, a21);
        return this->assign_diff(a11a22, a12a21);
    }

    expansion& expansion::assign_det3x3(
        const expansion& a11, const expansion& a12, const expansion& a13,
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        // Development w.r.t. first row
        const expansion& c11 = expansion_det2x2(a22, a23, a32, a33);
        const expansion& c12 = expansion_det2x2(a23, a21, a33, a31);
        const expansion& c13 = expansion_det2x2(a21, a22, a31, a32);
        const expansion& a11c11 = expansion_product(a11, c11);
        const expansion& a12c12 = expansion_product(a12, c12);
        const expansion& a13c13 = expansion_product(a13, c13);
        return this->assign_sum(a11c11, a12c12, a13c13);
    }

    expansion& expansion::assign_det_111_2x3(
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        const expansion& c11 = expansion_det2x2(a22, a23, a32, a33);
        const expansion& c12 = expansion_det2x2(a23, a21, a33, a31);
        const expansion& c13 = expansion_det2x2(a21, a22, a31, a32);
        return this->assign_sum(c11, c12, c13);
    }

    // =============  geometric operations ==================================

    expansion& expansion::assign_sq_dist(
        const double* p1, const double* p2, coord_index_t dim
    ) {
        geo_debug_assert(capacity() >= sq_dist_capacity(dim));
	geo_debug_assert(dim > 0);
        if(dim == 1) {
            double d0, d1;
            two_diff(p1[0], p2[0], d1, d0);
            two_square(d1, d0, x_);
            set_length(6);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            coord_index_t dim1 = dim / 2;
            coord_index_t dim2 = coord_index_t(dim - dim1);
            const double* p1_2 = p1 + dim1;
            const double* p2_2 = p2 + dim1;
            expansion& d1 = expansion_sq_dist(p1, p2, dim1);
            expansion& d2 = expansion_sq_dist(p1_2, p2_2, dim2);
            this->assign_sum(d1, d2);
        }
        return *this;
    }

    expansion& expansion::assign_dot_at(
        const double* p1, const double* p2, const double* p0,
        coord_index_t dim
    ) {
        geo_debug_assert(capacity() >= dot_at_capacity(dim));
        if(dim == 1) {

            double v[2];
            two_diff(p1[0], p0[0], v[1], v[0]);
            double w[2];
            two_diff(p2[0], p0[0], w[1], w[0]);
            two_two_product(v, w, x_);
            set_length(8);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            coord_index_t dim1 = dim / 2;
            coord_index_t dim2 = coord_index_t(dim - dim1);
            const double* p1_2 = p1 + dim1;
            const double* p2_2 = p2 + dim1;
            const double* p0_2 = p0 + dim1;
            expansion& d1 = expansion_dot_at(p1, p2, p0, dim1);
            expansion& d2 = expansion_dot_at(p1_2, p2_2, p0_2, dim2);
            this->assign_sum(d1, d2);
        }
        return *this;
    }

    expansion& expansion::assign_length2(
        const expansion& x, const expansion& y, const expansion& z
    ) {
        const expansion& x2 = expansion_square(x);
        const expansion& y2 = expansion_square(y);
        const expansion& z2 = expansion_square(z);
        this->assign_sum(x2,y2,z2);
        return *this;
    }



    Sign sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,
        const expansion& a10,const expansion& a11
    ) {
        const expansion& result = expansion_det2x2(a00, a01, a10, a11);
        return result.sign();
    }

    Sign sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,const expansion& a02,
        const expansion& a10,const expansion& a11,const expansion& a12,
        const expansion& a20,const expansion& a21,const expansion& a22
    ) {
        // First compute the det2x2
        const expansion& m01 =
            expansion_det2x2(a00, a10, a01, a11);
        const expansion& m02 =
            expansion_det2x2(a00, a20, a01, a21);
        const expansion& m12 =
            expansion_det2x2(a10, a20, a11, a21);

        // Now compute the minors of rank 3
        const expansion& z1 = expansion_product(m01,a22);
        const expansion& z2 = expansion_product(m02,a12).negate();
        const expansion& z3 = expansion_product(m12,a02);

        const expansion& result = expansion_sum3(z1,z2,z3);
        return result.sign();
    }

    Sign sign_of_expansion_determinant(
        const expansion& a00,const expansion& a01,
        const expansion& a02,const expansion& a03,
        const expansion& a10,const expansion& a11,
        const expansion& a12,const expansion& a13,
        const expansion& a20,const expansion& a21,
        const expansion& a22,const expansion& a23,
        const expansion& a30,const expansion& a31,
        const expansion& a32,const expansion& a33
    ) {

        // First compute the det2x2
        const expansion& m01 =
            expansion_det2x2(a10,a00,a11,a01);
        const expansion& m02 =
            expansion_det2x2(a20,a00,a21,a01);
        const expansion& m03 =
            expansion_det2x2(a30,a00,a31,a01);
        const expansion& m12 =
            expansion_det2x2(a20,a10,a21,a11);
        const expansion& m13 =
            expansion_det2x2(a30,a10,a31,a11);
        const expansion& m23 =
            expansion_det2x2(a30,a20,a31,a21);

        // Now compute the minors of rank 3
        const expansion& m012_1 = expansion_product(m12,a02);
        expansion& m012_2 = expansion_product(m02,a12); m012_2.negate();
        const expansion& m012_3 = expansion_product(m01,a22);
        const expansion& m012 = expansion_sum3(m012_1, m012_2, m012_3);

        const expansion& m013_1 = expansion_product(m13,a02);
        expansion& m013_2 = expansion_product(m03,a12); m013_2.negate();

        const expansion& m013_3 = expansion_product(m01,a32);
        const expansion& m013 = expansion_sum3(m013_1, m013_2, m013_3);

        const expansion& m023_1 = expansion_product(m23,a02);
        expansion& m023_2 = expansion_product(m03,a22); m023_2.negate();
        const expansion& m023_3 = expansion_product(m02,a32);
        const expansion& m023 = expansion_sum3(m023_1, m023_2, m023_3);

        const expansion& m123_1 = expansion_product(m23,a12);
        expansion& m123_2 = expansion_product(m13,a22); m123_2.negate();
        const expansion& m123_3 = expansion_product(m12,a32);
        const expansion& m123 = expansion_sum3(m123_1, m123_2, m123_3);

        // Now compute the minors of rank 4
        const expansion& m0123_1 = expansion_product(m123,a03);
        const expansion& m0123_2 = expansion_product(m023,a13);
        const expansion& m0123_3 = expansion_product(m013,a23);
        const expansion& m0123_4 = expansion_product(m012,a33);

        const expansion& z1 = expansion_sum(m0123_1, m0123_3);
        const expansion& z2 = expansion_sum(m0123_2, m0123_4);

        const expansion& result = expansion_diff(z1,z2);
        return result.sign();
    }



}


/******* extracted from predicates/side1.h *******/

inline int side1_3d_filter( const double* p0, const double* p1, const double* q0) {
    double p0_0_p1_0 = (p0[0] - p1[0]);
    double p0_1_p1_1 = (p0[1] - p1[1]);
    double p0_2_p1_2 = (p0[2] - p1[2]);
    double r;
    r = (1 * (((p0_0_p1_0 * p0_0_p1_0) + (p0_1_p1_1 * p0_1_p1_1)) + (p0_2_p1_2 * p0_2_p1_2)));
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    r = (r - (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0_0_p1_0);
    if( (max1 < fabs(p0_1_p1_1)) )
    {
        max1 = fabs(p0_1_p1_1);
    }
    if( (max1 < fabs(p0_2_p1_2)) )
    {
        max1 = fabs(p0_2_p1_2);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    double max2 = fabs(p0_0_p1_0);
    if( (max2 < fabs(p0_1_p1_1)) )
    {
        max2 = fabs(p0_1_p1_1);
    }
    if( (max2 < fabs(p0_2_p1_2)) )
    {
        max2 = fabs(p0_2_p1_2);
    }
    if( (max2 < fabs(q0_0_p0_0)) )
    {
        max2 = fabs(q0_0_p0_0);
    }
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 2.23755023300058943229e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.44425370757048798480e-15 * (max1 * max2));
        if( (r > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


inline int side1_4d_filter( const double* p0, const double* p1, const double* q0) {
    double p0_0_p1_0 = (p0[0] - p1[0]);
    double p0_1_p1_1 = (p0[1] - p1[1]);
    double p0_2_p1_2 = (p0[2] - p1[2]);
    double p0_3_p1_3 = (p0[3] - p1[3]);
    double r;
    r = (1 * ((((p0_0_p1_0 * p0_0_p1_0) + (p0_1_p1_1 * p0_1_p1_1)) + (p0_2_p1_2 * p0_2_p1_2)) + (p0_3_p1_3 * p0_3_p1_3)));
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    r = (r - (2 * ((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p1_0_p0_0);
    if( (max1 < fabs(p0_0_p1_0)) )
    {
        max1 = fabs(p0_0_p1_0);
    }
    if( (max1 < fabs(p0_1_p1_1)) )
    {
        max1 = fabs(p0_1_p1_1);
    }
    if( (max1 < fabs(p0_2_p1_2)) )
    {
        max1 = fabs(p0_2_p1_2);
    }
    if( (max1 < fabs(p0_3_p1_3)) )
    {
        max1 = fabs(p0_3_p1_3);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    double max2 = fabs(p0_0_p1_0);
    if( (max2 < fabs(p0_1_p1_1)) )
    {
        max2 = fabs(p0_1_p1_1);
    }
    if( (max2 < fabs(p0_2_p1_2)) )
    {
        max2 = fabs(p0_2_p1_2);
    }
    if( (max2 < fabs(p0_3_p1_3)) )
    {
        max2 = fabs(p0_3_p1_3);
    }
    if( (max2 < fabs(q0_0_p0_0)) )
    {
        max2 = fabs(q0_0_p0_0);
    }
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (lower_bound_1 < 1.85816790703293534018e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (6.44428177279185717888e-15 * (max1 * max2));
        if( (r > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


inline int side1_6d_filter( const double* p0, const double* p1, const double* q0) {
    double p0_0_p1_0 = (p0[0] - p1[0]);
    double p0_1_p1_1 = (p0[1] - p1[1]);
    double p0_2_p1_2 = (p0[2] - p1[2]);
    double p0_3_p1_3 = (p0[3] - p1[3]);
    double p0_4_p1_4 = (p0[4] - p1[4]);
    double p0_5_p1_5 = (p0[5] - p1[5]);
    double r;
    r = (1 * ((((((p0_0_p1_0 * p0_0_p1_0) + (p0_1_p1_1 * p0_1_p1_1)) + (p0_2_p1_2 * p0_2_p1_2)) + (p0_3_p1_3 * p0_3_p1_3)) + (p0_4_p1_4 * p0_4_p1_4)) + (p0_5_p1_5 * p0_5_p1_5)));
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    r = (r - (2 * ((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0_0_p1_0);
    if( (max1 < fabs(p0_1_p1_1)) )
    {
        max1 = fabs(p0_1_p1_1);
    }
    if( (max1 < fabs(p0_2_p1_2)) )
    {
        max1 = fabs(p0_2_p1_2);
    }
    if( (max1 < fabs(p0_3_p1_3)) )
    {
        max1 = fabs(p0_3_p1_3);
    }
    if( (max1 < fabs(p0_4_p1_4)) )
    {
        max1 = fabs(p0_4_p1_4);
    }
    if( (max1 < fabs(p0_5_p1_5)) )
    {
        max1 = fabs(p0_5_p1_5);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    double max2 = fabs(p0_0_p1_0);
    if( (max2 < fabs(p0_1_p1_1)) )
    {
        max2 = fabs(p0_1_p1_1);
    }
    if( (max2 < fabs(p0_2_p1_2)) )
    {
        max2 = fabs(p0_2_p1_2);
    }
    if( (max2 < fabs(p0_3_p1_3)) )
    {
        max2 = fabs(p0_3_p1_3);
    }
    if( (max2 < fabs(p0_4_p1_4)) )
    {
        max2 = fabs(p0_4_p1_4);
    }
    if( (max2 < fabs(p0_5_p1_5)) )
    {
        max2 = fabs(p0_5_p1_5);
    }
    if( (max2 < fabs(q0_0_p0_0)) )
    {
        max2 = fabs(q0_0_p0_0);
    }
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q0_4_p0_4)) )
    {
        max2 = fabs(q0_4_p0_4);
    }
    if( (max2 < fabs(q0_5_p0_5)) )
    {
        max2 = fabs(q0_5_p0_5);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.41511993781011659868e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.11111223981318615596e-14 * (max1 * max2));
        if( (r > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


inline int side1_7d_filter( const double* p0, const double* p1, const double* q0) {
    double p0_0_p1_0 = (p0[0] - p1[0]);
    double p0_1_p1_1 = (p0[1] - p1[1]);
    double p0_2_p1_2 = (p0[2] - p1[2]);
    double p0_3_p1_3 = (p0[3] - p1[3]);
    double p0_4_p1_4 = (p0[4] - p1[4]);
    double p0_5_p1_5 = (p0[5] - p1[5]);
    double p0_6_p1_6 = (p0[6] - p1[6]);
    double r;
    r = (1 * (((((((p0_0_p1_0 * p0_0_p1_0) + (p0_1_p1_1 * p0_1_p1_1)) + (p0_2_p1_2 * p0_2_p1_2)) + (p0_3_p1_3 * p0_3_p1_3)) + (p0_4_p1_4 * p0_4_p1_4)) + (p0_5_p1_5 * p0_5_p1_5)) + (p0_6_p1_6 * p0_6_p1_6)));
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    r = (r - (2 * (((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0_0_p1_0);
    if( (max1 < fabs(p0_1_p1_1)) )
    {
        max1 = fabs(p0_1_p1_1);
    }
    if( (max1 < fabs(p0_2_p1_2)) )
    {
        max1 = fabs(p0_2_p1_2);
    }
    if( (max1 < fabs(p0_3_p1_3)) )
    {
        max1 = fabs(p0_3_p1_3);
    }
    if( (max1 < fabs(p0_4_p1_4)) )
    {
        max1 = fabs(p0_4_p1_4);
    }
    if( (max1 < fabs(p0_5_p1_5)) )
    {
        max1 = fabs(p0_5_p1_5);
    }
    if( (max1 < fabs(p0_6_p1_6)) )
    {
        max1 = fabs(p0_6_p1_6);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    double max2 = fabs(p0_0_p1_0);
    if( (max2 < fabs(p0_1_p1_1)) )
    {
        max2 = fabs(p0_1_p1_1);
    }
    if( (max2 < fabs(p0_2_p1_2)) )
    {
        max2 = fabs(p0_2_p1_2);
    }
    if( (max2 < fabs(p0_3_p1_3)) )
    {
        max2 = fabs(p0_3_p1_3);
    }
    if( (max2 < fabs(p0_4_p1_4)) )
    {
        max2 = fabs(p0_4_p1_4);
    }
    if( (max2 < fabs(p0_5_p1_5)) )
    {
        max2 = fabs(p0_5_p1_5);
    }
    if( (max2 < fabs(p0_6_p1_6)) )
    {
        max2 = fabs(p0_6_p1_6);
    }
    if( (max2 < fabs(q0_0_p0_0)) )
    {
        max2 = fabs(q0_0_p0_0);
    }
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q0_4_p0_4)) )
    {
        max2 = fabs(q0_4_p0_4);
    }
    if( (max2 < fabs(q0_5_p0_5)) )
    {
        max2 = fabs(q0_5_p0_5);
    }
    if( (max2 < fabs(q0_6_p0_6)) )
    {
        max2 = fabs(q0_6_p0_6);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.27080861580266953580e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.37779349582504943796e-14 * (max1 * max2));
        if( (r > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


inline int side1_8d_filter( const double* p0, const double* p1, const double* q0) {
    double p0_0_p1_0 = (p0[0] - p1[0]);
    double p0_1_p1_1 = (p0[1] - p1[1]);
    double p0_2_p1_2 = (p0[2] - p1[2]);
    double p0_3_p1_3 = (p0[3] - p1[3]);
    double p0_4_p1_4 = (p0[4] - p1[4]);
    double p0_5_p1_5 = (p0[5] - p1[5]);
    double p0_6_p1_6 = (p0[6] - p1[6]);
    double p0_7_p1_7 = (p0[7] - p1[7]);
    double r;
    r = (1 * ((((((((p0_0_p1_0 * p0_0_p1_0) + (p0_1_p1_1 * p0_1_p1_1)) + (p0_2_p1_2 * p0_2_p1_2)) + (p0_3_p1_3 * p0_3_p1_3)) + (p0_4_p1_4 * p0_4_p1_4)) + (p0_5_p1_5 * p0_5_p1_5)) + (p0_6_p1_6 * p0_6_p1_6)) + (p0_7_p1_7 * p0_7_p1_7)));
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    r = (r - (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0_1_p1_1);
    if( (max1 < fabs(p0_2_p1_2)) )
    {
        max1 = fabs(p0_2_p1_2);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p0_0_p1_0)) )
    {
        max1 = fabs(p0_0_p1_0);
    }
    if( (max1 < fabs(p0_3_p1_3)) )
    {
        max1 = fabs(p0_3_p1_3);
    }
    if( (max1 < fabs(p0_4_p1_4)) )
    {
        max1 = fabs(p0_4_p1_4);
    }
    if( (max1 < fabs(p0_5_p1_5)) )
    {
        max1 = fabs(p0_5_p1_5);
    }
    if( (max1 < fabs(p0_6_p1_6)) )
    {
        max1 = fabs(p0_6_p1_6);
    }
    if( (max1 < fabs(p0_7_p1_7)) )
    {
        max1 = fabs(p0_7_p1_7);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    if( (max1 < fabs(p1_7_p0_7)) )
    {
        max1 = fabs(p1_7_p0_7);
    }
    double max2 = fabs(p0_1_p1_1);
    if( (max2 < fabs(p0_2_p1_2)) )
    {
        max2 = fabs(p0_2_p1_2);
    }
    if( (max2 < fabs(p0_0_p1_0)) )
    {
        max2 = fabs(p0_0_p1_0);
    }
    if( (max2 < fabs(p0_3_p1_3)) )
    {
        max2 = fabs(p0_3_p1_3);
    }
    if( (max2 < fabs(p0_4_p1_4)) )
    {
        max2 = fabs(p0_4_p1_4);
    }
    if( (max2 < fabs(p0_5_p1_5)) )
    {
        max2 = fabs(p0_5_p1_5);
    }
    if( (max2 < fabs(p0_6_p1_6)) )
    {
        max2 = fabs(p0_6_p1_6);
    }
    if( (max2 < fabs(p0_7_p1_7)) )
    {
        max2 = fabs(p0_7_p1_7);
    }
    if( (max2 < fabs(q0_0_p0_0)) )
    {
        max2 = fabs(q0_0_p0_0);
    }
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q0_4_p0_4)) )
    {
        max2 = fabs(q0_4_p0_4);
    }
    if( (max2 < fabs(q0_5_p0_5)) )
    {
        max2 = fabs(q0_5_p0_5);
    }
    if( (max2 < fabs(q0_6_p0_6)) )
    {
        max2 = fabs(q0_6_p0_6);
    }
    if( (max2 < fabs(q0_7_p0_7)) )
    {
        max2 = fabs(q0_7_p0_7);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.15542931091530087067e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.66670090166682227006e-14 * (max1 * max2));
        if( (r > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


/******* extracted from predicates/side2.h *******/

inline int side2_3d_filter( const double* p0, const double* p1, const double* p2, const double* q0, const double* q1) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double l1;
    l1 = (1 * (((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double l2;
    l2 = (1 * (((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double a10;
    a10 = (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double a11;
    a11 = (2 * (((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)));
    double a20;
    a20 = (2 * (((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)));
    double a21;
    a21 = (2 * (((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)));
    double Delta;
    Delta = (a11 - a10);
    double DeltaLambda0;
    DeltaLambda0 = (a11 - l1);
    double DeltaLambda1;
    DeltaLambda1 = (l1 - a10);
    double r;
    r = (((Delta * l2) - (a20 * DeltaLambda0)) - (a21 * DeltaLambda1));
    double eps;
    double max1 = fabs(p1_0_p0_0);
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    double max2 = fabs(q0_0_p0_0);
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    }
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    }
    if( (max2 < fabs(q1_2_p0_2)) )
    {
        max2 = fabs(q1_2_p0_2);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 2.23755023300058943229e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.74144419156711063983e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.44425370757048798480e-15 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max3 = max1;
    if( (max3 < max2) )
    {
        max3 = max2;
    }
    double max4 = max2;
    if( (max4 < fabs(p2_0_p0_0)) )
    {
        max4 = fabs(p2_0_p0_0);
    }
    if( (max4 < fabs(p2_1_p0_1)) )
    {
        max4 = fabs(p2_1_p0_1);
    }
    if( (max4 < fabs(p2_2_p0_2)) )
    {
        max4 = fabs(p2_2_p0_2);
    }
    if( (max3 < max4) )
    {
        max3 = max4;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 2.22985945097100191780e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.74144419156711063983e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.99983341597279045654e-14 * (((max1 * max4) * max4) * max3));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side2_4d_filter( const double* p0, const double* p1, const double* p2, const double* q0, const double* q1) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double l1;
    l1 = (1 * ((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double l2;
    l2 = (1 * ((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double a10;
    a10 = (2 * ((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double a11;
    a11 = (2 * ((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)));
    double a20;
    a20 = (2 * ((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)));
    double a21;
    a21 = (2 * ((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)));
    double Delta;
    Delta = (a11 - a10);
    double DeltaLambda0;
    DeltaLambda0 = (a11 - l1);
    double DeltaLambda1;
    DeltaLambda1 = (l1 - a10);
    double r;
    r = (((Delta * l2) - (a20 * DeltaLambda0)) - (a21 * DeltaLambda1));
    double eps;
    double max1 = fabs(p1_2_p0_2);
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    double max2 = fabs(q0_0_p0_0);
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    }
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    }
    if( (max2 < fabs(q1_2_p0_2)) )
    {
        max2 = fabs(q1_2_p0_2);
    }
    if( (max2 < fabs(q1_3_p0_3)) )
    {
        max2 = fabs(q1_3_p0_3);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.85816790703293534018e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (6.44428177279185717888e-15 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max3 = max1;
    if( (max3 < max2) )
    {
        max3 = max2;
    }
    double max4 = max2;
    if( (max4 < fabs(p2_0_p0_0)) )
    {
        max4 = fabs(p2_0_p0_0);
    }
    if( (max4 < fabs(p2_2_p0_2)) )
    {
        max4 = fabs(p2_2_p0_2);
    }
    if( (max4 < fabs(p2_1_p0_1)) )
    {
        max4 = fabs(p2_1_p0_1);
    }
    if( (max4 < fabs(p2_3_p0_3)) )
    {
        max4 = fabs(p2_3_p0_3);
    }
    if( (max3 < max4) )
    {
        max3 = max4;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 1.89528395402941802921e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.72443682410932010423e-13 * (((max1 * max4) * max4) * max3));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side2_6d_filter( const double* p0, const double* p1, const double* p2, const double* q0, const double* q1) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double l1;
    l1 = (1 * ((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double l2;
    l2 = (1 * ((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double a10;
    a10 = (2 * ((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double a11;
    a11 = (2 * ((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)));
    double a20;
    a20 = (2 * ((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)));
    double a21;
    a21 = (2 * ((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)));
    double Delta;
    Delta = (a11 - a10);
    double DeltaLambda0;
    DeltaLambda0 = (a11 - l1);
    double DeltaLambda1;
    DeltaLambda1 = (l1 - a10);
    double r;
    r = (((Delta * l2) - (a20 * DeltaLambda0)) - (a21 * DeltaLambda1));
    double eps;
    double max1 = fabs(p1_1_p0_1);
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    double max2 = fabs(q0_0_p0_0);
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q0_4_p0_4)) )
    {
        max2 = fabs(q0_4_p0_4);
    }
    if( (max2 < fabs(q0_5_p0_5)) )
    {
        max2 = fabs(q0_5_p0_5);
    }
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    }
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    }
    if( (max2 < fabs(q1_2_p0_2)) )
    {
        max2 = fabs(q1_2_p0_2);
    }
    if( (max2 < fabs(q1_3_p0_3)) )
    {
        max2 = fabs(q1_3_p0_3);
    }
    if( (max2 < fabs(q1_4_p0_4)) )
    {
        max2 = fabs(q1_4_p0_4);
    }
    if( (max2 < fabs(q1_5_p0_5)) )
    {
        max2 = fabs(q1_5_p0_5);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.41511993781011659868e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.11111223981318615596e-14 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max3 = max1;
    if( (max3 < max2) )
    {
        max3 = max2;
    }
    double max4 = max2;
    if( (max4 < fabs(p2_2_p0_2)) )
    {
        max4 = fabs(p2_2_p0_2);
    }
    if( (max4 < fabs(p2_4_p0_4)) )
    {
        max4 = fabs(p2_4_p0_4);
    }
    if( (max4 < fabs(p2_3_p0_3)) )
    {
        max4 = fabs(p2_3_p0_3);
    }
    if( (max4 < fabs(p2_0_p0_0)) )
    {
        max4 = fabs(p2_0_p0_0);
    }
    if( (max4 < fabs(p2_1_p0_1)) )
    {
        max4 = fabs(p2_1_p0_1);
    }
    if( (max4 < fabs(p2_5_p0_5)) )
    {
        max4 = fabs(p2_5_p0_5);
    }
    if( (max3 < max4) )
    {
        max3 = max4;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 1.49958502193059513986e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.40007476026584016994e-13 * (((max1 * max4) * max4) * max3));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side2_7d_filter( const double* p0, const double* p1, const double* p2, const double* q0, const double* q1) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double l1;
    l1 = (1 * (((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double l2;
    l2 = (1 * (((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double a10;
    a10 = (2 * (((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double a11;
    a11 = (2 * (((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)));
    double a20;
    a20 = (2 * (((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)));
    double a21;
    a21 = (2 * (((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)));
    double Delta;
    Delta = (a11 - a10);
    double DeltaLambda0;
    DeltaLambda0 = (a11 - l1);
    double DeltaLambda1;
    DeltaLambda1 = (l1 - a10);
    double r;
    r = (((Delta * l2) - (a20 * DeltaLambda0)) - (a21 * DeltaLambda1));
    double eps;
    double max1 = fabs(p1_2_p0_2);
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    double max2 = fabs(q0_0_p0_0);
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q0_4_p0_4)) )
    {
        max2 = fabs(q0_4_p0_4);
    }
    if( (max2 < fabs(q0_5_p0_5)) )
    {
        max2 = fabs(q0_5_p0_5);
    }
    if( (max2 < fabs(q0_6_p0_6)) )
    {
        max2 = fabs(q0_6_p0_6);
    }
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    }
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    }
    if( (max2 < fabs(q1_2_p0_2)) )
    {
        max2 = fabs(q1_2_p0_2);
    }
    if( (max2 < fabs(q1_3_p0_3)) )
    {
        max2 = fabs(q1_3_p0_3);
    }
    if( (max2 < fabs(q1_4_p0_4)) )
    {
        max2 = fabs(q1_4_p0_4);
    }
    if( (max2 < fabs(q1_5_p0_5)) )
    {
        max2 = fabs(q1_5_p0_5);
    }
    if( (max2 < fabs(q1_6_p0_6)) )
    {
        max2 = fabs(q1_6_p0_6);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.27080861580266953580e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.37779349582504943796e-14 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max3 = max1;
    if( (max3 < max2) )
    {
        max3 = max2;
    }
    double max4 = max2;
    if( (max4 < fabs(p2_0_p0_0)) )
    {
        max4 = fabs(p2_0_p0_0);
    }
    if( (max4 < fabs(p2_1_p0_1)) )
    {
        max4 = fabs(p2_1_p0_1);
    }
    if( (max4 < fabs(p2_2_p0_2)) )
    {
        max4 = fabs(p2_2_p0_2);
    }
    if( (max4 < fabs(p2_3_p0_3)) )
    {
        max4 = fabs(p2_3_p0_3);
    }
    if( (max4 < fabs(p2_4_p0_4)) )
    {
        max4 = fabs(p2_4_p0_4);
    }
    if( (max4 < fabs(p2_5_p0_5)) )
    {
        max4 = fabs(p2_5_p0_5);
    }
    if( (max4 < fabs(p2_6_p0_6)) )
    {
        max4 = fabs(p2_6_p0_6);
    }
    if( (max3 < max4) )
    {
        max3 = max4;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 1.36918881183883509035e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (6.33127335329798996022e-13 * (((max1 * max4) * max4) * max3));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side2_8d_filter( const double* p0, const double* p1, const double* p2, const double* q0, const double* q1) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double l1;
    l1 = (1 * ((((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)) + (p1_7_p0_7 * p1_7_p0_7)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double p2_7_p0_7 = (p2[7] - p0[7]);
    double l2;
    l2 = (1 * ((((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)) + (p2_7_p0_7 * p2_7_p0_7)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    double a10;
    a10 = (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double q1_7_p0_7 = (q1[7] - p0[7]);
    double a11;
    a11 = (2 * ((((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)) + (p1_7_p0_7 * q1_7_p0_7)));
    double a20;
    a20 = (2 * ((((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)) + (p2_7_p0_7 * q0_7_p0_7)));
    double a21;
    a21 = (2 * ((((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)) + (p2_7_p0_7 * q1_7_p0_7)));
    double Delta;
    Delta = (a11 - a10);
    double DeltaLambda0;
    DeltaLambda0 = (a11 - l1);
    double DeltaLambda1;
    DeltaLambda1 = (l1 - a10);
    double r;
    r = (((Delta * l2) - (a20 * DeltaLambda0)) - (a21 * DeltaLambda1));
    double eps;
    double max1 = fabs(p1_4_p0_4);
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_7_p0_7)) )
    {
        max1 = fabs(p1_7_p0_7);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    double max2 = fabs(q0_0_p0_0);
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q0_2_p0_2)) )
    {
        max2 = fabs(q0_2_p0_2);
    }
    if( (max2 < fabs(q0_3_p0_3)) )
    {
        max2 = fabs(q0_3_p0_3);
    }
    if( (max2 < fabs(q0_4_p0_4)) )
    {
        max2 = fabs(q0_4_p0_4);
    }
    if( (max2 < fabs(q0_5_p0_5)) )
    {
        max2 = fabs(q0_5_p0_5);
    }
    if( (max2 < fabs(q0_6_p0_6)) )
    {
        max2 = fabs(q0_6_p0_6);
    }
    if( (max2 < fabs(q0_7_p0_7)) )
    {
        max2 = fabs(q0_7_p0_7);
    }
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    }
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    }
    if( (max2 < fabs(q1_2_p0_2)) )
    {
        max2 = fabs(q1_2_p0_2);
    }
    if( (max2 < fabs(q1_3_p0_3)) )
    {
        max2 = fabs(q1_3_p0_3);
    }
    if( (max2 < fabs(q1_4_p0_4)) )
    {
        max2 = fabs(q1_4_p0_4);
    }
    if( (max2 < fabs(q1_5_p0_5)) )
    {
        max2 = fabs(q1_5_p0_5);
    }
    if( (max2 < fabs(q1_6_p0_6)) )
    {
        max2 = fabs(q1_6_p0_6);
    }
    if( (max2 < fabs(q1_7_p0_7)) )
    {
        max2 = fabs(q1_7_p0_7);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (lower_bound_1 < 1.15542931091530087067e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.66670090166682227006e-14 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max3 = max2;
    if( (max3 < max1) )
    {
        max3 = max1;
    }
    double max4 = max2;
    if( (max4 < fabs(p2_4_p0_4)) )
    {
        max4 = fabs(p2_4_p0_4);
    }
    if( (max4 < fabs(p2_2_p0_2)) )
    {
        max4 = fabs(p2_2_p0_2);
    }
    if( (max4 < fabs(p2_0_p0_0)) )
    {
        max4 = fabs(p2_0_p0_0);
    }
    if( (max4 < fabs(p2_1_p0_1)) )
    {
        max4 = fabs(p2_1_p0_1);
    }
    if( (max4 < fabs(p2_3_p0_3)) )
    {
        max4 = fabs(p2_3_p0_3);
    }
    if( (max4 < fabs(p2_5_p0_5)) )
    {
        max4 = fabs(p2_5_p0_5);
    }
    if( (max4 < fabs(p2_6_p0_6)) )
    {
        max4 = fabs(p2_6_p0_6);
    }
    if( (max4 < fabs(p2_7_p0_7)) )
    {
        max4 = fabs(p2_7_p0_7);
    }
    if( (max3 < max4) )
    {
        max3 = max4;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    if( (lower_bound_1 < 1.26419510663115923609e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355531992e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.71140112255785451890e-13 * (((max1 * max4) * max4) * max3));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


/******* extracted from predicates/side3.h *******/

inline int side3_2d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* q0, const double* q1, const double* q2) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double l1;
    l1 = (1 * ((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double l2;
    l2 = (1 * ((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double l3;
    l3 = (1 * ((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double a10;
    a10 = (2 * ((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double a11;
    a11 = (2 * ((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double a12;
    a12 = (2 * ((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)));
    double a20;
    a20 = (2 * ((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)));
    double a21;
    a21 = (2 * ((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)));
    double a22;
    a22 = (2 * ((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)));
    double a30;
    a30 = (2 * ((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)));
    double a31;
    a31 = (2 * ((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)));
    double a32;
    a32 = (2 * ((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)));
    double b00;
    b00 = ((a11 * a22) - (a12 * a21));
    double b01;
    b01 = (a21 - a22);
    double b02;
    b02 = (a12 - a11);
    double b10;
    b10 = ((a12 * a20) - (a10 * a22));
    double b11;
    b11 = (a22 - a20);
    double b12;
    b12 = (a10 - a12);
    double b20;
    b20 = ((a10 * a21) - (a11 * a20));
    double b21;
    b21 = (a20 - a21);
    double b22;
    b22 = (a11 - a10);
    double Delta;
    Delta = ((b00 + b10) + b20);
    double DeltaLambda0;
    DeltaLambda0 = (((b01 * l1) + (b02 * l2)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = (((b11 * l1) + (b12 * l2)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = (((b21 * l1) + (b22 * l2)) + b20);
    double r;
    r = ((Delta * l3) - (((a30 * DeltaLambda0) + (a31 * DeltaLambda1)) + (a32 * DeltaLambda2)));
    double eps;
    double max1 = fabs(p2_0_p0_0);
    if( (max1 < fabs(p2_1_p0_1)) )
    {
        max1 = fabs(p2_1_p0_1);
    }
    double max2 = fabs(q0_0_p0_0);
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    }
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    }
    if( (max2 < fabs(q2_0_p0_0)) )
    {
        max2 = fabs(q2_0_p0_0);
    }
    if( (max2 < fabs(q2_1_p0_1)) )
    {
        max2 = fabs(q2_1_p0_1);
    }
    double max3 = fabs(p1_0_p0_0);
    if( (max3 < fabs(p1_1_p0_1)) )
    {
        max3 = fabs(p1_1_p0_1);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 2.79532528033945620759e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 2.59614842926741294957e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (3.64430756537603111258e-14 * (((max3 * max2) * max1) * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max4 = max2;
    if( (max4 < max1) )
    {
        max4 = max1;
    }
    if( (max4 < max3) )
    {
        max4 = max3;
    }
    double max5 = max2;
    if( (max5 < max3) )
    {
        max5 = max3;
    }
    double max6 = max2;
    if( (max6 < max3) )
    {
        max6 = max3;
    }
    if( (max5 < max6) )
    {
        max5 = max6;
    }
    double max7 = max3;
    if( (max7 < fabs(p3_1_p0_1)) )
    {
        max7 = fabs(p3_1_p0_1);
    }
    if( (max7 < fabs(p3_0_p0_0)) )
    {
        max7 = fabs(p3_0_p0_0);
    }
    if( (max5 < max7) )
    {
        max5 = max7;
    }
    if( (max4 < max5) )
    {
        max4 = max5;
    }
    if( (max4 < max6) )
    {
        max4 = max6;
    }
    if( (max4 < max7) )
    {
        max4 = max7;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    if( (lower_bound_1 < 6.01986729486167248087e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 2.59614842926741294957e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.67544471613800658534e-13 * (((((max7 * max2) * max1) * max6) * max5) * max4));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side3_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* q0, const double* q1, const double* q2) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double l1;
    l1 = (1 * (((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double l2;
    l2 = (1 * (((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double l3;
    l3 = (1 * (((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double a10;
    a10 = (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double a11;
    a11 = (2 * (((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double a12;
    a12 = (2 * (((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)));
    double a20;
    a20 = (2 * (((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)));
    double a21;
    a21 = (2 * (((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)));
    double a22;
    a22 = (2 * (((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)));
    double a30;
    a30 = (2 * (((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)));
    double a31;
    a31 = (2 * (((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)));
    double a32;
    a32 = (2 * (((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)));
    double b00;
    b00 = ((a11 * a22) - (a12 * a21));
    double b01;
    b01 = (a21 - a22);
    double b02;
    b02 = (a12 - a11);
    double b10;
    b10 = ((a12 * a20) - (a10 * a22));
    double b11;
    b11 = (a22 - a20);
    double b12;
    b12 = (a10 - a12);
    double b20;
    b20 = ((a10 * a21) - (a11 * a20));
    double b21;
    b21 = (a20 - a21);
    double b22;
    b22 = (a11 - a10);
    double Delta;
    Delta = ((b00 + b10) + b20);
    double DeltaLambda0;
    DeltaLambda0 = (((b01 * l1) + (b02 * l2)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = (((b11 * l1) + (b12 * l2)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = (((b21 * l1) + (b22 * l2)) + b20);
    double r;
    r = ((Delta * l3) - (((a30 * DeltaLambda0) + (a31 * DeltaLambda1)) + (a32 * DeltaLambda2)));
    double eps;
    double max1 = fabs(p1_1_p0_1);
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    double max2 = fabs(q0_2_p0_2);
    if( (max2 < fabs(q0_0_p0_0)) )
    {
        max2 = fabs(q0_0_p0_0);
    }
    if( (max2 < fabs(q0_1_p0_1)) )
    {
        max2 = fabs(q0_1_p0_1);
    }
    if( (max2 < fabs(q1_0_p0_0)) )
    {
        max2 = fabs(q1_0_p0_0);
    }
    if( (max2 < fabs(q1_1_p0_1)) )
    {
        max2 = fabs(q1_1_p0_1);
    }
    if( (max2 < fabs(q1_2_p0_2)) )
    {
        max2 = fabs(q1_2_p0_2);
    }
    if( (max2 < fabs(q2_0_p0_0)) )
    {
        max2 = fabs(q2_0_p0_0);
    }
    if( (max2 < fabs(q2_1_p0_1)) )
    {
        max2 = fabs(q2_1_p0_1);
    }
    if( (max2 < fabs(q2_2_p0_2)) )
    {
        max2 = fabs(q2_2_p0_2);
    }
    double max3 = fabs(p2_2_p0_2);
    if( (max3 < fabs(p2_0_p0_0)) )
    {
        max3 = fabs(p2_0_p0_0);
    }
    if( (max3 < fabs(p2_1_p0_1)) )
    {
        max3 = fabs(p2_1_p0_1);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 2.22985945097100191780e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 2.59614842926741294957e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.99983341597279045654e-14 * (((max1 * max2) * max3) * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max4 = max1;
    if( (max4 < max2) )
    {
        max4 = max2;
    }
    if( (max4 < max3) )
    {
        max4 = max3;
    }
    double max5 = max1;
    if( (max5 < max2) )
    {
        max5 = max2;
    }
    double max6 = max1;
    if( (max6 < fabs(p3_0_p0_0)) )
    {
        max6 = fabs(p3_0_p0_0);
    }
    if( (max6 < fabs(p3_1_p0_1)) )
    {
        max6 = fabs(p3_1_p0_1);
    }
    if( (max6 < fabs(p3_2_p0_2)) )
    {
        max6 = fabs(p3_2_p0_2);
    }
    if( (max5 < max6) )
    {
        max5 = max6;
    }
    double max7 = max1;
    if( (max7 < max2) )
    {
        max7 = max2;
    }
    if( (max5 < max7) )
    {
        max5 = max7;
    }
    if( (max4 < max5) )
    {
        max4 = max5;
    }
    if( (max4 < max6) )
    {
        max4 = max6;
    }
    if( (max4 < max7) )
    {
        max4 = max7;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    if( (lower_bound_1 < 4.84416636653081796592e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 2.59614842926741294957e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.72198804259438718181e-12 * (((((max6 * max2) * max3) * max7) * max5) * max4));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side3_4d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* q0, const double* q1, const double* q2) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double l1;
    l1 = (1 * ((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double l2;
    l2 = (1 * ((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double l3;
    l3 = (1 * ((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double a10;
    a10 = (2 * ((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double a11;
    a11 = (2 * ((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double a12;
    a12 = (2 * ((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)));
    double a20;
    a20 = (2 * ((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)));
    double a21;
    a21 = (2 * ((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)));
    double a22;
    a22 = (2 * ((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)));
    double a30;
    a30 = (2 * ((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)));
    double a31;
    a31 = (2 * ((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)));
    double a32;
    a32 = (2 * ((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)));
    double b00;
    b00 = ((a11 * a22) - (a12 * a21));
    double b01;
    b01 = (a21 - a22);
    double b02;
    b02 = (a12 - a11);
    double b10;
    b10 = ((a12 * a20) - (a10 * a22));
    double b11;
    b11 = (a22 - a20);
    double b12;
    b12 = (a10 - a12);
    double b20;
    b20 = ((a10 * a21) - (a11 * a20));
    double b21;
    b21 = (a20 - a21);
    double b22;
    b22 = (a11 - a10);
    double Delta;
    Delta = ((b00 + b10) + b20);
    double DeltaLambda0;
    DeltaLambda0 = (((b01 * l1) + (b02 * l2)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = (((b11 * l1) + (b12 * l2)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = (((b21 * l1) + (b22 * l2)) + b20);
    double r;
    r = ((Delta * l3) - (((a30 * DeltaLambda0) + (a31 * DeltaLambda1)) + (a32 * DeltaLambda2)));
    double eps;
    double max1 = fabs(p1_3_p0_3);
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    double max2 = fabs(p2_3_p0_3);
    if( (max2 < fabs(p2_2_p0_2)) )
    {
        max2 = fabs(p2_2_p0_2);
    }
    if( (max2 < fabs(p2_0_p0_0)) )
    {
        max2 = fabs(p2_0_p0_0);
    }
    if( (max2 < fabs(p2_1_p0_1)) )
    {
        max2 = fabs(p2_1_p0_1);
    }
    double max3 = fabs(q0_1_p0_1);
    if( (max3 < fabs(q0_0_p0_0)) )
    {
        max3 = fabs(q0_0_p0_0);
    }
    if( (max3 < fabs(q0_2_p0_2)) )
    {
        max3 = fabs(q0_2_p0_2);
    }
    if( (max3 < fabs(q0_3_p0_3)) )
    {
        max3 = fabs(q0_3_p0_3);
    }
    if( (max3 < fabs(q1_0_p0_0)) )
    {
        max3 = fabs(q1_0_p0_0);
    }
    if( (max3 < fabs(q1_1_p0_1)) )
    {
        max3 = fabs(q1_1_p0_1);
    }
    if( (max3 < fabs(q1_2_p0_2)) )
    {
        max3 = fabs(q1_2_p0_2);
    }
    if( (max3 < fabs(q1_3_p0_3)) )
    {
        max3 = fabs(q1_3_p0_3);
    }
    if( (max3 < fabs(q2_0_p0_0)) )
    {
        max3 = fabs(q2_0_p0_0);
    }
    if( (max3 < fabs(q2_1_p0_1)) )
    {
        max3 = fabs(q2_1_p0_1);
    }
    if( (max3 < fabs(q2_2_p0_2)) )
    {
        max3 = fabs(q2_2_p0_2);
    }
    if( (max3 < fabs(q2_3_p0_3)) )
    {
        max3 = fabs(q2_3_p0_3);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 1.89528395402941802921e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.72443682410931985179e-13 * (((max1 * max3) * max2) * max3));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max4 = max1;
    double max5 = max1;
    double max6 = max1;
    if( (max6 < fabs(p3_0_p0_0)) )
    {
        max6 = fabs(p3_0_p0_0);
    }
    if( (max6 < fabs(p3_3_p0_3)) )
    {
        max6 = fabs(p3_3_p0_3);
    }
    if( (max6 < fabs(p3_2_p0_2)) )
    {
        max6 = fabs(p3_2_p0_2);
    }
    if( (max6 < fabs(p3_1_p0_1)) )
    {
        max6 = fabs(p3_1_p0_1);
    }
    if( (max5 < max6) )
    {
        max5 = max6;
    }
    if( (max5 < max3) )
    {
        max5 = max3;
    }
    double max7 = max1;
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    if( (max5 < max7) )
    {
        max5 = max7;
    }
    if( (max4 < max5) )
    {
        max4 = max5;
    }
    if( (max4 < max6) )
    {
        max4 = max6;
    }
    if( (max4 < max2) )
    {
        max4 = max2;
    }
    if( (max4 < max3) )
    {
        max4 = max3;
    }
    if( (max4 < max7) )
    {
        max4 = max7;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max5;
    upper_bound_1 = max5;
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    if( (lower_bound_1 < 4.14607644401726239868e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.38046888801178809320e-12 * (((((max6 * max3) * max2) * max7) * max5) * max4));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side3_6d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* q0, const double* q1, const double* q2) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double l1;
    l1 = (1 * ((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double l2;
    l2 = (1 * ((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double l3;
    l3 = (1 * ((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double a10;
    a10 = (2 * ((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double a11;
    a11 = (2 * ((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double a12;
    a12 = (2 * ((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)));
    double a20;
    a20 = (2 * ((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)));
    double a21;
    a21 = (2 * ((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)));
    double a22;
    a22 = (2 * ((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)));
    double a30;
    a30 = (2 * ((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)));
    double a31;
    a31 = (2 * ((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)));
    double a32;
    a32 = (2 * ((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)));
    double b00;
    b00 = ((a11 * a22) - (a12 * a21));
    double b01;
    b01 = (a21 - a22);
    double b02;
    b02 = (a12 - a11);
    double b10;
    b10 = ((a12 * a20) - (a10 * a22));
    double b11;
    b11 = (a22 - a20);
    double b12;
    b12 = (a10 - a12);
    double b20;
    b20 = ((a10 * a21) - (a11 * a20));
    double b21;
    b21 = (a20 - a21);
    double b22;
    b22 = (a11 - a10);
    double Delta;
    Delta = ((b00 + b10) + b20);
    double DeltaLambda0;
    DeltaLambda0 = (((b01 * l1) + (b02 * l2)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = (((b11 * l1) + (b12 * l2)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = (((b21 * l1) + (b22 * l2)) + b20);
    double r;
    r = ((Delta * l3) - (((a30 * DeltaLambda0) + (a31 * DeltaLambda1)) + (a32 * DeltaLambda2)));
    double eps;
    double max1 = fabs(p1_0_p0_0);
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    double max2 = fabs(p2_0_p0_0);
    if( (max2 < fabs(p2_1_p0_1)) )
    {
        max2 = fabs(p2_1_p0_1);
    }
    if( (max2 < fabs(p2_2_p0_2)) )
    {
        max2 = fabs(p2_2_p0_2);
    }
    if( (max2 < fabs(p2_3_p0_3)) )
    {
        max2 = fabs(p2_3_p0_3);
    }
    if( (max2 < fabs(p2_4_p0_4)) )
    {
        max2 = fabs(p2_4_p0_4);
    }
    if( (max2 < fabs(p2_5_p0_5)) )
    {
        max2 = fabs(p2_5_p0_5);
    }
    double max3 = fabs(q0_0_p0_0);
    if( (max3 < fabs(q0_1_p0_1)) )
    {
        max3 = fabs(q0_1_p0_1);
    }
    if( (max3 < fabs(q0_2_p0_2)) )
    {
        max3 = fabs(q0_2_p0_2);
    }
    if( (max3 < fabs(q0_3_p0_3)) )
    {
        max3 = fabs(q0_3_p0_3);
    }
    if( (max3 < fabs(q0_4_p0_4)) )
    {
        max3 = fabs(q0_4_p0_4);
    }
    if( (max3 < fabs(q0_5_p0_5)) )
    {
        max3 = fabs(q0_5_p0_5);
    }
    if( (max3 < fabs(q1_0_p0_0)) )
    {
        max3 = fabs(q1_0_p0_0);
    }
    if( (max3 < fabs(q1_1_p0_1)) )
    {
        max3 = fabs(q1_1_p0_1);
    }
    if( (max3 < fabs(q1_2_p0_2)) )
    {
        max3 = fabs(q1_2_p0_2);
    }
    if( (max3 < fabs(q1_3_p0_3)) )
    {
        max3 = fabs(q1_3_p0_3);
    }
    if( (max3 < fabs(q1_4_p0_4)) )
    {
        max3 = fabs(q1_4_p0_4);
    }
    if( (max3 < fabs(q1_5_p0_5)) )
    {
        max3 = fabs(q1_5_p0_5);
    }
    if( (max3 < fabs(q2_0_p0_0)) )
    {
        max3 = fabs(q2_0_p0_0);
    }
    if( (max3 < fabs(q2_1_p0_1)) )
    {
        max3 = fabs(q2_1_p0_1);
    }
    if( (max3 < fabs(q2_2_p0_2)) )
    {
        max3 = fabs(q2_2_p0_2);
    }
    if( (max3 < fabs(q2_3_p0_3)) )
    {
        max3 = fabs(q2_3_p0_3);
    }
    if( (max3 < fabs(q2_4_p0_4)) )
    {
        max3 = fabs(q2_4_p0_4);
    }
    if( (max3 < fabs(q2_5_p0_5)) )
    {
        max3 = fabs(q2_5_p0_5);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 1.49958502193059513986e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.40007476026583916019e-13 * (((max1 * max3) * max2) * max3));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max4 = max1;
    if( (max4 < max2) )
    {
        max4 = max2;
    }
    if( (max4 < max3) )
    {
        max4 = max3;
    }
    double max5 = max1;
    if( (max5 < fabs(p3_1_p0_1)) )
    {
        max5 = fabs(p3_1_p0_1);
    }
    if( (max5 < fabs(p3_2_p0_2)) )
    {
        max5 = fabs(p3_2_p0_2);
    }
    if( (max5 < fabs(p3_0_p0_0)) )
    {
        max5 = fabs(p3_0_p0_0);
    }
    if( (max5 < fabs(p3_3_p0_3)) )
    {
        max5 = fabs(p3_3_p0_3);
    }
    if( (max5 < fabs(p3_4_p0_4)) )
    {
        max5 = fabs(p3_4_p0_4);
    }
    if( (max5 < fabs(p3_5_p0_5)) )
    {
        max5 = fabs(p3_5_p0_5);
    }
    if( (max4 < max5) )
    {
        max4 = max5;
    }
    double max6 = max1;
    if( (max6 < max3) )
    {
        max6 = max3;
    }
    if( (max6 < max5) )
    {
        max6 = max5;
    }
    double max7 = max1;
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    if( (max6 < max7) )
    {
        max6 = max7;
    }
    if( (max4 < max6) )
    {
        max4 = max6;
    }
    if( (max4 < max7) )
    {
        max4 = max7;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    if( (lower_bound_1 < 3.31864264949884013629e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.66564133587113197628e-11 * (((((max5 * max3) * max2) * max7) * max6) * max4));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side3_7d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* q0, const double* q1, const double* q2) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double l1;
    l1 = (1 * (((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double l2;
    l2 = (1 * (((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double p3_6_p0_6 = (p3[6] - p0[6]);
    double l3;
    l3 = (1 * (((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)) + (p3_6_p0_6 * p3_6_p0_6)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double a10;
    a10 = (2 * (((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double a11;
    a11 = (2 * (((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double q2_6_p0_6 = (q2[6] - p0[6]);
    double a12;
    a12 = (2 * (((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)) + (p1_6_p0_6 * q2_6_p0_6)));
    double a20;
    a20 = (2 * (((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)));
    double a21;
    a21 = (2 * (((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)));
    double a22;
    a22 = (2 * (((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)) + (p2_6_p0_6 * q2_6_p0_6)));
    double a30;
    a30 = (2 * (((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)) + (p3_6_p0_6 * q0_6_p0_6)));
    double a31;
    a31 = (2 * (((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)) + (p3_6_p0_6 * q1_6_p0_6)));
    double a32;
    a32 = (2 * (((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)) + (p3_6_p0_6 * q2_6_p0_6)));
    double b00;
    b00 = ((a11 * a22) - (a12 * a21));
    double b01;
    b01 = (a21 - a22);
    double b02;
    b02 = (a12 - a11);
    double b10;
    b10 = ((a12 * a20) - (a10 * a22));
    double b11;
    b11 = (a22 - a20);
    double b12;
    b12 = (a10 - a12);
    double b20;
    b20 = ((a10 * a21) - (a11 * a20));
    double b21;
    b21 = (a20 - a21);
    double b22;
    b22 = (a11 - a10);
    double Delta;
    Delta = ((b00 + b10) + b20);
    double DeltaLambda0;
    DeltaLambda0 = (((b01 * l1) + (b02 * l2)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = (((b11 * l1) + (b12 * l2)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = (((b21 * l1) + (b22 * l2)) + b20);
    double r;
    r = ((Delta * l3) - (((a30 * DeltaLambda0) + (a31 * DeltaLambda1)) + (a32 * DeltaLambda2)));
    double eps;
    double max1 = fabs(p1_1_p0_1);
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    double max2 = fabs(p2_0_p0_0);
    if( (max2 < fabs(p2_1_p0_1)) )
    {
        max2 = fabs(p2_1_p0_1);
    }
    if( (max2 < fabs(p2_2_p0_2)) )
    {
        max2 = fabs(p2_2_p0_2);
    }
    if( (max2 < fabs(p2_3_p0_3)) )
    {
        max2 = fabs(p2_3_p0_3);
    }
    if( (max2 < fabs(p2_4_p0_4)) )
    {
        max2 = fabs(p2_4_p0_4);
    }
    if( (max2 < fabs(p2_5_p0_5)) )
    {
        max2 = fabs(p2_5_p0_5);
    }
    if( (max2 < fabs(p2_6_p0_6)) )
    {
        max2 = fabs(p2_6_p0_6);
    }
    double max3 = fabs(q0_0_p0_0);
    if( (max3 < fabs(q0_1_p0_1)) )
    {
        max3 = fabs(q0_1_p0_1);
    }
    if( (max3 < fabs(q0_2_p0_2)) )
    {
        max3 = fabs(q0_2_p0_2);
    }
    if( (max3 < fabs(q0_3_p0_3)) )
    {
        max3 = fabs(q0_3_p0_3);
    }
    if( (max3 < fabs(q0_4_p0_4)) )
    {
        max3 = fabs(q0_4_p0_4);
    }
    if( (max3 < fabs(q0_5_p0_5)) )
    {
        max3 = fabs(q0_5_p0_5);
    }
    if( (max3 < fabs(q0_6_p0_6)) )
    {
        max3 = fabs(q0_6_p0_6);
    }
    if( (max3 < fabs(q1_0_p0_0)) )
    {
        max3 = fabs(q1_0_p0_0);
    }
    if( (max3 < fabs(q1_1_p0_1)) )
    {
        max3 = fabs(q1_1_p0_1);
    }
    if( (max3 < fabs(q1_2_p0_2)) )
    {
        max3 = fabs(q1_2_p0_2);
    }
    if( (max3 < fabs(q1_3_p0_3)) )
    {
        max3 = fabs(q1_3_p0_3);
    }
    if( (max3 < fabs(q1_4_p0_4)) )
    {
        max3 = fabs(q1_4_p0_4);
    }
    if( (max3 < fabs(q1_5_p0_5)) )
    {
        max3 = fabs(q1_5_p0_5);
    }
    if( (max3 < fabs(q1_6_p0_6)) )
    {
        max3 = fabs(q1_6_p0_6);
    }
    if( (max3 < fabs(q2_0_p0_0)) )
    {
        max3 = fabs(q2_0_p0_0);
    }
    if( (max3 < fabs(q2_1_p0_1)) )
    {
        max3 = fabs(q2_1_p0_1);
    }
    if( (max3 < fabs(q2_2_p0_2)) )
    {
        max3 = fabs(q2_2_p0_2);
    }
    if( (max3 < fabs(q2_3_p0_3)) )
    {
        max3 = fabs(q2_3_p0_3);
    }
    if( (max3 < fabs(q2_4_p0_4)) )
    {
        max3 = fabs(q2_4_p0_4);
    }
    if( (max3 < fabs(q2_5_p0_5)) )
    {
        max3 = fabs(q2_5_p0_5);
    }
    if( (max3 < fabs(q2_6_p0_6)) )
    {
        max3 = fabs(q2_6_p0_6);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max2;
    upper_bound_1 = max2;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 1.36918881183883509035e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (6.33127335329798996022e-13 * (((max1 * max3) * max2) * max3));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max4;
    double max7 = max1;
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    max4 = max7;
    if( (max4 < max2) )
    {
        max4 = max2;
    }
    double max5 = max1;
    if( (max5 < fabs(p3_0_p0_0)) )
    {
        max5 = fabs(p3_0_p0_0);
    }
    if( (max5 < fabs(p3_1_p0_1)) )
    {
        max5 = fabs(p3_1_p0_1);
    }
    if( (max5 < fabs(p3_2_p0_2)) )
    {
        max5 = fabs(p3_2_p0_2);
    }
    if( (max5 < fabs(p3_3_p0_3)) )
    {
        max5 = fabs(p3_3_p0_3);
    }
    if( (max5 < fabs(p3_4_p0_4)) )
    {
        max5 = fabs(p3_4_p0_4);
    }
    if( (max5 < fabs(p3_5_p0_5)) )
    {
        max5 = fabs(p3_5_p0_5);
    }
    if( (max5 < fabs(p3_6_p0_6)) )
    {
        max5 = fabs(p3_6_p0_6);
    }
    if( (max4 < max5) )
    {
        max4 = max5;
    }
    if( (max4 < max1) )
    {
        max4 = max1;
    }
    if( (max4 < max3) )
    {
        max4 = max3;
    }
    double max6 = max7;
    if( (max6 < max5) )
    {
        max6 = max5;
    }
    if( (max6 < max1) )
    {
        max6 = max1;
    }
    if( (max6 < max3) )
    {
        max6 = max3;
    }
    if( (max4 < max6) )
    {
        max4 = max6;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max7;
    upper_bound_1 = max7;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    if( (lower_bound_1 < 3.04548303565602498901e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.78873548804336160566e-11 * (((((max5 * max3) * max2) * max7) * max6) * max4));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


inline int side3_8d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* q0, const double* q1, const double* q2) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double l1;
    l1 = (1 * ((((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)) + (p1_7_p0_7 * p1_7_p0_7)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double p2_7_p0_7 = (p2[7] - p0[7]);
    double l2;
    l2 = (1 * ((((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)) + (p2_7_p0_7 * p2_7_p0_7)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double p3_6_p0_6 = (p3[6] - p0[6]);
    double p3_7_p0_7 = (p3[7] - p0[7]);
    double l3;
    l3 = (1 * ((((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)) + (p3_6_p0_6 * p3_6_p0_6)) + (p3_7_p0_7 * p3_7_p0_7)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    double a10;
    a10 = (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double q1_7_p0_7 = (q1[7] - p0[7]);
    double a11;
    a11 = (2 * ((((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)) + (p1_7_p0_7 * q1_7_p0_7)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double q2_6_p0_6 = (q2[6] - p0[6]);
    double q2_7_p0_7 = (q2[7] - p0[7]);
    double a12;
    a12 = (2 * ((((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)) + (p1_6_p0_6 * q2_6_p0_6)) + (p1_7_p0_7 * q2_7_p0_7)));
    double a20;
    a20 = (2 * ((((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)) + (p2_7_p0_7 * q0_7_p0_7)));
    double a21;
    a21 = (2 * ((((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)) + (p2_7_p0_7 * q1_7_p0_7)));
    double a22;
    a22 = (2 * ((((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)) + (p2_6_p0_6 * q2_6_p0_6)) + (p2_7_p0_7 * q2_7_p0_7)));
    double a30;
    a30 = (2 * ((((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)) + (p3_6_p0_6 * q0_6_p0_6)) + (p3_7_p0_7 * q0_7_p0_7)));
    double a31;
    a31 = (2 * ((((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)) + (p3_6_p0_6 * q1_6_p0_6)) + (p3_7_p0_7 * q1_7_p0_7)));
    double a32;
    a32 = (2 * ((((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)) + (p3_6_p0_6 * q2_6_p0_6)) + (p3_7_p0_7 * q2_7_p0_7)));
    double b00;
    b00 = ((a11 * a22) - (a12 * a21));
    double b01;
    b01 = (a21 - a22);
    double b02;
    b02 = (a12 - a11);
    double b10;
    b10 = ((a12 * a20) - (a10 * a22));
    double b11;
    b11 = (a22 - a20);
    double b12;
    b12 = (a10 - a12);
    double b20;
    b20 = ((a10 * a21) - (a11 * a20));
    double b21;
    b21 = (a20 - a21);
    double b22;
    b22 = (a11 - a10);
    double Delta;
    Delta = ((b00 + b10) + b20);
    double DeltaLambda0;
    DeltaLambda0 = (((b01 * l1) + (b02 * l2)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = (((b11 * l1) + (b12 * l2)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = (((b21 * l1) + (b22 * l2)) + b20);
    double r;
    r = ((Delta * l3) - (((a30 * DeltaLambda0) + (a31 * DeltaLambda1)) + (a32 * DeltaLambda2)));
    double eps;
    double max1 = fabs(p2_1_p0_1);
    if( (max1 < fabs(p2_0_p0_0)) )
    {
        max1 = fabs(p2_0_p0_0);
    }
    if( (max1 < fabs(p2_3_p0_3)) )
    {
        max1 = fabs(p2_3_p0_3);
    }
    if( (max1 < fabs(p2_2_p0_2)) )
    {
        max1 = fabs(p2_2_p0_2);
    }
    if( (max1 < fabs(p2_4_p0_4)) )
    {
        max1 = fabs(p2_4_p0_4);
    }
    if( (max1 < fabs(p2_7_p0_7)) )
    {
        max1 = fabs(p2_7_p0_7);
    }
    if( (max1 < fabs(p2_5_p0_5)) )
    {
        max1 = fabs(p2_5_p0_5);
    }
    if( (max1 < fabs(p2_6_p0_6)) )
    {
        max1 = fabs(p2_6_p0_6);
    }
    double max2 = fabs(p1_4_p0_4);
    if( (max2 < fabs(p1_1_p0_1)) )
    {
        max2 = fabs(p1_1_p0_1);
    }
    if( (max2 < fabs(p1_0_p0_0)) )
    {
        max2 = fabs(p1_0_p0_0);
    }
    if( (max2 < fabs(p1_3_p0_3)) )
    {
        max2 = fabs(p1_3_p0_3);
    }
    if( (max2 < fabs(p1_2_p0_2)) )
    {
        max2 = fabs(p1_2_p0_2);
    }
    if( (max2 < fabs(p1_5_p0_5)) )
    {
        max2 = fabs(p1_5_p0_5);
    }
    if( (max2 < fabs(p1_6_p0_6)) )
    {
        max2 = fabs(p1_6_p0_6);
    }
    if( (max2 < fabs(p1_7_p0_7)) )
    {
        max2 = fabs(p1_7_p0_7);
    }
    double max3 = fabs(q0_0_p0_0);
    if( (max3 < fabs(q0_1_p0_1)) )
    {
        max3 = fabs(q0_1_p0_1);
    }
    if( (max3 < fabs(q0_2_p0_2)) )
    {
        max3 = fabs(q0_2_p0_2);
    }
    if( (max3 < fabs(q0_3_p0_3)) )
    {
        max3 = fabs(q0_3_p0_3);
    }
    if( (max3 < fabs(q0_4_p0_4)) )
    {
        max3 = fabs(q0_4_p0_4);
    }
    if( (max3 < fabs(q0_5_p0_5)) )
    {
        max3 = fabs(q0_5_p0_5);
    }
    if( (max3 < fabs(q0_6_p0_6)) )
    {
        max3 = fabs(q0_6_p0_6);
    }
    if( (max3 < fabs(q0_7_p0_7)) )
    {
        max3 = fabs(q0_7_p0_7);
    }
    if( (max3 < fabs(q1_0_p0_0)) )
    {
        max3 = fabs(q1_0_p0_0);
    }
    if( (max3 < fabs(q1_1_p0_1)) )
    {
        max3 = fabs(q1_1_p0_1);
    }
    if( (max3 < fabs(q1_2_p0_2)) )
    {
        max3 = fabs(q1_2_p0_2);
    }
    if( (max3 < fabs(q1_3_p0_3)) )
    {
        max3 = fabs(q1_3_p0_3);
    }
    if( (max3 < fabs(q1_4_p0_4)) )
    {
        max3 = fabs(q1_4_p0_4);
    }
    if( (max3 < fabs(q1_5_p0_5)) )
    {
        max3 = fabs(q1_5_p0_5);
    }
    if( (max3 < fabs(q1_6_p0_6)) )
    {
        max3 = fabs(q1_6_p0_6);
    }
    if( (max3 < fabs(q1_7_p0_7)) )
    {
        max3 = fabs(q1_7_p0_7);
    }
    if( (max3 < fabs(q2_0_p0_0)) )
    {
        max3 = fabs(q2_0_p0_0);
    }
    if( (max3 < fabs(q2_1_p0_1)) )
    {
        max3 = fabs(q2_1_p0_1);
    }
    if( (max3 < fabs(q2_2_p0_2)) )
    {
        max3 = fabs(q2_2_p0_2);
    }
    if( (max3 < fabs(q2_3_p0_3)) )
    {
        max3 = fabs(q2_3_p0_3);
    }
    if( (max3 < fabs(q2_4_p0_4)) )
    {
        max3 = fabs(q2_4_p0_4);
    }
    if( (max3 < fabs(q2_5_p0_5)) )
    {
        max3 = fabs(q2_5_p0_5);
    }
    if( (max3 < fabs(q2_6_p0_6)) )
    {
        max3 = fabs(q2_6_p0_6);
    }
    if( (max3 < fabs(q2_7_p0_7)) )
    {
        max3 = fabs(q2_7_p0_7);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 1.26419510663115923609e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.71140112255785451890e-13 * (((max2 * max3) * max1) * max3));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max4 = max1;
    if( (max4 < max2) )
    {
        max4 = max2;
    }
    if( (max4 < max3) )
    {
        max4 = max3;
    }
    double max5 = max2;
    if( (max5 < fabs(p3_0_p0_0)) )
    {
        max5 = fabs(p3_0_p0_0);
    }
    if( (max5 < fabs(p3_1_p0_1)) )
    {
        max5 = fabs(p3_1_p0_1);
    }
    if( (max5 < fabs(p3_2_p0_2)) )
    {
        max5 = fabs(p3_2_p0_2);
    }
    if( (max5 < fabs(p3_3_p0_3)) )
    {
        max5 = fabs(p3_3_p0_3);
    }
    if( (max5 < fabs(p3_4_p0_4)) )
    {
        max5 = fabs(p3_4_p0_4);
    }
    if( (max5 < fabs(p3_5_p0_5)) )
    {
        max5 = fabs(p3_5_p0_5);
    }
    if( (max5 < fabs(p3_6_p0_6)) )
    {
        max5 = fabs(p3_6_p0_6);
    }
    if( (max5 < fabs(p3_7_p0_7)) )
    {
        max5 = fabs(p3_7_p0_7);
    }
    if( (max4 < max5) )
    {
        max4 = max5;
    }
    double max6 = max2;
    if( (max6 < max3) )
    {
        max6 = max3;
    }
    if( (max6 < max5) )
    {
        max6 = max5;
    }
    double max7 = max2;
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    if( (max6 < max7) )
    {
        max6 = max7;
    }
    if( (max4 < max6) )
    {
        max4 = max6;
    }
    if( (max4 < max7) )
    {
        max4 = max7;
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    if( (lower_bound_1 < 2.82528483194754087282e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.29807421463370647479e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.37492894694731169807e-11 * (((((max5 * max3) * max1) * max7) * max6) * max4));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


/******* extracted from predicates/side3h.h *******/

inline int side3h_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3, double h0, double h1, double h2, double h3, const double* q0, const double* q1, const double* q2) {
    double l1;
    l1 = (h1 - h0);
    double l2;
    l2 = (h2 - h0);
    double l3;
    l3 = (h3 - h0);
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double a10;
    a10 = (2 * (((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double a11;
    a11 = (2 * (((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double a12;
    a12 = (2 * (((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double a20;
    a20 = (2 * (((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)));
    double a21;
    a21 = (2 * (((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)));
    double a22;
    a22 = (2 * (((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double a30;
    a30 = (2 * (((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)));
    double a31;
    a31 = (2 * (((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)));
    double a32;
    a32 = (2 * (((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)));
    double b00;
    b00 = ((a11 * a22) - (a12 * a21));
    double b01;
    b01 = (a21 - a22);
    double b02;
    b02 = (a12 - a11);
    double b10;
    b10 = ((a12 * a20) - (a10 * a22));
    double b11;
    b11 = (a22 - a20);
    double b12;
    b12 = (a10 - a12);
    double b20;
    b20 = ((a10 * a21) - (a11 * a20));
    double b21;
    b21 = (a20 - a21);
    double b22;
    b22 = (a11 - a10);
    double Delta;
    Delta = ((b00 + b10) + b20);
    double DeltaLambda0;
    DeltaLambda0 = (((b01 * l1) + (b02 * l2)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = (((b11 * l1) + (b12 * l2)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = (((b21 * l1) + (b22 * l2)) + b20);
    double r;
    r = ((Delta * l3) - (((a30 * DeltaLambda0) + (a31 * DeltaLambda1)) + (a32 * DeltaLambda2)));
    double eps;
    double max1 = fabs(q2_2_p0_2);
    if( (max1 < fabs(q0_0_p0_0)) )
    {
        max1 = fabs(q0_0_p0_0);
    }
    if( (max1 < fabs(q0_1_p0_1)) )
    {
        max1 = fabs(q0_1_p0_1);
    }
    if( (max1 < fabs(q0_2_p0_2)) )
    {
        max1 = fabs(q0_2_p0_2);
    }
    if( (max1 < fabs(q1_0_p0_0)) )
    {
        max1 = fabs(q1_0_p0_0);
    }
    if( (max1 < fabs(q1_1_p0_1)) )
    {
        max1 = fabs(q1_1_p0_1);
    }
    if( (max1 < fabs(q1_2_p0_2)) )
    {
        max1 = fabs(q1_2_p0_2);
    }
    if( (max1 < fabs(q2_0_p0_0)) )
    {
        max1 = fabs(q2_0_p0_0);
    }
    if( (max1 < fabs(q2_1_p0_1)) )
    {
        max1 = fabs(q2_1_p0_1);
    }
    double max2 = fabs(p2_0_p0_0);
    if( (max2 < fabs(p2_1_p0_1)) )
    {
        max2 = fabs(p2_1_p0_1);
    }
    if( (max2 < fabs(p2_2_p0_2)) )
    {
        max2 = fabs(p2_2_p0_2);
    }
    double max3 = fabs(p1_0_p0_0);
    if( (max3 < fabs(p1_1_p0_1)) )
    {
        max3 = fabs(p1_1_p0_1);
    }
    if( (max3 < fabs(p1_2_p0_2)) )
    {
        max3 = fabs(p1_2_p0_2);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 2.22985945097100191780e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 2.59614842926741294957e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.99983341597279045654e-14 * (((max3 * max1) * max2) * max1));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    double max4 = max2;
    if( (max4 < fabs(l1)) )
    {
        max4 = fabs(l1);
    }
    if( (max4 < fabs(l2)) )
    {
        max4 = fabs(l2);
    }
    double max5 = max2;
    if( (max5 < max3) )
    {
        max5 = max3;
    }
    if( (max5 < fabs(l3)) )
    {
        max5 = fabs(l3);
    }
    double max6 = max2;
    if( (max6 < fabs(q2_2_p0_2)) )
    {
        max6 = fabs(q2_2_p0_2);
    }
    if( (max6 < fabs(q0_0_p0_0)) )
    {
        max6 = fabs(q0_0_p0_0);
    }
    if( (max6 < fabs(q0_1_p0_1)) )
    {
        max6 = fabs(q0_1_p0_1);
    }
    if( (max6 < fabs(q0_2_p0_2)) )
    {
        max6 = fabs(q0_2_p0_2);
    }
    if( (max6 < fabs(q2_0_p0_0)) )
    {
        max6 = fabs(q2_0_p0_0);
    }
    if( (max6 < fabs(q2_1_p0_1)) )
    {
        max6 = fabs(q2_1_p0_1);
    }
    double max7 = max3;
    if( (max7 < fabs(p3_0_p0_0)) )
    {
        max7 = fabs(p3_0_p0_0);
    }
    if( (max7 < fabs(p3_1_p0_1)) )
    {
        max7 = fabs(p3_1_p0_1);
    }
    if( (max7 < fabs(p3_2_p0_2)) )
    {
        max7 = fabs(p3_2_p0_2);
    }
    int r_sign;
    int int_tmp_result_FFWKCAA;
    lower_bound_1 = max6;
    upper_bound_1 = max6;
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    else
    {
        if( (max7 > upper_bound_1) )
        {
            upper_bound_1 = max7;
        }
    }
    if( (lower_bound_1 < 5.53478725478149652989e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 2.59614842926741294957e+33) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (7.73996217364502738018e-13 * (((((max7 * max1) * max6) * max1) * max5) * max4));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    r_sign = int_tmp_result_FFWKCAA;
    return (Delta_sign * r_sign);
}


/******* extracted from predicates/side3_2dlifted.h *******/

inline int side3_2dlifted_2d_filter( const double* p0, const double* p1, const double* p2, const double* p3, double h0, double h1, double h2, double h3) {
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (h0 - h1);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (h0 - h2);
    double a31;
    a31 = (p3[0] - p0[0]);
    double a32;
    a32 = (p3[1] - p0[1]);
    double a33;
    a33 = (h0 - h3);
    double Delta1;
    Delta1 = ((a21 * a32) - (a22 * a31));
    double Delta2;
    Delta2 = ((a11 * a32) - (a12 * a31));
    double Delta3;
    Delta3 = ((a11 * a22) - (a12 * a21));
    double r;
    r = (((Delta1 * a13) - (Delta2 * a23)) + (Delta3 * a33));
    double eps;
    double max1 = fabs(a11);
    if( (max1 < fabs(a12)) )
    {
        max1 = fabs(a12);
    }
    double max2 = fabs(a21);
    if( (max2 < fabs(a22)) )
    {
        max2 = fabs(a22);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta3_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 5.00368081960964635413e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.88720573725927976811e-16 * (max1 * max2));
        if( (Delta3 > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta3 < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta3_sign = int_tmp_result;
    int int_tmp_result_FFWKCAA;
    double max3 = max1;
    if( (max3 < max2) )
    {
        max3 = max2;
    }
    double max4 = fabs(a13);
    if( (max4 < fabs(a23)) )
    {
        max4 = fabs(a23);
    }
    if( (max4 < fabs(a33)) )
    {
        max4 = fabs(a33);
    }
    double max5 = max2;
    if( (max5 < fabs(a31)) )
    {
        max5 = fabs(a31);
    }
    if( (max5 < fabs(a32)) )
    {
        max5 = fabs(a32);
    }
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (lower_bound_1 < 1.63288018496748314939e-98) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 5.59936185544450928309e+101) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (5.11071278299732992696e-15 * ((max3 * max5) * max4));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta3_sign * int_tmp_result_FFWKCAA);
}

/******* extracted from predicates/side4.h *******/

inline int side4_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4) {
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double a14;
    a14 = -(((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2));
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double a24;
    a24 = -(((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2));
    double a31;
    a31 = (p3[0] - p0[0]);
    double a32;
    a32 = (p3[1] - p0[1]);
    double a33;
    a33 = (p3[2] - p0[2]);
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double a34;
    a34 = -(((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2));
    double a41;
    a41 = (p4[0] - p0[0]);
    double a42;
    a42 = (p4[1] - p0[1]);
    double a43;
    a43 = (p4[2] - p0[2]);
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double a44;
    a44 = -(((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2));
    double Delta1;
    Delta1 = (((a21 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a22 * a43) - (a23 * a42)))) + (a41 * ((a22 * a33) - (a23 * a32))));
    double Delta2;
    Delta2 = (((a11 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a12 * a43) - (a13 * a42)))) + (a41 * ((a12 * a33) - (a13 * a32))));
    double Delta3;
    Delta3 = (((a11 * ((a22 * a43) - (a23 * a42))) - (a21 * ((a12 * a43) - (a13 * a42)))) + (a41 * ((a12 * a23) - (a13 * a22))));
    double Delta4;
    Delta4 = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
    double r;
    r = ((((Delta1 * a14) - (Delta2 * a24)) + (Delta3 * a34)) - (Delta4 * a44));
    double eps;
    double max1 = fabs(a11);
    if( (max1 < fabs(a21)) )
    {
        max1 = fabs(a21);
    }
    if( (max1 < fabs(a31)) )
    {
        max1 = fabs(a31);
    }
    double max2 = fabs(a12);
    if( (max2 < fabs(a13)) )
    {
        max2 = fabs(a13);
    }
    if( (max2 < fabs(a22)) )
    {
        max2 = fabs(a22);
    }
    if( (max2 < fabs(a23)) )
    {
        max2 = fabs(a23);
    }
    double max3 = fabs(a22);
    if( (max3 < fabs(a23)) )
    {
        max3 = fabs(a23);
    }
    if( (max3 < fabs(a32)) )
    {
        max3 = fabs(a32);
    }
    if( (max3 < fabs(a33)) )
    {
        max3 = fabs(a33);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta4_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 1.63288018496748314939e-98) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.21387608851797948065e+60) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (5.11071278299732992696e-15 * ((max2 * max3) * max1));
        if( (Delta4 > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta4 < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta4_sign = int_tmp_result;
    int int_tmp_result_FFWKCAA;
    double max4 = max1;
    if( (max4 < fabs(a41)) )
    {
        max4 = fabs(a41);
    }
    double max5 = max3;
    if( (max5 < max2) )
    {
        max5 = max2;
    }
    double max6 = max3;
    if( (max6 < fabs(a42)) )
    {
        max6 = fabs(a42);
    }
    if( (max6 < fabs(a43)) )
    {
        max6 = fabs(a43);
    }
    double max7 = fabs(p1_0_p0_0);
    if( (max7 < fabs(p1_1_p0_1)) )
    {
        max7 = fabs(p1_1_p0_1);
    }
    if( (max7 < fabs(p1_2_p0_2)) )
    {
        max7 = fabs(p1_2_p0_2);
    }
    if( (max7 < fabs(p2_0_p0_0)) )
    {
        max7 = fabs(p2_0_p0_0);
    }
    if( (max7 < fabs(p2_2_p0_2)) )
    {
        max7 = fabs(p2_2_p0_2);
    }
    if( (max7 < fabs(p2_1_p0_1)) )
    {
        max7 = fabs(p2_1_p0_1);
    }
    if( (max7 < fabs(p3_0_p0_0)) )
    {
        max7 = fabs(p3_0_p0_0);
    }
    if( (max7 < fabs(p3_1_p0_1)) )
    {
        max7 = fabs(p3_1_p0_1);
    }
    if( (max7 < fabs(p3_2_p0_2)) )
    {
        max7 = fabs(p3_2_p0_2);
    }
    if( (max7 < fabs(p4_0_p0_0)) )
    {
        max7 = fabs(p4_0_p0_0);
    }
    if( (max7 < fabs(p4_1_p0_1)) )
    {
        max7 = fabs(p4_1_p0_1);
    }
    if( (max7 < fabs(p4_2_p0_2)) )
    {
        max7 = fabs(p4_2_p0_2);
    }
    lower_bound_1 = max7;
    upper_bound_1 = max7;
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (lower_bound_1 < 1.12285198342304832993e-59) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 3.21387608851797948065e+60) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.24661365310273025710e-13 * ((((max5 * max6) * max4) * max7) * max7));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta4_sign * int_tmp_result_FFWKCAA);
}


inline int side4_4d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, const double* q0, const double* q1, const double* q2, const double* q3) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double l1;
    l1 = (1 * ((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double l2;
    l2 = (1 * ((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double l3;
    l3 = (1 * ((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)));
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double p4_3_p0_3 = (p4[3] - p0[3]);
    double l4;
    l4 = (1 * ((((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2)) + (p4_3_p0_3 * p4_3_p0_3)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double a10;
    a10 = (2 * ((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double a11;
    a11 = (2 * ((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double a12;
    a12 = (2 * ((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)));
    double q3_0_p0_0 = (q3[0] - p0[0]);
    double q3_1_p0_1 = (q3[1] - p0[1]);
    double q3_2_p0_2 = (q3[2] - p0[2]);
    double q3_3_p0_3 = (q3[3] - p0[3]);
    double a13;
    a13 = (2 * ((((p1_0_p0_0 * q3_0_p0_0) + (p1_1_p0_1 * q3_1_p0_1)) + (p1_2_p0_2 * q3_2_p0_2)) + (p1_3_p0_3 * q3_3_p0_3)));
    double a20;
    a20 = (2 * ((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)));
    double a21;
    a21 = (2 * ((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)));
    double a22;
    a22 = (2 * ((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)));
    double a23;
    a23 = (2 * ((((p2_0_p0_0 * q3_0_p0_0) + (p2_1_p0_1 * q3_1_p0_1)) + (p2_2_p0_2 * q3_2_p0_2)) + (p2_3_p0_3 * q3_3_p0_3)));
    double a30;
    a30 = (2 * ((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)));
    double a31;
    a31 = (2 * ((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)));
    double a32;
    a32 = (2 * ((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)));
    double a33;
    a33 = (2 * ((((p3_0_p0_0 * q3_0_p0_0) + (p3_1_p0_1 * q3_1_p0_1)) + (p3_2_p0_2 * q3_2_p0_2)) + (p3_3_p0_3 * q3_3_p0_3)));
    double a40;
    a40 = (2 * ((((p4_0_p0_0 * q0_0_p0_0) + (p4_1_p0_1 * q0_1_p0_1)) + (p4_2_p0_2 * q0_2_p0_2)) + (p4_3_p0_3 * q0_3_p0_3)));
    double a41;
    a41 = (2 * ((((p4_0_p0_0 * q1_0_p0_0) + (p4_1_p0_1 * q1_1_p0_1)) + (p4_2_p0_2 * q1_2_p0_2)) + (p4_3_p0_3 * q1_3_p0_3)));
    double a42;
    a42 = (2 * ((((p4_0_p0_0 * q2_0_p0_0) + (p4_1_p0_1 * q2_1_p0_1)) + (p4_2_p0_2 * q2_2_p0_2)) + (p4_3_p0_3 * q2_3_p0_3)));
    double a43;
    a43 = (2 * ((((p4_0_p0_0 * q3_0_p0_0) + (p4_1_p0_1 * q3_1_p0_1)) + (p4_2_p0_2 * q3_2_p0_2)) + (p4_3_p0_3 * q3_3_p0_3)));
    double b00;
    b00 = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
    double b01;
    b01 = -((((a22 * a33) - (a23 * a32)) + ((a23 * a31) - (a21 * a33))) + ((a21 * a32) - (a22 * a31)));
    double b02;
    b02 = ((((a12 * a33) - (a13 * a32)) + ((a13 * a31) - (a11 * a33))) + ((a11 * a32) - (a12 * a31)));
    double b03;
    b03 = -((((a12 * a23) - (a13 * a22)) + ((a13 * a21) - (a11 * a23))) + ((a11 * a22) - (a12 * a21)));
    double b10;
    b10 = -(((a10 * ((a22 * a33) - (a23 * a32))) - (a20 * ((a12 * a33) - (a13 * a32)))) + (a30 * ((a12 * a23) - (a13 * a22))));
    double b11;
    b11 = ((((a22 * a33) - (a23 * a32)) + ((a23 * a30) - (a20 * a33))) + ((a20 * a32) - (a22 * a30)));
    double b12;
    b12 = -((((a12 * a33) - (a13 * a32)) + ((a13 * a30) - (a10 * a33))) + ((a10 * a32) - (a12 * a30)));
    double b13;
    b13 = ((((a12 * a23) - (a13 * a22)) + ((a13 * a20) - (a10 * a23))) + ((a10 * a22) - (a12 * a20)));
    double b20;
    b20 = (((a10 * ((a21 * a33) - (a23 * a31))) - (a20 * ((a11 * a33) - (a13 * a31)))) + (a30 * ((a11 * a23) - (a13 * a21))));
    double b21;
    b21 = -((((a21 * a33) - (a23 * a31)) + ((a23 * a30) - (a20 * a33))) + ((a20 * a31) - (a21 * a30)));
    double b22;
    b22 = ((((a11 * a33) - (a13 * a31)) + ((a13 * a30) - (a10 * a33))) + ((a10 * a31) - (a11 * a30)));
    double b23;
    b23 = -((((a11 * a23) - (a13 * a21)) + ((a13 * a20) - (a10 * a23))) + ((a10 * a21) - (a11 * a20)));
    double b30;
    b30 = -(((a10 * ((a21 * a32) - (a22 * a31))) - (a20 * ((a11 * a32) - (a12 * a31)))) + (a30 * ((a11 * a22) - (a12 * a21))));
    double b31;
    b31 = ((((a21 * a32) - (a22 * a31)) + ((a22 * a30) - (a20 * a32))) + ((a20 * a31) - (a21 * a30)));
    double b32;
    b32 = -((((a11 * a32) - (a12 * a31)) + ((a12 * a30) - (a10 * a32))) + ((a10 * a31) - (a11 * a30)));
    double b33;
    b33 = ((((a11 * a22) - (a12 * a21)) + ((a12 * a20) - (a10 * a22))) + ((a10 * a21) - (a11 * a20)));
    double Delta;
    Delta = (((b00 + b10) + b20) + b30);
    double DeltaLambda0;
    DeltaLambda0 = ((((b01 * l1) + (b02 * l2)) + (b03 * l3)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = ((((b11 * l1) + (b12 * l2)) + (b13 * l3)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = ((((b21 * l1) + (b22 * l2)) + (b23 * l3)) + b20);
    double DeltaLambda3;
    DeltaLambda3 = ((((b31 * l1) + (b32 * l2)) + (b33 * l3)) + b30);
    double r;
    r = ((Delta * l4) - ((((a40 * DeltaLambda0) + (a41 * DeltaLambda1)) + (a42 * DeltaLambda2)) + (a43 * DeltaLambda3)));
    double eps;
    double max1 = fabs(p1_3_p0_3);
    if( (max1 < fabs(p1_0_p0_0)) )
    {
        max1 = fabs(p1_0_p0_0);
    }
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    double max2 = fabs(p2_2_p0_2);
    if( (max2 < fabs(p2_1_p0_1)) )
    {
        max2 = fabs(p2_1_p0_1);
    }
    if( (max2 < fabs(p2_3_p0_3)) )
    {
        max2 = fabs(p2_3_p0_3);
    }
    if( (max2 < fabs(p2_0_p0_0)) )
    {
        max2 = fabs(p2_0_p0_0);
    }
    double max3 = fabs(p3_0_p0_0);
    if( (max3 < fabs(p3_1_p0_1)) )
    {
        max3 = fabs(p3_1_p0_1);
    }
    if( (max3 < fabs(p3_3_p0_3)) )
    {
        max3 = fabs(p3_3_p0_3);
    }
    if( (max3 < fabs(p3_2_p0_2)) )
    {
        max3 = fabs(p3_2_p0_2);
    }
    double max4 = fabs(q0_3_p0_3);
    if( (max4 < fabs(q0_0_p0_0)) )
    {
        max4 = fabs(q0_0_p0_0);
    }
    if( (max4 < fabs(q0_1_p0_1)) )
    {
        max4 = fabs(q0_1_p0_1);
    }
    if( (max4 < fabs(q0_2_p0_2)) )
    {
        max4 = fabs(q0_2_p0_2);
    }
    if( (max4 < fabs(q1_0_p0_0)) )
    {
        max4 = fabs(q1_0_p0_0);
    }
    if( (max4 < fabs(q1_1_p0_1)) )
    {
        max4 = fabs(q1_1_p0_1);
    }
    if( (max4 < fabs(q1_2_p0_2)) )
    {
        max4 = fabs(q1_2_p0_2);
    }
    if( (max4 < fabs(q1_3_p0_3)) )
    {
        max4 = fabs(q1_3_p0_3);
    }
    double max5 = fabs(q1_0_p0_0);
    if( (max5 < fabs(q1_1_p0_1)) )
    {
        max5 = fabs(q1_1_p0_1);
    }
    if( (max5 < fabs(q1_2_p0_2)) )
    {
        max5 = fabs(q1_2_p0_2);
    }
    if( (max5 < fabs(q1_3_p0_3)) )
    {
        max5 = fabs(q1_3_p0_3);
    }
    if( (max5 < fabs(q2_0_p0_0)) )
    {
        max5 = fabs(q2_0_p0_0);
    }
    if( (max5 < fabs(q2_1_p0_1)) )
    {
        max5 = fabs(q2_1_p0_1);
    }
    if( (max5 < fabs(q2_2_p0_2)) )
    {
        max5 = fabs(q2_2_p0_2);
    }
    if( (max5 < fabs(q2_3_p0_3)) )
    {
        max5 = fabs(q2_3_p0_3);
    }
    double max6 = fabs(q2_0_p0_0);
    if( (max6 < fabs(q2_1_p0_1)) )
    {
        max6 = fabs(q2_1_p0_1);
    }
    if( (max6 < fabs(q2_2_p0_2)) )
    {
        max6 = fabs(q2_2_p0_2);
    }
    if( (max6 < fabs(q2_3_p0_3)) )
    {
        max6 = fabs(q2_3_p0_3);
    }
    if( (max6 < fabs(q3_0_p0_0)) )
    {
        max6 = fabs(q3_0_p0_0);
    }
    if( (max6 < fabs(q3_1_p0_1)) )
    {
        max6 = fabs(q3_1_p0_1);
    }
    if( (max6 < fabs(q3_2_p0_2)) )
    {
        max6 = fabs(q3_2_p0_2);
    }
    if( (max6 < fabs(q3_3_p0_3)) )
    {
        max6 = fabs(q3_3_p0_3);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (lower_bound_1 < 4.14607644401726239868e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.38046888801178809320e-12 * (((((max1 * max4) * max2) * max5) * max3) * max6));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    int int_tmp_result_FFWKCAA;
    double max7 = max1;
    if( (max7 < max2) )
    {
        max7 = max2;
    }
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    if( (max7 < max6) )
    {
        max7 = max6;
    }
    double max8 = max1;
    if( (max8 < fabs(p4_1_p0_1)) )
    {
        max8 = fabs(p4_1_p0_1);
    }
    if( (max8 < fabs(p4_2_p0_2)) )
    {
        max8 = fabs(p4_2_p0_2);
    }
    if( (max8 < fabs(p4_0_p0_0)) )
    {
        max8 = fabs(p4_0_p0_0);
    }
    if( (max8 < fabs(p4_3_p0_3)) )
    {
        max8 = fabs(p4_3_p0_3);
    }
    if( (max7 < max8) )
    {
        max7 = max8;
    }
    double max9 = max1;
    if( (max9 < max5) )
    {
        max9 = max5;
    }
    if( (max9 < max8) )
    {
        max9 = max8;
    }
    double max10;
    double max11 = max4;
    if( (max11 < max5) )
    {
        max11 = max5;
    }
    max10 = max11;
    if( (max10 < max1) )
    {
        max10 = max1;
    }
    if( (max10 < max4) )
    {
        max10 = max4;
    }
    if( (max10 < max5) )
    {
        max10 = max5;
    }
    if( (max10 < max6) )
    {
        max10 = max6;
    }
    lower_bound_1 = max10;
    upper_bound_1 = max10;
    if( (max11 < lower_bound_1) )
    {
        lower_bound_1 = max11;
    }
    else
    {
        if( (max11 > upper_bound_1) )
        {
            upper_bound_1 = max11;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    else
    {
        if( (max7 > upper_bound_1) )
        {
            upper_bound_1 = max7;
        }
    }
    if( (max8 < lower_bound_1) )
    {
        lower_bound_1 = max8;
    }
    if( (max9 < lower_bound_1) )
    {
        lower_bound_1 = max9;
    }
    else
    {
        if( (max9 > upper_bound_1) )
        {
            upper_bound_1 = max9;
        }
    }
    if( (lower_bound_1 < 6.06263132863556750071e-38) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.21914442286604163181e-10 * (((((((max8 * max11) * max2) * max10) * max3) * max10) * max9) * max7));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta_sign * int_tmp_result_FFWKCAA);
}


inline int side4_6d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, const double* q0, const double* q1, const double* q2, const double* q3) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double l1;
    l1 = (1 * ((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double l2;
    l2 = (1 * ((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double l3;
    l3 = (1 * ((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)));
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double p4_3_p0_3 = (p4[3] - p0[3]);
    double p4_4_p0_4 = (p4[4] - p0[4]);
    double p4_5_p0_5 = (p4[5] - p0[5]);
    double l4;
    l4 = (1 * ((((((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2)) + (p4_3_p0_3 * p4_3_p0_3)) + (p4_4_p0_4 * p4_4_p0_4)) + (p4_5_p0_5 * p4_5_p0_5)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double a10;
    a10 = (2 * ((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double a11;
    a11 = (2 * ((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double a12;
    a12 = (2 * ((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)));
    double q3_0_p0_0 = (q3[0] - p0[0]);
    double q3_1_p0_1 = (q3[1] - p0[1]);
    double q3_2_p0_2 = (q3[2] - p0[2]);
    double q3_3_p0_3 = (q3[3] - p0[3]);
    double q3_4_p0_4 = (q3[4] - p0[4]);
    double q3_5_p0_5 = (q3[5] - p0[5]);
    double a13;
    a13 = (2 * ((((((p1_0_p0_0 * q3_0_p0_0) + (p1_1_p0_1 * q3_1_p0_1)) + (p1_2_p0_2 * q3_2_p0_2)) + (p1_3_p0_3 * q3_3_p0_3)) + (p1_4_p0_4 * q3_4_p0_4)) + (p1_5_p0_5 * q3_5_p0_5)));
    double a20;
    a20 = (2 * ((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)));
    double a21;
    a21 = (2 * ((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)));
    double a22;
    a22 = (2 * ((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)));
    double a23;
    a23 = (2 * ((((((p2_0_p0_0 * q3_0_p0_0) + (p2_1_p0_1 * q3_1_p0_1)) + (p2_2_p0_2 * q3_2_p0_2)) + (p2_3_p0_3 * q3_3_p0_3)) + (p2_4_p0_4 * q3_4_p0_4)) + (p2_5_p0_5 * q3_5_p0_5)));
    double a30;
    a30 = (2 * ((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)));
    double a31;
    a31 = (2 * ((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)));
    double a32;
    a32 = (2 * ((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)));
    double a33;
    a33 = (2 * ((((((p3_0_p0_0 * q3_0_p0_0) + (p3_1_p0_1 * q3_1_p0_1)) + (p3_2_p0_2 * q3_2_p0_2)) + (p3_3_p0_3 * q3_3_p0_3)) + (p3_4_p0_4 * q3_4_p0_4)) + (p3_5_p0_5 * q3_5_p0_5)));
    double a40;
    a40 = (2 * ((((((p4_0_p0_0 * q0_0_p0_0) + (p4_1_p0_1 * q0_1_p0_1)) + (p4_2_p0_2 * q0_2_p0_2)) + (p4_3_p0_3 * q0_3_p0_3)) + (p4_4_p0_4 * q0_4_p0_4)) + (p4_5_p0_5 * q0_5_p0_5)));
    double a41;
    a41 = (2 * ((((((p4_0_p0_0 * q1_0_p0_0) + (p4_1_p0_1 * q1_1_p0_1)) + (p4_2_p0_2 * q1_2_p0_2)) + (p4_3_p0_3 * q1_3_p0_3)) + (p4_4_p0_4 * q1_4_p0_4)) + (p4_5_p0_5 * q1_5_p0_5)));
    double a42;
    a42 = (2 * ((((((p4_0_p0_0 * q2_0_p0_0) + (p4_1_p0_1 * q2_1_p0_1)) + (p4_2_p0_2 * q2_2_p0_2)) + (p4_3_p0_3 * q2_3_p0_3)) + (p4_4_p0_4 * q2_4_p0_4)) + (p4_5_p0_5 * q2_5_p0_5)));
    double a43;
    a43 = (2 * ((((((p4_0_p0_0 * q3_0_p0_0) + (p4_1_p0_1 * q3_1_p0_1)) + (p4_2_p0_2 * q3_2_p0_2)) + (p4_3_p0_3 * q3_3_p0_3)) + (p4_4_p0_4 * q3_4_p0_4)) + (p4_5_p0_5 * q3_5_p0_5)));
    double b00;
    b00 = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
    double b01;
    b01 = -((((a22 * a33) - (a23 * a32)) + ((a23 * a31) - (a21 * a33))) + ((a21 * a32) - (a22 * a31)));
    double b02;
    b02 = ((((a12 * a33) - (a13 * a32)) + ((a13 * a31) - (a11 * a33))) + ((a11 * a32) - (a12 * a31)));
    double b03;
    b03 = -((((a12 * a23) - (a13 * a22)) + ((a13 * a21) - (a11 * a23))) + ((a11 * a22) - (a12 * a21)));
    double b10;
    b10 = -(((a10 * ((a22 * a33) - (a23 * a32))) - (a20 * ((a12 * a33) - (a13 * a32)))) + (a30 * ((a12 * a23) - (a13 * a22))));
    double b11;
    b11 = ((((a22 * a33) - (a23 * a32)) + ((a23 * a30) - (a20 * a33))) + ((a20 * a32) - (a22 * a30)));
    double b12;
    b12 = -((((a12 * a33) - (a13 * a32)) + ((a13 * a30) - (a10 * a33))) + ((a10 * a32) - (a12 * a30)));
    double b13;
    b13 = ((((a12 * a23) - (a13 * a22)) + ((a13 * a20) - (a10 * a23))) + ((a10 * a22) - (a12 * a20)));
    double b20;
    b20 = (((a10 * ((a21 * a33) - (a23 * a31))) - (a20 * ((a11 * a33) - (a13 * a31)))) + (a30 * ((a11 * a23) - (a13 * a21))));
    double b21;
    b21 = -((((a21 * a33) - (a23 * a31)) + ((a23 * a30) - (a20 * a33))) + ((a20 * a31) - (a21 * a30)));
    double b22;
    b22 = ((((a11 * a33) - (a13 * a31)) + ((a13 * a30) - (a10 * a33))) + ((a10 * a31) - (a11 * a30)));
    double b23;
    b23 = -((((a11 * a23) - (a13 * a21)) + ((a13 * a20) - (a10 * a23))) + ((a10 * a21) - (a11 * a20)));
    double b30;
    b30 = -(((a10 * ((a21 * a32) - (a22 * a31))) - (a20 * ((a11 * a32) - (a12 * a31)))) + (a30 * ((a11 * a22) - (a12 * a21))));
    double b31;
    b31 = ((((a21 * a32) - (a22 * a31)) + ((a22 * a30) - (a20 * a32))) + ((a20 * a31) - (a21 * a30)));
    double b32;
    b32 = -((((a11 * a32) - (a12 * a31)) + ((a12 * a30) - (a10 * a32))) + ((a10 * a31) - (a11 * a30)));
    double b33;
    b33 = ((((a11 * a22) - (a12 * a21)) + ((a12 * a20) - (a10 * a22))) + ((a10 * a21) - (a11 * a20)));
    double Delta;
    Delta = (((b00 + b10) + b20) + b30);
    double DeltaLambda0;
    DeltaLambda0 = ((((b01 * l1) + (b02 * l2)) + (b03 * l3)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = ((((b11 * l1) + (b12 * l2)) + (b13 * l3)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = ((((b21 * l1) + (b22 * l2)) + (b23 * l3)) + b20);
    double DeltaLambda3;
    DeltaLambda3 = ((((b31 * l1) + (b32 * l2)) + (b33 * l3)) + b30);
    double r;
    r = ((Delta * l4) - ((((a40 * DeltaLambda0) + (a41 * DeltaLambda1)) + (a42 * DeltaLambda2)) + (a43 * DeltaLambda3)));
    double eps;
    double max1 = fabs(p3_2_p0_2);
    if( (max1 < fabs(p3_0_p0_0)) )
    {
        max1 = fabs(p3_0_p0_0);
    }
    if( (max1 < fabs(p3_3_p0_3)) )
    {
        max1 = fabs(p3_3_p0_3);
    }
    if( (max1 < fabs(p3_4_p0_4)) )
    {
        max1 = fabs(p3_4_p0_4);
    }
    if( (max1 < fabs(p3_1_p0_1)) )
    {
        max1 = fabs(p3_1_p0_1);
    }
    if( (max1 < fabs(p3_5_p0_5)) )
    {
        max1 = fabs(p3_5_p0_5);
    }
    double max2 = fabs(p2_1_p0_1);
    if( (max2 < fabs(p2_4_p0_4)) )
    {
        max2 = fabs(p2_4_p0_4);
    }
    if( (max2 < fabs(p2_2_p0_2)) )
    {
        max2 = fabs(p2_2_p0_2);
    }
    if( (max2 < fabs(p2_0_p0_0)) )
    {
        max2 = fabs(p2_0_p0_0);
    }
    if( (max2 < fabs(p2_3_p0_3)) )
    {
        max2 = fabs(p2_3_p0_3);
    }
    if( (max2 < fabs(p2_5_p0_5)) )
    {
        max2 = fabs(p2_5_p0_5);
    }
    double max3 = fabs(p1_0_p0_0);
    if( (max3 < fabs(p1_1_p0_1)) )
    {
        max3 = fabs(p1_1_p0_1);
    }
    if( (max3 < fabs(p1_2_p0_2)) )
    {
        max3 = fabs(p1_2_p0_2);
    }
    if( (max3 < fabs(p1_3_p0_3)) )
    {
        max3 = fabs(p1_3_p0_3);
    }
    if( (max3 < fabs(p1_4_p0_4)) )
    {
        max3 = fabs(p1_4_p0_4);
    }
    if( (max3 < fabs(p1_5_p0_5)) )
    {
        max3 = fabs(p1_5_p0_5);
    }
    double max4 = fabs(q0_0_p0_0);
    if( (max4 < fabs(q0_1_p0_1)) )
    {
        max4 = fabs(q0_1_p0_1);
    }
    if( (max4 < fabs(q0_2_p0_2)) )
    {
        max4 = fabs(q0_2_p0_2);
    }
    if( (max4 < fabs(q0_3_p0_3)) )
    {
        max4 = fabs(q0_3_p0_3);
    }
    if( (max4 < fabs(q0_4_p0_4)) )
    {
        max4 = fabs(q0_4_p0_4);
    }
    if( (max4 < fabs(q0_5_p0_5)) )
    {
        max4 = fabs(q0_5_p0_5);
    }
    if( (max4 < fabs(q1_0_p0_0)) )
    {
        max4 = fabs(q1_0_p0_0);
    }
    if( (max4 < fabs(q1_1_p0_1)) )
    {
        max4 = fabs(q1_1_p0_1);
    }
    if( (max4 < fabs(q1_2_p0_2)) )
    {
        max4 = fabs(q1_2_p0_2);
    }
    if( (max4 < fabs(q1_3_p0_3)) )
    {
        max4 = fabs(q1_3_p0_3);
    }
    if( (max4 < fabs(q1_4_p0_4)) )
    {
        max4 = fabs(q1_4_p0_4);
    }
    if( (max4 < fabs(q1_5_p0_5)) )
    {
        max4 = fabs(q1_5_p0_5);
    }
    double max5 = fabs(q1_0_p0_0);
    if( (max5 < fabs(q1_1_p0_1)) )
    {
        max5 = fabs(q1_1_p0_1);
    }
    if( (max5 < fabs(q1_2_p0_2)) )
    {
        max5 = fabs(q1_2_p0_2);
    }
    if( (max5 < fabs(q1_3_p0_3)) )
    {
        max5 = fabs(q1_3_p0_3);
    }
    if( (max5 < fabs(q1_4_p0_4)) )
    {
        max5 = fabs(q1_4_p0_4);
    }
    if( (max5 < fabs(q1_5_p0_5)) )
    {
        max5 = fabs(q1_5_p0_5);
    }
    if( (max5 < fabs(q2_0_p0_0)) )
    {
        max5 = fabs(q2_0_p0_0);
    }
    if( (max5 < fabs(q2_1_p0_1)) )
    {
        max5 = fabs(q2_1_p0_1);
    }
    if( (max5 < fabs(q2_2_p0_2)) )
    {
        max5 = fabs(q2_2_p0_2);
    }
    if( (max5 < fabs(q2_3_p0_3)) )
    {
        max5 = fabs(q2_3_p0_3);
    }
    if( (max5 < fabs(q2_4_p0_4)) )
    {
        max5 = fabs(q2_4_p0_4);
    }
    if( (max5 < fabs(q2_5_p0_5)) )
    {
        max5 = fabs(q2_5_p0_5);
    }
    double max6 = fabs(q2_0_p0_0);
    if( (max6 < fabs(q2_1_p0_1)) )
    {
        max6 = fabs(q2_1_p0_1);
    }
    if( (max6 < fabs(q2_2_p0_2)) )
    {
        max6 = fabs(q2_2_p0_2);
    }
    if( (max6 < fabs(q2_3_p0_3)) )
    {
        max6 = fabs(q2_3_p0_3);
    }
    if( (max6 < fabs(q2_4_p0_4)) )
    {
        max6 = fabs(q2_4_p0_4);
    }
    if( (max6 < fabs(q2_5_p0_5)) )
    {
        max6 = fabs(q2_5_p0_5);
    }
    if( (max6 < fabs(q3_0_p0_0)) )
    {
        max6 = fabs(q3_0_p0_0);
    }
    if( (max6 < fabs(q3_1_p0_1)) )
    {
        max6 = fabs(q3_1_p0_1);
    }
    if( (max6 < fabs(q3_2_p0_2)) )
    {
        max6 = fabs(q3_2_p0_2);
    }
    if( (max6 < fabs(q3_3_p0_3)) )
    {
        max6 = fabs(q3_3_p0_3);
    }
    if( (max6 < fabs(q3_4_p0_4)) )
    {
        max6 = fabs(q3_4_p0_4);
    }
    if( (max6 < fabs(q3_5_p0_5)) )
    {
        max6 = fabs(q3_5_p0_5);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max4;
    upper_bound_1 = max4;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (lower_bound_1 < 3.31864264949884013629e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.66564133587113165316e-11 * (((((max3 * max4) * max2) * max5) * max1) * max6));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    int int_tmp_result_FFWKCAA;
    double max7 = max1;
    if( (max7 < max2) )
    {
        max7 = max2;
    }
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    double max8 = max3;
    if( (max8 < fabs(p4_5_p0_5)) )
    {
        max8 = fabs(p4_5_p0_5);
    }
    if( (max8 < fabs(p4_0_p0_0)) )
    {
        max8 = fabs(p4_0_p0_0);
    }
    if( (max8 < fabs(p4_1_p0_1)) )
    {
        max8 = fabs(p4_1_p0_1);
    }
    if( (max8 < fabs(p4_2_p0_2)) )
    {
        max8 = fabs(p4_2_p0_2);
    }
    if( (max8 < fabs(p4_3_p0_3)) )
    {
        max8 = fabs(p4_3_p0_3);
    }
    if( (max8 < fabs(p4_4_p0_4)) )
    {
        max8 = fabs(p4_4_p0_4);
    }
    if( (max7 < max8) )
    {
        max7 = max8;
    }
    if( (max7 < max6) )
    {
        max7 = max6;
    }
    double max9 = max3;
    if( (max9 < max8) )
    {
        max9 = max8;
    }
    if( (max9 < max5) )
    {
        max9 = max5;
    }
    double max10 = max4;
    if( (max10 < max3) )
    {
        max10 = max3;
    }
    if( (max10 < max5) )
    {
        max10 = max5;
    }
    double max11 = max4;
    if( (max11 < max5) )
    {
        max11 = max5;
    }
    if( (max10 < max11) )
    {
        max10 = max11;
    }
    if( (max10 < max6) )
    {
        max10 = max6;
    }
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    else
    {
        if( (max7 > upper_bound_1) )
        {
            upper_bound_1 = max7;
        }
    }
    if( (max8 < lower_bound_1) )
    {
        lower_bound_1 = max8;
    }
    if( (max9 < lower_bound_1) )
    {
        lower_bound_1 = max9;
    }
    else
    {
        if( (max9 > upper_bound_1) )
        {
            upper_bound_1 = max9;
        }
    }
    if( (max10 < lower_bound_1) )
    {
        lower_bound_1 = max10;
    }
    else
    {
        if( (max10 > upper_bound_1) )
        {
            upper_bound_1 = max10;
        }
    }
    if( (max11 < lower_bound_1) )
    {
        lower_bound_1 = max11;
    }
    else
    {
        if( (max11 > upper_bound_1) )
        {
            upper_bound_1 = max11;
        }
    }
    if( (lower_bound_1 < 4.87975611107819181771e-38) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (6.92085098542795335117e-10 * (((((((max8 * max11) * max2) * max10) * max1) * max10) * max9) * max7));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta_sign * int_tmp_result_FFWKCAA);
}


inline int side4_7d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, const double* q0, const double* q1, const double* q2, const double* q3) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double l1;
    l1 = (1 * (((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double l2;
    l2 = (1 * (((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double p3_6_p0_6 = (p3[6] - p0[6]);
    double l3;
    l3 = (1 * (((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)) + (p3_6_p0_6 * p3_6_p0_6)));
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double p4_3_p0_3 = (p4[3] - p0[3]);
    double p4_4_p0_4 = (p4[4] - p0[4]);
    double p4_5_p0_5 = (p4[5] - p0[5]);
    double p4_6_p0_6 = (p4[6] - p0[6]);
    double l4;
    l4 = (1 * (((((((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2)) + (p4_3_p0_3 * p4_3_p0_3)) + (p4_4_p0_4 * p4_4_p0_4)) + (p4_5_p0_5 * p4_5_p0_5)) + (p4_6_p0_6 * p4_6_p0_6)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double a10;
    a10 = (2 * (((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double a11;
    a11 = (2 * (((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double q2_6_p0_6 = (q2[6] - p0[6]);
    double a12;
    a12 = (2 * (((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)) + (p1_6_p0_6 * q2_6_p0_6)));
    double q3_0_p0_0 = (q3[0] - p0[0]);
    double q3_1_p0_1 = (q3[1] - p0[1]);
    double q3_2_p0_2 = (q3[2] - p0[2]);
    double q3_3_p0_3 = (q3[3] - p0[3]);
    double q3_4_p0_4 = (q3[4] - p0[4]);
    double q3_5_p0_5 = (q3[5] - p0[5]);
    double q3_6_p0_6 = (q3[6] - p0[6]);
    double a13;
    a13 = (2 * (((((((p1_0_p0_0 * q3_0_p0_0) + (p1_1_p0_1 * q3_1_p0_1)) + (p1_2_p0_2 * q3_2_p0_2)) + (p1_3_p0_3 * q3_3_p0_3)) + (p1_4_p0_4 * q3_4_p0_4)) + (p1_5_p0_5 * q3_5_p0_5)) + (p1_6_p0_6 * q3_6_p0_6)));
    double a20;
    a20 = (2 * (((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)));
    double a21;
    a21 = (2 * (((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)));
    double a22;
    a22 = (2 * (((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)) + (p2_6_p0_6 * q2_6_p0_6)));
    double a23;
    a23 = (2 * (((((((p2_0_p0_0 * q3_0_p0_0) + (p2_1_p0_1 * q3_1_p0_1)) + (p2_2_p0_2 * q3_2_p0_2)) + (p2_3_p0_3 * q3_3_p0_3)) + (p2_4_p0_4 * q3_4_p0_4)) + (p2_5_p0_5 * q3_5_p0_5)) + (p2_6_p0_6 * q3_6_p0_6)));
    double a30;
    a30 = (2 * (((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)) + (p3_6_p0_6 * q0_6_p0_6)));
    double a31;
    a31 = (2 * (((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)) + (p3_6_p0_6 * q1_6_p0_6)));
    double a32;
    a32 = (2 * (((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)) + (p3_6_p0_6 * q2_6_p0_6)));
    double a33;
    a33 = (2 * (((((((p3_0_p0_0 * q3_0_p0_0) + (p3_1_p0_1 * q3_1_p0_1)) + (p3_2_p0_2 * q3_2_p0_2)) + (p3_3_p0_3 * q3_3_p0_3)) + (p3_4_p0_4 * q3_4_p0_4)) + (p3_5_p0_5 * q3_5_p0_5)) + (p3_6_p0_6 * q3_6_p0_6)));
    double a40;
    a40 = (2 * (((((((p4_0_p0_0 * q0_0_p0_0) + (p4_1_p0_1 * q0_1_p0_1)) + (p4_2_p0_2 * q0_2_p0_2)) + (p4_3_p0_3 * q0_3_p0_3)) + (p4_4_p0_4 * q0_4_p0_4)) + (p4_5_p0_5 * q0_5_p0_5)) + (p4_6_p0_6 * q0_6_p0_6)));
    double a41;
    a41 = (2 * (((((((p4_0_p0_0 * q1_0_p0_0) + (p4_1_p0_1 * q1_1_p0_1)) + (p4_2_p0_2 * q1_2_p0_2)) + (p4_3_p0_3 * q1_3_p0_3)) + (p4_4_p0_4 * q1_4_p0_4)) + (p4_5_p0_5 * q1_5_p0_5)) + (p4_6_p0_6 * q1_6_p0_6)));
    double a42;
    a42 = (2 * (((((((p4_0_p0_0 * q2_0_p0_0) + (p4_1_p0_1 * q2_1_p0_1)) + (p4_2_p0_2 * q2_2_p0_2)) + (p4_3_p0_3 * q2_3_p0_3)) + (p4_4_p0_4 * q2_4_p0_4)) + (p4_5_p0_5 * q2_5_p0_5)) + (p4_6_p0_6 * q2_6_p0_6)));
    double a43;
    a43 = (2 * (((((((p4_0_p0_0 * q3_0_p0_0) + (p4_1_p0_1 * q3_1_p0_1)) + (p4_2_p0_2 * q3_2_p0_2)) + (p4_3_p0_3 * q3_3_p0_3)) + (p4_4_p0_4 * q3_4_p0_4)) + (p4_5_p0_5 * q3_5_p0_5)) + (p4_6_p0_6 * q3_6_p0_6)));
    double b00;
    b00 = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
    double b01;
    b01 = -((((a22 * a33) - (a23 * a32)) + ((a23 * a31) - (a21 * a33))) + ((a21 * a32) - (a22 * a31)));
    double b02;
    b02 = ((((a12 * a33) - (a13 * a32)) + ((a13 * a31) - (a11 * a33))) + ((a11 * a32) - (a12 * a31)));
    double b03;
    b03 = -((((a12 * a23) - (a13 * a22)) + ((a13 * a21) - (a11 * a23))) + ((a11 * a22) - (a12 * a21)));
    double b10;
    b10 = -(((a10 * ((a22 * a33) - (a23 * a32))) - (a20 * ((a12 * a33) - (a13 * a32)))) + (a30 * ((a12 * a23) - (a13 * a22))));
    double b11;
    b11 = ((((a22 * a33) - (a23 * a32)) + ((a23 * a30) - (a20 * a33))) + ((a20 * a32) - (a22 * a30)));
    double b12;
    b12 = -((((a12 * a33) - (a13 * a32)) + ((a13 * a30) - (a10 * a33))) + ((a10 * a32) - (a12 * a30)));
    double b13;
    b13 = ((((a12 * a23) - (a13 * a22)) + ((a13 * a20) - (a10 * a23))) + ((a10 * a22) - (a12 * a20)));
    double b20;
    b20 = (((a10 * ((a21 * a33) - (a23 * a31))) - (a20 * ((a11 * a33) - (a13 * a31)))) + (a30 * ((a11 * a23) - (a13 * a21))));
    double b21;
    b21 = -((((a21 * a33) - (a23 * a31)) + ((a23 * a30) - (a20 * a33))) + ((a20 * a31) - (a21 * a30)));
    double b22;
    b22 = ((((a11 * a33) - (a13 * a31)) + ((a13 * a30) - (a10 * a33))) + ((a10 * a31) - (a11 * a30)));
    double b23;
    b23 = -((((a11 * a23) - (a13 * a21)) + ((a13 * a20) - (a10 * a23))) + ((a10 * a21) - (a11 * a20)));
    double b30;
    b30 = -(((a10 * ((a21 * a32) - (a22 * a31))) - (a20 * ((a11 * a32) - (a12 * a31)))) + (a30 * ((a11 * a22) - (a12 * a21))));
    double b31;
    b31 = ((((a21 * a32) - (a22 * a31)) + ((a22 * a30) - (a20 * a32))) + ((a20 * a31) - (a21 * a30)));
    double b32;
    b32 = -((((a11 * a32) - (a12 * a31)) + ((a12 * a30) - (a10 * a32))) + ((a10 * a31) - (a11 * a30)));
    double b33;
    b33 = ((((a11 * a22) - (a12 * a21)) + ((a12 * a20) - (a10 * a22))) + ((a10 * a21) - (a11 * a20)));
    double Delta;
    Delta = (((b00 + b10) + b20) + b30);
    double DeltaLambda0;
    DeltaLambda0 = ((((b01 * l1) + (b02 * l2)) + (b03 * l3)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = ((((b11 * l1) + (b12 * l2)) + (b13 * l3)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = ((((b21 * l1) + (b22 * l2)) + (b23 * l3)) + b20);
    double DeltaLambda3;
    DeltaLambda3 = ((((b31 * l1) + (b32 * l2)) + (b33 * l3)) + b30);
    double r;
    r = ((Delta * l4) - ((((a40 * DeltaLambda0) + (a41 * DeltaLambda1)) + (a42 * DeltaLambda2)) + (a43 * DeltaLambda3)));
    double eps;
    double max1 = fabs(p1_0_p0_0);
    if( (max1 < fabs(p1_1_p0_1)) )
    {
        max1 = fabs(p1_1_p0_1);
    }
    if( (max1 < fabs(p1_2_p0_2)) )
    {
        max1 = fabs(p1_2_p0_2);
    }
    if( (max1 < fabs(p1_3_p0_3)) )
    {
        max1 = fabs(p1_3_p0_3);
    }
    if( (max1 < fabs(p1_4_p0_4)) )
    {
        max1 = fabs(p1_4_p0_4);
    }
    if( (max1 < fabs(p1_5_p0_5)) )
    {
        max1 = fabs(p1_5_p0_5);
    }
    if( (max1 < fabs(p1_6_p0_6)) )
    {
        max1 = fabs(p1_6_p0_6);
    }
    double max2 = fabs(p3_0_p0_0);
    if( (max2 < fabs(p3_4_p0_4)) )
    {
        max2 = fabs(p3_4_p0_4);
    }
    if( (max2 < fabs(p3_2_p0_2)) )
    {
        max2 = fabs(p3_2_p0_2);
    }
    if( (max2 < fabs(p3_1_p0_1)) )
    {
        max2 = fabs(p3_1_p0_1);
    }
    if( (max2 < fabs(p3_3_p0_3)) )
    {
        max2 = fabs(p3_3_p0_3);
    }
    if( (max2 < fabs(p3_5_p0_5)) )
    {
        max2 = fabs(p3_5_p0_5);
    }
    if( (max2 < fabs(p3_6_p0_6)) )
    {
        max2 = fabs(p3_6_p0_6);
    }
    double max3 = fabs(p2_5_p0_5);
    if( (max3 < fabs(p2_2_p0_2)) )
    {
        max3 = fabs(p2_2_p0_2);
    }
    if( (max3 < fabs(p2_3_p0_3)) )
    {
        max3 = fabs(p2_3_p0_3);
    }
    if( (max3 < fabs(p2_0_p0_0)) )
    {
        max3 = fabs(p2_0_p0_0);
    }
    if( (max3 < fabs(p2_1_p0_1)) )
    {
        max3 = fabs(p2_1_p0_1);
    }
    if( (max3 < fabs(p2_6_p0_6)) )
    {
        max3 = fabs(p2_6_p0_6);
    }
    if( (max3 < fabs(p2_4_p0_4)) )
    {
        max3 = fabs(p2_4_p0_4);
    }
    double max4 = fabs(q0_0_p0_0);
    if( (max4 < fabs(q0_1_p0_1)) )
    {
        max4 = fabs(q0_1_p0_1);
    }
    if( (max4 < fabs(q0_2_p0_2)) )
    {
        max4 = fabs(q0_2_p0_2);
    }
    if( (max4 < fabs(q0_3_p0_3)) )
    {
        max4 = fabs(q0_3_p0_3);
    }
    if( (max4 < fabs(q0_4_p0_4)) )
    {
        max4 = fabs(q0_4_p0_4);
    }
    if( (max4 < fabs(q0_5_p0_5)) )
    {
        max4 = fabs(q0_5_p0_5);
    }
    if( (max4 < fabs(q0_6_p0_6)) )
    {
        max4 = fabs(q0_6_p0_6);
    }
    if( (max4 < fabs(q1_0_p0_0)) )
    {
        max4 = fabs(q1_0_p0_0);
    }
    if( (max4 < fabs(q1_1_p0_1)) )
    {
        max4 = fabs(q1_1_p0_1);
    }
    if( (max4 < fabs(q1_2_p0_2)) )
    {
        max4 = fabs(q1_2_p0_2);
    }
    if( (max4 < fabs(q1_3_p0_3)) )
    {
        max4 = fabs(q1_3_p0_3);
    }
    if( (max4 < fabs(q1_4_p0_4)) )
    {
        max4 = fabs(q1_4_p0_4);
    }
    if( (max4 < fabs(q1_5_p0_5)) )
    {
        max4 = fabs(q1_5_p0_5);
    }
    if( (max4 < fabs(q1_6_p0_6)) )
    {
        max4 = fabs(q1_6_p0_6);
    }
    double max5 = fabs(q1_0_p0_0);
    if( (max5 < fabs(q1_1_p0_1)) )
    {
        max5 = fabs(q1_1_p0_1);
    }
    if( (max5 < fabs(q1_2_p0_2)) )
    {
        max5 = fabs(q1_2_p0_2);
    }
    if( (max5 < fabs(q1_3_p0_3)) )
    {
        max5 = fabs(q1_3_p0_3);
    }
    if( (max5 < fabs(q1_4_p0_4)) )
    {
        max5 = fabs(q1_4_p0_4);
    }
    if( (max5 < fabs(q1_5_p0_5)) )
    {
        max5 = fabs(q1_5_p0_5);
    }
    if( (max5 < fabs(q1_6_p0_6)) )
    {
        max5 = fabs(q1_6_p0_6);
    }
    if( (max5 < fabs(q2_0_p0_0)) )
    {
        max5 = fabs(q2_0_p0_0);
    }
    if( (max5 < fabs(q2_1_p0_1)) )
    {
        max5 = fabs(q2_1_p0_1);
    }
    if( (max5 < fabs(q2_2_p0_2)) )
    {
        max5 = fabs(q2_2_p0_2);
    }
    if( (max5 < fabs(q2_3_p0_3)) )
    {
        max5 = fabs(q2_3_p0_3);
    }
    if( (max5 < fabs(q2_4_p0_4)) )
    {
        max5 = fabs(q2_4_p0_4);
    }
    if( (max5 < fabs(q2_5_p0_5)) )
    {
        max5 = fabs(q2_5_p0_5);
    }
    if( (max5 < fabs(q2_6_p0_6)) )
    {
        max5 = fabs(q2_6_p0_6);
    }
    double max6 = fabs(q2_0_p0_0);
    if( (max6 < fabs(q2_1_p0_1)) )
    {
        max6 = fabs(q2_1_p0_1);
    }
    if( (max6 < fabs(q2_2_p0_2)) )
    {
        max6 = fabs(q2_2_p0_2);
    }
    if( (max6 < fabs(q2_3_p0_3)) )
    {
        max6 = fabs(q2_3_p0_3);
    }
    if( (max6 < fabs(q2_4_p0_4)) )
    {
        max6 = fabs(q2_4_p0_4);
    }
    if( (max6 < fabs(q2_5_p0_5)) )
    {
        max6 = fabs(q2_5_p0_5);
    }
    if( (max6 < fabs(q2_6_p0_6)) )
    {
        max6 = fabs(q2_6_p0_6);
    }
    if( (max6 < fabs(q3_0_p0_0)) )
    {
        max6 = fabs(q3_0_p0_0);
    }
    if( (max6 < fabs(q3_1_p0_1)) )
    {
        max6 = fabs(q3_1_p0_1);
    }
    if( (max6 < fabs(q3_2_p0_2)) )
    {
        max6 = fabs(q3_2_p0_2);
    }
    if( (max6 < fabs(q3_3_p0_3)) )
    {
        max6 = fabs(q3_3_p0_3);
    }
    if( (max6 < fabs(q3_4_p0_4)) )
    {
        max6 = fabs(q3_4_p0_4);
    }
    if( (max6 < fabs(q3_5_p0_5)) )
    {
        max6 = fabs(q3_5_p0_5);
    }
    if( (max6 < fabs(q3_6_p0_6)) )
    {
        max6 = fabs(q3_6_p0_6);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (lower_bound_1 < 3.04548303565602498901e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.78873548804336160566e-11 * (((((max1 * max4) * max3) * max5) * max2) * max6));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    int int_tmp_result_FFWKCAA;
    double max7 = max3;
    if( (max7 < max1) )
    {
        max7 = max1;
    }
    if( (max7 < max2) )
    {
        max7 = max2;
    }
    if( (max7 < max6) )
    {
        max7 = max6;
    }
    double max8 = max1;
    if( (max8 < fabs(p4_4_p0_4)) )
    {
        max8 = fabs(p4_4_p0_4);
    }
    if( (max8 < fabs(p4_1_p0_1)) )
    {
        max8 = fabs(p4_1_p0_1);
    }
    if( (max8 < fabs(p4_5_p0_5)) )
    {
        max8 = fabs(p4_5_p0_5);
    }
    if( (max8 < fabs(p4_0_p0_0)) )
    {
        max8 = fabs(p4_0_p0_0);
    }
    if( (max8 < fabs(p4_2_p0_2)) )
    {
        max8 = fabs(p4_2_p0_2);
    }
    if( (max8 < fabs(p4_3_p0_3)) )
    {
        max8 = fabs(p4_3_p0_3);
    }
    if( (max8 < fabs(p4_6_p0_6)) )
    {
        max8 = fabs(p4_6_p0_6);
    }
    if( (max7 < max8) )
    {
        max7 = max8;
    }
    double max9 = max1;
    if( (max9 < max5) )
    {
        max9 = max5;
    }
    if( (max9 < max8) )
    {
        max9 = max8;
    }
    double max10 = max4;
    if( (max10 < max1) )
    {
        max10 = max1;
    }
    if( (max10 < max5) )
    {
        max10 = max5;
    }
    if( (max10 < max6) )
    {
        max10 = max6;
    }
    double max11 = max4;
    if( (max11 < max5) )
    {
        max11 = max5;
    }
    if( (max10 < max11) )
    {
        max10 = max11;
    }
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    else
    {
        if( (max7 > upper_bound_1) )
        {
            upper_bound_1 = max7;
        }
    }
    if( (max8 < lower_bound_1) )
    {
        lower_bound_1 = max8;
    }
    if( (max9 < lower_bound_1) )
    {
        lower_bound_1 = max9;
    }
    else
    {
        if( (max9 > upper_bound_1) )
        {
            upper_bound_1 = max9;
        }
    }
    if( (max10 < lower_bound_1) )
    {
        lower_bound_1 = max10;
    }
    else
    {
        if( (max10 > upper_bound_1) )
        {
            upper_bound_1 = max10;
        }
    }
    if( (max11 < lower_bound_1) )
    {
        lower_bound_1 = max11;
    }
    else
    {
        if( (max11 > upper_bound_1) )
        {
            upper_bound_1 = max11;
        }
    }
    if( (lower_bound_1 < 4.48906690519700369396e-38) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.34926049830188433875e-09 * (((((((max8 * max11) * max3) * max10) * max2) * max10) * max9) * max7));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta_sign * int_tmp_result_FFWKCAA);
}


inline int side4_8d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, const double* q0, const double* q1, const double* q2, const double* q3) {
    double p1_0_p0_0 = (p1[0] - p0[0]);
    double p1_1_p0_1 = (p1[1] - p0[1]);
    double p1_2_p0_2 = (p1[2] - p0[2]);
    double p1_3_p0_3 = (p1[3] - p0[3]);
    double p1_4_p0_4 = (p1[4] - p0[4]);
    double p1_5_p0_5 = (p1[5] - p0[5]);
    double p1_6_p0_6 = (p1[6] - p0[6]);
    double p1_7_p0_7 = (p1[7] - p0[7]);
    double l1;
    l1 = (1 * ((((((((p1_0_p0_0 * p1_0_p0_0) + (p1_1_p0_1 * p1_1_p0_1)) + (p1_2_p0_2 * p1_2_p0_2)) + (p1_3_p0_3 * p1_3_p0_3)) + (p1_4_p0_4 * p1_4_p0_4)) + (p1_5_p0_5 * p1_5_p0_5)) + (p1_6_p0_6 * p1_6_p0_6)) + (p1_7_p0_7 * p1_7_p0_7)));
    double p2_0_p0_0 = (p2[0] - p0[0]);
    double p2_1_p0_1 = (p2[1] - p0[1]);
    double p2_2_p0_2 = (p2[2] - p0[2]);
    double p2_3_p0_3 = (p2[3] - p0[3]);
    double p2_4_p0_4 = (p2[4] - p0[4]);
    double p2_5_p0_5 = (p2[5] - p0[5]);
    double p2_6_p0_6 = (p2[6] - p0[6]);
    double p2_7_p0_7 = (p2[7] - p0[7]);
    double l2;
    l2 = (1 * ((((((((p2_0_p0_0 * p2_0_p0_0) + (p2_1_p0_1 * p2_1_p0_1)) + (p2_2_p0_2 * p2_2_p0_2)) + (p2_3_p0_3 * p2_3_p0_3)) + (p2_4_p0_4 * p2_4_p0_4)) + (p2_5_p0_5 * p2_5_p0_5)) + (p2_6_p0_6 * p2_6_p0_6)) + (p2_7_p0_7 * p2_7_p0_7)));
    double p3_0_p0_0 = (p3[0] - p0[0]);
    double p3_1_p0_1 = (p3[1] - p0[1]);
    double p3_2_p0_2 = (p3[2] - p0[2]);
    double p3_3_p0_3 = (p3[3] - p0[3]);
    double p3_4_p0_4 = (p3[4] - p0[4]);
    double p3_5_p0_5 = (p3[5] - p0[5]);
    double p3_6_p0_6 = (p3[6] - p0[6]);
    double p3_7_p0_7 = (p3[7] - p0[7]);
    double l3;
    l3 = (1 * ((((((((p3_0_p0_0 * p3_0_p0_0) + (p3_1_p0_1 * p3_1_p0_1)) + (p3_2_p0_2 * p3_2_p0_2)) + (p3_3_p0_3 * p3_3_p0_3)) + (p3_4_p0_4 * p3_4_p0_4)) + (p3_5_p0_5 * p3_5_p0_5)) + (p3_6_p0_6 * p3_6_p0_6)) + (p3_7_p0_7 * p3_7_p0_7)));
    double p4_0_p0_0 = (p4[0] - p0[0]);
    double p4_1_p0_1 = (p4[1] - p0[1]);
    double p4_2_p0_2 = (p4[2] - p0[2]);
    double p4_3_p0_3 = (p4[3] - p0[3]);
    double p4_4_p0_4 = (p4[4] - p0[4]);
    double p4_5_p0_5 = (p4[5] - p0[5]);
    double p4_6_p0_6 = (p4[6] - p0[6]);
    double p4_7_p0_7 = (p4[7] - p0[7]);
    double l4;
    l4 = (1 * ((((((((p4_0_p0_0 * p4_0_p0_0) + (p4_1_p0_1 * p4_1_p0_1)) + (p4_2_p0_2 * p4_2_p0_2)) + (p4_3_p0_3 * p4_3_p0_3)) + (p4_4_p0_4 * p4_4_p0_4)) + (p4_5_p0_5 * p4_5_p0_5)) + (p4_6_p0_6 * p4_6_p0_6)) + (p4_7_p0_7 * p4_7_p0_7)));
    double q0_0_p0_0 = (q0[0] - p0[0]);
    double q0_1_p0_1 = (q0[1] - p0[1]);
    double q0_2_p0_2 = (q0[2] - p0[2]);
    double q0_3_p0_3 = (q0[3] - p0[3]);
    double q0_4_p0_4 = (q0[4] - p0[4]);
    double q0_5_p0_5 = (q0[5] - p0[5]);
    double q0_6_p0_6 = (q0[6] - p0[6]);
    double q0_7_p0_7 = (q0[7] - p0[7]);
    double a10;
    a10 = (2 * ((((((((p1_0_p0_0 * q0_0_p0_0) + (p1_1_p0_1 * q0_1_p0_1)) + (p1_2_p0_2 * q0_2_p0_2)) + (p1_3_p0_3 * q0_3_p0_3)) + (p1_4_p0_4 * q0_4_p0_4)) + (p1_5_p0_5 * q0_5_p0_5)) + (p1_6_p0_6 * q0_6_p0_6)) + (p1_7_p0_7 * q0_7_p0_7)));
    double q1_0_p0_0 = (q1[0] - p0[0]);
    double q1_1_p0_1 = (q1[1] - p0[1]);
    double q1_2_p0_2 = (q1[2] - p0[2]);
    double q1_3_p0_3 = (q1[3] - p0[3]);
    double q1_4_p0_4 = (q1[4] - p0[4]);
    double q1_5_p0_5 = (q1[5] - p0[5]);
    double q1_6_p0_6 = (q1[6] - p0[6]);
    double q1_7_p0_7 = (q1[7] - p0[7]);
    double a11;
    a11 = (2 * ((((((((p1_0_p0_0 * q1_0_p0_0) + (p1_1_p0_1 * q1_1_p0_1)) + (p1_2_p0_2 * q1_2_p0_2)) + (p1_3_p0_3 * q1_3_p0_3)) + (p1_4_p0_4 * q1_4_p0_4)) + (p1_5_p0_5 * q1_5_p0_5)) + (p1_6_p0_6 * q1_6_p0_6)) + (p1_7_p0_7 * q1_7_p0_7)));
    double q2_0_p0_0 = (q2[0] - p0[0]);
    double q2_1_p0_1 = (q2[1] - p0[1]);
    double q2_2_p0_2 = (q2[2] - p0[2]);
    double q2_3_p0_3 = (q2[3] - p0[3]);
    double q2_4_p0_4 = (q2[4] - p0[4]);
    double q2_5_p0_5 = (q2[5] - p0[5]);
    double q2_6_p0_6 = (q2[6] - p0[6]);
    double q2_7_p0_7 = (q2[7] - p0[7]);
    double a12;
    a12 = (2 * ((((((((p1_0_p0_0 * q2_0_p0_0) + (p1_1_p0_1 * q2_1_p0_1)) + (p1_2_p0_2 * q2_2_p0_2)) + (p1_3_p0_3 * q2_3_p0_3)) + (p1_4_p0_4 * q2_4_p0_4)) + (p1_5_p0_5 * q2_5_p0_5)) + (p1_6_p0_6 * q2_6_p0_6)) + (p1_7_p0_7 * q2_7_p0_7)));
    double q3_0_p0_0 = (q3[0] - p0[0]);
    double q3_1_p0_1 = (q3[1] - p0[1]);
    double q3_2_p0_2 = (q3[2] - p0[2]);
    double q3_3_p0_3 = (q3[3] - p0[3]);
    double q3_4_p0_4 = (q3[4] - p0[4]);
    double q3_5_p0_5 = (q3[5] - p0[5]);
    double q3_6_p0_6 = (q3[6] - p0[6]);
    double q3_7_p0_7 = (q3[7] - p0[7]);
    double a13;
    a13 = (2 * ((((((((p1_0_p0_0 * q3_0_p0_0) + (p1_1_p0_1 * q3_1_p0_1)) + (p1_2_p0_2 * q3_2_p0_2)) + (p1_3_p0_3 * q3_3_p0_3)) + (p1_4_p0_4 * q3_4_p0_4)) + (p1_5_p0_5 * q3_5_p0_5)) + (p1_6_p0_6 * q3_6_p0_6)) + (p1_7_p0_7 * q3_7_p0_7)));
    double a20;
    a20 = (2 * ((((((((p2_0_p0_0 * q0_0_p0_0) + (p2_1_p0_1 * q0_1_p0_1)) + (p2_2_p0_2 * q0_2_p0_2)) + (p2_3_p0_3 * q0_3_p0_3)) + (p2_4_p0_4 * q0_4_p0_4)) + (p2_5_p0_5 * q0_5_p0_5)) + (p2_6_p0_6 * q0_6_p0_6)) + (p2_7_p0_7 * q0_7_p0_7)));
    double a21;
    a21 = (2 * ((((((((p2_0_p0_0 * q1_0_p0_0) + (p2_1_p0_1 * q1_1_p0_1)) + (p2_2_p0_2 * q1_2_p0_2)) + (p2_3_p0_3 * q1_3_p0_3)) + (p2_4_p0_4 * q1_4_p0_4)) + (p2_5_p0_5 * q1_5_p0_5)) + (p2_6_p0_6 * q1_6_p0_6)) + (p2_7_p0_7 * q1_7_p0_7)));
    double a22;
    a22 = (2 * ((((((((p2_0_p0_0 * q2_0_p0_0) + (p2_1_p0_1 * q2_1_p0_1)) + (p2_2_p0_2 * q2_2_p0_2)) + (p2_3_p0_3 * q2_3_p0_3)) + (p2_4_p0_4 * q2_4_p0_4)) + (p2_5_p0_5 * q2_5_p0_5)) + (p2_6_p0_6 * q2_6_p0_6)) + (p2_7_p0_7 * q2_7_p0_7)));
    double a23;
    a23 = (2 * ((((((((p2_0_p0_0 * q3_0_p0_0) + (p2_1_p0_1 * q3_1_p0_1)) + (p2_2_p0_2 * q3_2_p0_2)) + (p2_3_p0_3 * q3_3_p0_3)) + (p2_4_p0_4 * q3_4_p0_4)) + (p2_5_p0_5 * q3_5_p0_5)) + (p2_6_p0_6 * q3_6_p0_6)) + (p2_7_p0_7 * q3_7_p0_7)));
    double a30;
    a30 = (2 * ((((((((p3_0_p0_0 * q0_0_p0_0) + (p3_1_p0_1 * q0_1_p0_1)) + (p3_2_p0_2 * q0_2_p0_2)) + (p3_3_p0_3 * q0_3_p0_3)) + (p3_4_p0_4 * q0_4_p0_4)) + (p3_5_p0_5 * q0_5_p0_5)) + (p3_6_p0_6 * q0_6_p0_6)) + (p3_7_p0_7 * q0_7_p0_7)));
    double a31;
    a31 = (2 * ((((((((p3_0_p0_0 * q1_0_p0_0) + (p3_1_p0_1 * q1_1_p0_1)) + (p3_2_p0_2 * q1_2_p0_2)) + (p3_3_p0_3 * q1_3_p0_3)) + (p3_4_p0_4 * q1_4_p0_4)) + (p3_5_p0_5 * q1_5_p0_5)) + (p3_6_p0_6 * q1_6_p0_6)) + (p3_7_p0_7 * q1_7_p0_7)));
    double a32;
    a32 = (2 * ((((((((p3_0_p0_0 * q2_0_p0_0) + (p3_1_p0_1 * q2_1_p0_1)) + (p3_2_p0_2 * q2_2_p0_2)) + (p3_3_p0_3 * q2_3_p0_3)) + (p3_4_p0_4 * q2_4_p0_4)) + (p3_5_p0_5 * q2_5_p0_5)) + (p3_6_p0_6 * q2_6_p0_6)) + (p3_7_p0_7 * q2_7_p0_7)));
    double a33;
    a33 = (2 * ((((((((p3_0_p0_0 * q3_0_p0_0) + (p3_1_p0_1 * q3_1_p0_1)) + (p3_2_p0_2 * q3_2_p0_2)) + (p3_3_p0_3 * q3_3_p0_3)) + (p3_4_p0_4 * q3_4_p0_4)) + (p3_5_p0_5 * q3_5_p0_5)) + (p3_6_p0_6 * q3_6_p0_6)) + (p3_7_p0_7 * q3_7_p0_7)));
    double a40;
    a40 = (2 * ((((((((p4_0_p0_0 * q0_0_p0_0) + (p4_1_p0_1 * q0_1_p0_1)) + (p4_2_p0_2 * q0_2_p0_2)) + (p4_3_p0_3 * q0_3_p0_3)) + (p4_4_p0_4 * q0_4_p0_4)) + (p4_5_p0_5 * q0_5_p0_5)) + (p4_6_p0_6 * q0_6_p0_6)) + (p4_7_p0_7 * q0_7_p0_7)));
    double a41;
    a41 = (2 * ((((((((p4_0_p0_0 * q1_0_p0_0) + (p4_1_p0_1 * q1_1_p0_1)) + (p4_2_p0_2 * q1_2_p0_2)) + (p4_3_p0_3 * q1_3_p0_3)) + (p4_4_p0_4 * q1_4_p0_4)) + (p4_5_p0_5 * q1_5_p0_5)) + (p4_6_p0_6 * q1_6_p0_6)) + (p4_7_p0_7 * q1_7_p0_7)));
    double a42;
    a42 = (2 * ((((((((p4_0_p0_0 * q2_0_p0_0) + (p4_1_p0_1 * q2_1_p0_1)) + (p4_2_p0_2 * q2_2_p0_2)) + (p4_3_p0_3 * q2_3_p0_3)) + (p4_4_p0_4 * q2_4_p0_4)) + (p4_5_p0_5 * q2_5_p0_5)) + (p4_6_p0_6 * q2_6_p0_6)) + (p4_7_p0_7 * q2_7_p0_7)));
    double a43;
    a43 = (2 * ((((((((p4_0_p0_0 * q3_0_p0_0) + (p4_1_p0_1 * q3_1_p0_1)) + (p4_2_p0_2 * q3_2_p0_2)) + (p4_3_p0_3 * q3_3_p0_3)) + (p4_4_p0_4 * q3_4_p0_4)) + (p4_5_p0_5 * q3_5_p0_5)) + (p4_6_p0_6 * q3_6_p0_6)) + (p4_7_p0_7 * q3_7_p0_7)));
    double b00;
    b00 = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
    double b01;
    b01 = -((((a22 * a33) - (a23 * a32)) + ((a23 * a31) - (a21 * a33))) + ((a21 * a32) - (a22 * a31)));
    double b02;
    b02 = ((((a12 * a33) - (a13 * a32)) + ((a13 * a31) - (a11 * a33))) + ((a11 * a32) - (a12 * a31)));
    double b03;
    b03 = -((((a12 * a23) - (a13 * a22)) + ((a13 * a21) - (a11 * a23))) + ((a11 * a22) - (a12 * a21)));
    double b10;
    b10 = -(((a10 * ((a22 * a33) - (a23 * a32))) - (a20 * ((a12 * a33) - (a13 * a32)))) + (a30 * ((a12 * a23) - (a13 * a22))));
    double b11;
    b11 = ((((a22 * a33) - (a23 * a32)) + ((a23 * a30) - (a20 * a33))) + ((a20 * a32) - (a22 * a30)));
    double b12;
    b12 = -((((a12 * a33) - (a13 * a32)) + ((a13 * a30) - (a10 * a33))) + ((a10 * a32) - (a12 * a30)));
    double b13;
    b13 = ((((a12 * a23) - (a13 * a22)) + ((a13 * a20) - (a10 * a23))) + ((a10 * a22) - (a12 * a20)));
    double b20;
    b20 = (((a10 * ((a21 * a33) - (a23 * a31))) - (a20 * ((a11 * a33) - (a13 * a31)))) + (a30 * ((a11 * a23) - (a13 * a21))));
    double b21;
    b21 = -((((a21 * a33) - (a23 * a31)) + ((a23 * a30) - (a20 * a33))) + ((a20 * a31) - (a21 * a30)));
    double b22;
    b22 = ((((a11 * a33) - (a13 * a31)) + ((a13 * a30) - (a10 * a33))) + ((a10 * a31) - (a11 * a30)));
    double b23;
    b23 = -((((a11 * a23) - (a13 * a21)) + ((a13 * a20) - (a10 * a23))) + ((a10 * a21) - (a11 * a20)));
    double b30;
    b30 = -(((a10 * ((a21 * a32) - (a22 * a31))) - (a20 * ((a11 * a32) - (a12 * a31)))) + (a30 * ((a11 * a22) - (a12 * a21))));
    double b31;
    b31 = ((((a21 * a32) - (a22 * a31)) + ((a22 * a30) - (a20 * a32))) + ((a20 * a31) - (a21 * a30)));
    double b32;
    b32 = -((((a11 * a32) - (a12 * a31)) + ((a12 * a30) - (a10 * a32))) + ((a10 * a31) - (a11 * a30)));
    double b33;
    b33 = ((((a11 * a22) - (a12 * a21)) + ((a12 * a20) - (a10 * a22))) + ((a10 * a21) - (a11 * a20)));
    double Delta;
    Delta = (((b00 + b10) + b20) + b30);
    double DeltaLambda0;
    DeltaLambda0 = ((((b01 * l1) + (b02 * l2)) + (b03 * l3)) + b00);
    double DeltaLambda1;
    DeltaLambda1 = ((((b11 * l1) + (b12 * l2)) + (b13 * l3)) + b10);
    double DeltaLambda2;
    DeltaLambda2 = ((((b21 * l1) + (b22 * l2)) + (b23 * l3)) + b20);
    double DeltaLambda3;
    DeltaLambda3 = ((((b31 * l1) + (b32 * l2)) + (b33 * l3)) + b30);
    double r;
    r = ((Delta * l4) - ((((a40 * DeltaLambda0) + (a41 * DeltaLambda1)) + (a42 * DeltaLambda2)) + (a43 * DeltaLambda3)));
    double eps;
    double max1 = fabs(p2_5_p0_5);
    if( (max1 < fabs(p2_3_p0_3)) )
    {
        max1 = fabs(p2_3_p0_3);
    }
    if( (max1 < fabs(p2_0_p0_0)) )
    {
        max1 = fabs(p2_0_p0_0);
    }
    if( (max1 < fabs(p2_1_p0_1)) )
    {
        max1 = fabs(p2_1_p0_1);
    }
    if( (max1 < fabs(p2_6_p0_6)) )
    {
        max1 = fabs(p2_6_p0_6);
    }
    if( (max1 < fabs(p2_2_p0_2)) )
    {
        max1 = fabs(p2_2_p0_2);
    }
    if( (max1 < fabs(p2_4_p0_4)) )
    {
        max1 = fabs(p2_4_p0_4);
    }
    if( (max1 < fabs(p2_7_p0_7)) )
    {
        max1 = fabs(p2_7_p0_7);
    }
    double max2 = fabs(p1_4_p0_4);
    if( (max2 < fabs(p1_3_p0_3)) )
    {
        max2 = fabs(p1_3_p0_3);
    }
    if( (max2 < fabs(p1_7_p0_7)) )
    {
        max2 = fabs(p1_7_p0_7);
    }
    if( (max2 < fabs(p1_0_p0_0)) )
    {
        max2 = fabs(p1_0_p0_0);
    }
    if( (max2 < fabs(p1_2_p0_2)) )
    {
        max2 = fabs(p1_2_p0_2);
    }
    if( (max2 < fabs(p1_5_p0_5)) )
    {
        max2 = fabs(p1_5_p0_5);
    }
    if( (max2 < fabs(p1_1_p0_1)) )
    {
        max2 = fabs(p1_1_p0_1);
    }
    if( (max2 < fabs(p1_6_p0_6)) )
    {
        max2 = fabs(p1_6_p0_6);
    }
    double max3 = fabs(p3_3_p0_3);
    if( (max3 < fabs(p3_0_p0_0)) )
    {
        max3 = fabs(p3_0_p0_0);
    }
    if( (max3 < fabs(p3_1_p0_1)) )
    {
        max3 = fabs(p3_1_p0_1);
    }
    if( (max3 < fabs(p3_2_p0_2)) )
    {
        max3 = fabs(p3_2_p0_2);
    }
    if( (max3 < fabs(p3_4_p0_4)) )
    {
        max3 = fabs(p3_4_p0_4);
    }
    if( (max3 < fabs(p3_5_p0_5)) )
    {
        max3 = fabs(p3_5_p0_5);
    }
    if( (max3 < fabs(p3_6_p0_6)) )
    {
        max3 = fabs(p3_6_p0_6);
    }
    if( (max3 < fabs(p3_7_p0_7)) )
    {
        max3 = fabs(p3_7_p0_7);
    }
    double max4 = fabs(q0_0_p0_0);
    if( (max4 < fabs(q0_1_p0_1)) )
    {
        max4 = fabs(q0_1_p0_1);
    }
    if( (max4 < fabs(q0_2_p0_2)) )
    {
        max4 = fabs(q0_2_p0_2);
    }
    if( (max4 < fabs(q0_3_p0_3)) )
    {
        max4 = fabs(q0_3_p0_3);
    }
    if( (max4 < fabs(q0_4_p0_4)) )
    {
        max4 = fabs(q0_4_p0_4);
    }
    if( (max4 < fabs(q0_5_p0_5)) )
    {
        max4 = fabs(q0_5_p0_5);
    }
    if( (max4 < fabs(q0_6_p0_6)) )
    {
        max4 = fabs(q0_6_p0_6);
    }
    if( (max4 < fabs(q0_7_p0_7)) )
    {
        max4 = fabs(q0_7_p0_7);
    }
    if( (max4 < fabs(q1_0_p0_0)) )
    {
        max4 = fabs(q1_0_p0_0);
    }
    if( (max4 < fabs(q1_1_p0_1)) )
    {
        max4 = fabs(q1_1_p0_1);
    }
    if( (max4 < fabs(q1_2_p0_2)) )
    {
        max4 = fabs(q1_2_p0_2);
    }
    if( (max4 < fabs(q1_3_p0_3)) )
    {
        max4 = fabs(q1_3_p0_3);
    }
    if( (max4 < fabs(q1_4_p0_4)) )
    {
        max4 = fabs(q1_4_p0_4);
    }
    if( (max4 < fabs(q1_5_p0_5)) )
    {
        max4 = fabs(q1_5_p0_5);
    }
    if( (max4 < fabs(q1_6_p0_6)) )
    {
        max4 = fabs(q1_6_p0_6);
    }
    if( (max4 < fabs(q1_7_p0_7)) )
    {
        max4 = fabs(q1_7_p0_7);
    }
    double max5 = fabs(q1_0_p0_0);
    if( (max5 < fabs(q1_1_p0_1)) )
    {
        max5 = fabs(q1_1_p0_1);
    }
    if( (max5 < fabs(q1_2_p0_2)) )
    {
        max5 = fabs(q1_2_p0_2);
    }
    if( (max5 < fabs(q1_3_p0_3)) )
    {
        max5 = fabs(q1_3_p0_3);
    }
    if( (max5 < fabs(q1_4_p0_4)) )
    {
        max5 = fabs(q1_4_p0_4);
    }
    if( (max5 < fabs(q1_5_p0_5)) )
    {
        max5 = fabs(q1_5_p0_5);
    }
    if( (max5 < fabs(q1_6_p0_6)) )
    {
        max5 = fabs(q1_6_p0_6);
    }
    if( (max5 < fabs(q1_7_p0_7)) )
    {
        max5 = fabs(q1_7_p0_7);
    }
    if( (max5 < fabs(q2_0_p0_0)) )
    {
        max5 = fabs(q2_0_p0_0);
    }
    if( (max5 < fabs(q2_1_p0_1)) )
    {
        max5 = fabs(q2_1_p0_1);
    }
    if( (max5 < fabs(q2_2_p0_2)) )
    {
        max5 = fabs(q2_2_p0_2);
    }
    if( (max5 < fabs(q2_3_p0_3)) )
    {
        max5 = fabs(q2_3_p0_3);
    }
    if( (max5 < fabs(q2_4_p0_4)) )
    {
        max5 = fabs(q2_4_p0_4);
    }
    if( (max5 < fabs(q2_5_p0_5)) )
    {
        max5 = fabs(q2_5_p0_5);
    }
    if( (max5 < fabs(q2_6_p0_6)) )
    {
        max5 = fabs(q2_6_p0_6);
    }
    if( (max5 < fabs(q2_7_p0_7)) )
    {
        max5 = fabs(q2_7_p0_7);
    }
    double max6 = fabs(q2_0_p0_0);
    if( (max6 < fabs(q2_1_p0_1)) )
    {
        max6 = fabs(q2_1_p0_1);
    }
    if( (max6 < fabs(q2_2_p0_2)) )
    {
        max6 = fabs(q2_2_p0_2);
    }
    if( (max6 < fabs(q2_3_p0_3)) )
    {
        max6 = fabs(q2_3_p0_3);
    }
    if( (max6 < fabs(q2_4_p0_4)) )
    {
        max6 = fabs(q2_4_p0_4);
    }
    if( (max6 < fabs(q2_5_p0_5)) )
    {
        max6 = fabs(q2_5_p0_5);
    }
    if( (max6 < fabs(q2_6_p0_6)) )
    {
        max6 = fabs(q2_6_p0_6);
    }
    if( (max6 < fabs(q2_7_p0_7)) )
    {
        max6 = fabs(q2_7_p0_7);
    }
    if( (max6 < fabs(q3_0_p0_0)) )
    {
        max6 = fabs(q3_0_p0_0);
    }
    if( (max6 < fabs(q3_1_p0_1)) )
    {
        max6 = fabs(q3_1_p0_1);
    }
    if( (max6 < fabs(q3_2_p0_2)) )
    {
        max6 = fabs(q3_2_p0_2);
    }
    if( (max6 < fabs(q3_3_p0_3)) )
    {
        max6 = fabs(q3_3_p0_3);
    }
    if( (max6 < fabs(q3_4_p0_4)) )
    {
        max6 = fabs(q3_4_p0_4);
    }
    if( (max6 < fabs(q3_5_p0_5)) )
    {
        max6 = fabs(q3_5_p0_5);
    }
    if( (max6 < fabs(q3_6_p0_6)) )
    {
        max6 = fabs(q3_6_p0_6);
    }
    if( (max6 < fabs(q3_7_p0_7)) )
    {
        max6 = fabs(q3_7_p0_7);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta_sign;
    int int_tmp_result;
    lower_bound_1 = max4;
    upper_bound_1 = max4;
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (lower_bound_1 < 2.82528483194754087282e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (4.37492894694731040560e-11 * (((((max2 * max4) * max1) * max5) * max3) * max6));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta_sign = int_tmp_result;
    int int_tmp_result_FFWKCAA;
    double max7 = max1;
    if( (max7 < max3) )
    {
        max7 = max3;
    }
    double max8 = max2;
    if( (max8 < fabs(p4_3_p0_3)) )
    {
        max8 = fabs(p4_3_p0_3);
    }
    if( (max8 < fabs(p4_1_p0_1)) )
    {
        max8 = fabs(p4_1_p0_1);
    }
    if( (max8 < fabs(p4_4_p0_4)) )
    {
        max8 = fabs(p4_4_p0_4);
    }
    if( (max8 < fabs(p4_0_p0_0)) )
    {
        max8 = fabs(p4_0_p0_0);
    }
    if( (max8 < fabs(p4_2_p0_2)) )
    {
        max8 = fabs(p4_2_p0_2);
    }
    if( (max8 < fabs(p4_5_p0_5)) )
    {
        max8 = fabs(p4_5_p0_5);
    }
    if( (max8 < fabs(p4_6_p0_6)) )
    {
        max8 = fabs(p4_6_p0_6);
    }
    if( (max8 < fabs(p4_7_p0_7)) )
    {
        max8 = fabs(p4_7_p0_7);
    }
    if( (max7 < max8) )
    {
        max7 = max8;
    }
    if( (max7 < max2) )
    {
        max7 = max2;
    }
    if( (max7 < max6) )
    {
        max7 = max6;
    }
    double max9 = max8;
    if( (max9 < max2) )
    {
        max9 = max2;
    }
    if( (max9 < max5) )
    {
        max9 = max5;
    }
    double max10 = max4;
    double max11 = max4;
    if( (max11 < max5) )
    {
        max11 = max5;
    }
    if( (max10 < max11) )
    {
        max10 = max11;
    }
    if( (max10 < max2) )
    {
        max10 = max2;
    }
    if( (max10 < max6) )
    {
        max10 = max6;
    }
    if( (max10 < max5) )
    {
        max10 = max5;
    }
    lower_bound_1 = max11;
    upper_bound_1 = max11;
    if( (max9 < lower_bound_1) )
    {
        lower_bound_1 = max9;
    }
    else
    {
        if( (max9 > upper_bound_1) )
        {
            upper_bound_1 = max9;
        }
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    else
    {
        if( (max7 > upper_bound_1) )
        {
            upper_bound_1 = max7;
        }
    }
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    if( (max8 < lower_bound_1) )
    {
        lower_bound_1 = max8;
    }
    if( (max10 < lower_bound_1) )
    {
        lower_bound_1 = max10;
    }
    else
    {
        if( (max10 > upper_bound_1) )
        {
            upper_bound_1 = max10;
        }
    }
    if( (lower_bound_1 < 4.17402518597284772324e-38) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 4.83570327845851562508e+24) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.41492645607254025015e-09 * (((((((max8 * max11) * max1) * max10) * max3) * max10) * max9) * max7));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta_sign * int_tmp_result_FFWKCAA);
}


/******* extracted from predicates/side4h.h *******/

inline int side4h_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, double h0, double h1, double h2, double h3, double h4) {
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double a14;
    a14 = (h0 - h1);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double a24;
    a24 = (h0 - h2);
    double a31;
    a31 = (p3[0] - p0[0]);
    double a32;
    a32 = (p3[1] - p0[1]);
    double a33;
    a33 = (p3[2] - p0[2]);
    double a34;
    a34 = (h0 - h3);
    double a41;
    a41 = (p4[0] - p0[0]);
    double a42;
    a42 = (p4[1] - p0[1]);
    double a43;
    a43 = (p4[2] - p0[2]);
    double a44;
    a44 = (h0 - h4);
    double Delta1;
    Delta1 = (((a21 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a22 * a43) - (a23 * a42)))) + (a41 * ((a22 * a33) - (a23 * a32))));
    double Delta2;
    Delta2 = (((a11 * ((a32 * a43) - (a33 * a42))) - (a31 * ((a12 * a43) - (a13 * a42)))) + (a41 * ((a12 * a33) - (a13 * a32))));
    double Delta3;
    Delta3 = (((a11 * ((a22 * a43) - (a23 * a42))) - (a21 * ((a12 * a43) - (a13 * a42)))) + (a41 * ((a12 * a23) - (a13 * a22))));
    double Delta4;
    Delta4 = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
    double r;
    r = ((((Delta1 * a14) - (Delta2 * a24)) + (Delta3 * a34)) - (Delta4 * a44));
    double eps;
    double max1 = fabs(a11);
    if( (max1 < fabs(a21)) )
    {
        max1 = fabs(a21);
    }
    if( (max1 < fabs(a31)) )
    {
        max1 = fabs(a31);
    }
    double max2 = fabs(a12);
    if( (max2 < fabs(a13)) )
    {
        max2 = fabs(a13);
    }
    if( (max2 < fabs(a22)) )
    {
        max2 = fabs(a22);
    }
    if( (max2 < fabs(a23)) )
    {
        max2 = fabs(a23);
    }
    double max3 = fabs(a22);
    if( (max3 < fabs(a23)) )
    {
        max3 = fabs(a23);
    }
    if( (max3 < fabs(a32)) )
    {
        max3 = fabs(a32);
    }
    if( (max3 < fabs(a33)) )
    {
        max3 = fabs(a33);
    }
    double lower_bound_1;
    double upper_bound_1;
    int Delta4_sign;
    int int_tmp_result;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 1.63288018496748314939e-98) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 7.23700557733225980357e+75) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (5.11071278299732992696e-15 * ((max2 * max3) * max1));
        if( (Delta4 > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta4 < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    Delta4_sign = int_tmp_result;
    int int_tmp_result_FFWKCAA;
    double max4 = max1;
    if( (max4 < fabs(a41)) )
    {
        max4 = fabs(a41);
    }
    double max5 = max2;
    if( (max5 < max3) )
    {
        max5 = max3;
    }
    double max6 = fabs(a14);
    if( (max6 < fabs(a24)) )
    {
        max6 = fabs(a24);
    }
    if( (max6 < fabs(a34)) )
    {
        max6 = fabs(a34);
    }
    if( (max6 < fabs(a44)) )
    {
        max6 = fabs(a44);
    }
    double max7 = max3;
    if( (max7 < fabs(a42)) )
    {
        max7 = fabs(a42);
    }
    if( (max7 < fabs(a43)) )
    {
        max7 = fabs(a43);
    }
    lower_bound_1 = max4;
    upper_bound_1 = max4;
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max6 < lower_bound_1) )
    {
        lower_bound_1 = max6;
    }
    else
    {
        if( (max6 > upper_bound_1) )
        {
            upper_bound_1 = max6;
        }
    }
    if( (max7 < lower_bound_1) )
    {
        lower_bound_1 = max7;
    }
    else
    {
        if( (max7 > upper_bound_1) )
        {
            upper_bound_1 = max7;
        }
    }
    if( (lower_bound_1 < 2.89273249588395194294e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 7.23700557733225980357e+75) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (3.17768858673611390687e-14 * (((max5 * max7) * max4) * max6));
        if( (r > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (r < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return (Delta4_sign * int_tmp_result_FFWKCAA);
}


/******* extracted from predicates/orient2d.h *******/

inline int orient_2d_filter( const double* p0, const double* p1, const double* p2) {
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double Delta;
    Delta = ((a11 * a22) - (a12 * a21));
    int int_tmp_result;
    double eps;
    double max1 = fabs(a11);
    if( (max1 < fabs(a12)) )
    {
        max1 = fabs(a12);
    }
    double max2 = fabs(a21);
    if( (max2 < fabs(a22)) )
    {
        max2 = fabs(a22);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 5.00368081960964635413e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.67597599124282407923e+153) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.88720573725927976811e-16 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


/******* extracted from predicates/orient3d.h *******/

inline int orient_3d_filter(const double* p0, const double* p1, const double* p2, const double* p3) {
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double a31;
    a31 = (p3[0] - p0[0]);
    double a32;
    a32 = (p3[1] - p0[1]);
    double a33;
    a33 = (p3[2] - p0[2]);
    double Delta;
    Delta = (((a11 * ((a22 * a33) - (a23 * a32))) - (a21 * ((a12 * a33) - (a13 * a32)))) + (a31 * ((a12 * a23) - (a13 * a22))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(a11);
    if((max1 < fabs(a21)))
    {
        max1 = fabs(a21);
    }
    if((max1 < fabs(a31)))
    {
        max1 = fabs(a31);
    }
    double max2 = fabs(a12);
    if((max2 < fabs(a13)))
    {
        max2 = fabs(a13);
    }
    if((max2 < fabs(a22)))
    {
        max2 = fabs(a22);
    }
    if((max2 < fabs(a23)))
    {
        max2 = fabs(a23);
    }
    double max3 = fabs(a22);
    if((max3 < fabs(a23)))
    {
        max3 = fabs(a23);
    }
    if((max3 < fabs(a32)))
    {
        max3 = fabs(a32);
    }
    if((max3 < fabs(a33)))
    {
        max3 = fabs(a33);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if((max2 < lower_bound_1))
    {
        lower_bound_1 = max2;
    }
    else
    {
        if((max2 > upper_bound_1))
        {
            upper_bound_1 = max2;
        }
    }
    if((max3 < lower_bound_1))
    {
        lower_bound_1 = max3;
    }
    else
    {
        if((max3 > upper_bound_1))
        {
            upper_bound_1 = max3;
        }
    }
    if((lower_bound_1 < 1.63288018496748314939e-98))
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if((upper_bound_1 > 5.59936185544450928309e+101))
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (5.11071278299732992696e-15 * ((max2 * max3) * max1));
        if((Delta > eps))
        {
            int_tmp_result = 1;
        }
        else
        {
            if((Delta < -eps))
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


/******* extracted from predicates/dot3d.h *******/

inline int dot_3d_filter( const double* p0, const double* p1, const double* p2) {
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double Delta;
    Delta = (((a11 * a21) + (a12 * a22)) + (a13 * a23));
    int int_tmp_result;
    double eps;
    double max1 = fabs(a11);
    if( (max1 < fabs(a12)) )
    {
        max1 = fabs(a12);
    }
    if( (max1 < fabs(a13)) )
    {
        max1 = fabs(a13);
    }
    double max2 = fabs(a21);
    if( (max2 < fabs(a22)) )
    {
        max2 = fabs(a22);
    }
    if( (max2 < fabs(a23)) )
    {
        max2 = fabs(a23);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 3.78232824369468524638e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.67597599124282407923e+153) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.55534235888797977480e-15 * (max1 * max2));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}

/******* extracted from predicates/dot_compare_3d.h *******/

inline int dot_compare_3d_filter( const double* p0, const double* p1, const double* p2) {
    double d1;
    d1 = (((p0[0] * p1[0]) + (p0[1] * p1[1])) + (p0[2] * p1[2]));
    double d2;
    d2 = (((p0[0] * p2[0]) + (p0[1] * p2[1])) + (p0[2] * p2[2]));
    int int_tmp_result;
    double double_tmp_result;
    double eps;
    double_tmp_result = (d1 - d2);
    double max1 = fabs(p0[0]);
    if( (max1 < fabs(p0[1])) )
    {
        max1 = fabs(p0[1]);
    }
    if( (max1 < fabs(p0[2])) )
    {
        max1 = fabs(p0[2]);
    }
    double max2 = fabs(p1[0]);
    if( (max2 < fabs(p1[1])) )
    {
        max2 = fabs(p1[1]);
    }
    if( (max2 < fabs(p1[2])) )
    {
        max2 = fabs(p1[2]);
    }
    if( (max2 < fabs(p2[0])) )
    {
        max2 = fabs(p2[0]);
    }
    if( (max2 < fabs(p2[1])) )
    {
        max2 = fabs(p2[1]);
    }
    if( (max2 < fabs(p2[2])) )
    {
        max2 = fabs(p2[2]);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 3.01698158319050667656e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.67597599124282407923e+153) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.44455106181954323552e-15 * (max1 * max2));
        if( (double_tmp_result > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (double_tmp_result < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


/******* extracted from predicates/det_compare_4d.h *******/

inline int det_compare_4d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4) {
    double a3_0;
    a3_0 = (p4[0] - p3[0]);
    double a3_1;
    a3_1 = (p4[1] - p3[1]);
    double a3_2;
    a3_2 = (p4[2] - p3[2]);
    double a3_3;
    a3_3 = (p4[3] - p3[3]);
    double m12;
    m12 = ((p1[0] * p0[1]) - (p0[0] * p1[1]));
    double m13;
    m13 = ((p2[0] * p0[1]) - (p0[0] * p2[1]));
    double m14;
    m14 = ((a3_0 * p0[1]) - (p0[0] * a3_1));
    double m23;
    m23 = ((p2[0] * p1[1]) - (p1[0] * p2[1]));
    double m24;
    m24 = ((a3_0 * p1[1]) - (p1[0] * a3_1));
    double m34;
    m34 = ((a3_0 * p2[1]) - (p2[0] * a3_1));
    double m123;
    m123 = (((m23 * p0[2]) - (m13 * p1[2])) + (m12 * p2[2]));
    double m124;
    m124 = (((m24 * p0[2]) - (m14 * p1[2])) + (m12 * a3_2));
    double m134;
    m134 = (((m34 * p0[2]) - (m14 * p2[2])) + (m13 * a3_2));
    double m234;
    m234 = (((m34 * p1[2]) - (m24 * p2[2])) + (m23 * a3_2));
    double Delta;
    Delta = ((((m234 * p0[3]) - (m134 * p1[3])) + (m124 * p2[3])) - (m123 * a3_3));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0[0]);
    if( (max1 < fabs(p1[0])) )
    {
        max1 = fabs(p1[0]);
    }
    if( (max1 < fabs(p2[0])) )
    {
        max1 = fabs(p2[0]);
    }
    if( (max1 < fabs(a3_0)) )
    {
        max1 = fabs(a3_0);
    }
    double max2 = fabs(p0[1]);
    if( (max2 < fabs(p1[1])) )
    {
        max2 = fabs(p1[1]);
    }
    if( (max2 < fabs(p2[1])) )
    {
        max2 = fabs(p2[1]);
    }
    if( (max2 < fabs(a3_1)) )
    {
        max2 = fabs(a3_1);
    }
    double max3 = fabs(p0[2]);
    if( (max3 < fabs(p1[2])) )
    {
        max3 = fabs(p1[2]);
    }
    if( (max3 < fabs(p2[2])) )
    {
        max3 = fabs(p2[2]);
    }
    if( (max3 < fabs(a3_2)) )
    {
        max3 = fabs(a3_2);
    }
    double max4 = fabs(p0[3]);
    if( (max4 < fabs(p1[3])) )
    {
        max4 = fabs(p1[3]);
    }
    if( (max4 < fabs(p2[3])) )
    {
        max4 = fabs(p2[3]);
    }
    if( (max4 < fabs(a3_3)) )
    {
        max4 = fabs(a3_3);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (lower_bound_1 < 3.11018333467425326847e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.44740111546645196071e+76) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.37793769622390420735e-14 * (((max1 * max2) * max3) * max4));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


/******* extracted from predicates/det3d.h *******/

inline int det_3d_filter( const double* p0, const double* p1, const double* p2) {
    double Delta;
    Delta = (((p0[0] * ((p1[1] * p2[2]) - (p1[2] * p2[1]))) - (p1[0] * ((p0[1] * p2[2]) - (p0[2] * p2[1])))) + (p2[0] * ((p0[1] * p1[2]) - (p0[2] * p1[1]))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0[0]);
    if( (max1 < fabs(p1[0])) )
    {
        max1 = fabs(p1[0]);
    }
    if( (max1 < fabs(p2[0])) )
    {
        max1 = fabs(p2[0]);
    }
    double max2 = fabs(p0[1]);
    if( (max2 < fabs(p0[2])) )
    {
        max2 = fabs(p0[2]);
    }
    if( (max2 < fabs(p1[1])) )
    {
        max2 = fabs(p1[1]);
    }
    if( (max2 < fabs(p1[2])) )
    {
        max2 = fabs(p1[2]);
    }
    double max3 = fabs(p1[1]);
    if( (max3 < fabs(p1[2])) )
    {
        max3 = fabs(p1[2]);
    }
    if( (max3 < fabs(p2[1])) )
    {
        max3 = fabs(p2[1]);
    }
    if( (max3 < fabs(p2[2])) )
    {
        max3 = fabs(p2[2]);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 1.92663387981871579179e-98) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.11987237108890185662e+102) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (3.11133555671680765034e-15 * ((max2 * max3) * max1));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}

/******* extracted from predicates/det4d.h *******/

inline int det_4d_filter( const double* p0, const double* p1, const double* p2, const double* p3) {
    double m12;
    m12 = ((p1[0] * p0[1]) - (p0[0] * p1[1]));
    double m13;
    m13 = ((p2[0] * p0[1]) - (p0[0] * p2[1]));
    double m14;
    m14 = ((p3[0] * p0[1]) - (p0[0] * p3[1]));
    double m23;
    m23 = ((p2[0] * p1[1]) - (p1[0] * p2[1]));
    double m24;
    m24 = ((p3[0] * p1[1]) - (p1[0] * p3[1]));
    double m34;
    m34 = ((p3[0] * p2[1]) - (p2[0] * p3[1]));
    double m123;
    m123 = (((m23 * p0[2]) - (m13 * p1[2])) + (m12 * p2[2]));
    double m124;
    m124 = (((m24 * p0[2]) - (m14 * p1[2])) + (m12 * p3[2]));
    double m134;
    m134 = (((m34 * p0[2]) - (m14 * p2[2])) + (m13 * p3[2]));
    double m234;
    m234 = (((m34 * p1[2]) - (m24 * p2[2])) + (m23 * p3[2]));
    double Delta;
    Delta = ((((m234 * p0[3]) - (m134 * p1[3])) + (m124 * p2[3])) - (m123 * p3[3]));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0[0]);
    if( (max1 < fabs(p1[0])) )
    {
        max1 = fabs(p1[0]);
    }
    if( (max1 < fabs(p2[0])) )
    {
        max1 = fabs(p2[0]);
    }
    if( (max1 < fabs(p3[0])) )
    {
        max1 = fabs(p3[0]);
    }
    double max2 = fabs(p0[1]);
    if( (max2 < fabs(p1[1])) )
    {
        max2 = fabs(p1[1]);
    }
    if( (max2 < fabs(p2[1])) )
    {
        max2 = fabs(p2[1]);
    }
    if( (max2 < fabs(p3[1])) )
    {
        max2 = fabs(p3[1]);
    }
    double max3 = fabs(p0[2]);
    if( (max3 < fabs(p1[2])) )
    {
        max3 = fabs(p1[2]);
    }
    if( (max3 < fabs(p2[2])) )
    {
        max3 = fabs(p2[2]);
    }
    if( (max3 < fabs(p3[2])) )
    {
        max3 = fabs(p3[2]);
    }
    double max4 = fabs(p0[3]);
    if( (max4 < fabs(p1[3])) )
    {
        max4 = fabs(p1[3]);
    }
    if( (max4 < fabs(p2[3])) )
    {
        max4 = fabs(p2[3]);
    }
    if( (max4 < fabs(p3[3])) )
    {
        max4 = fabs(p3[3]);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (lower_bound_1 < 3.20402459074399025456e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.44740111546645196071e+76) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.11135406605316806158e-14 * (((max1 * max2) * max3) * max4));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}


/******* extracted from predicates/aligned3d.h *******/

inline int aligned_3d_filter( const double* p0, const double* p1, const double* p2) {
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double delta1;
    delta1 = ((a12 * a23) - (a22 * a13));
    double delta2;
    delta2 = ((a13 * a21) - (a23 * a11));
    double delta3;
    delta3 = ((a11 * a22) - (a21 * a12));
    int int_tmp_result;
    double eps;
    int int_tmp_result_FFWKCAA;
    int int_tmp_result_k60Ocge;
    double max1 = fabs(a12);
    if( (max1 < fabs(a22)) )
    {
        max1 = fabs(a22);
    }
    double max2 = fabs(a13);
    if( (max2 < fabs(a23)) )
    {
        max2 = fabs(a23);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 5.00368081960964635413e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.67597599124282407923e+153) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.88720573725927976811e-16 * (max1 * max2));
        if( (delta1 > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (delta1 < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    double max3 = fabs(a11);
    if( (max3 < fabs(a21)) )
    {
        max3 = fabs(a21);
    }
    lower_bound_1 = max3;
    upper_bound_1 = max3;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (lower_bound_1 < 5.00368081960964635413e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.67597599124282407923e+153) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.88720573725927976811e-16 * (max2 * max3));
        if( (delta2 > eps) )
        {
            int_tmp_result_FFWKCAA = 1;
        }
        else
        {
            if( (delta2 < -eps) )
            {
                int_tmp_result_FFWKCAA = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (lower_bound_1 < 5.00368081960964635413e-147) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.67597599124282407923e+153) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.88720573725927976811e-16 * (max3 * max1));
        if( (delta3 > eps) )
        {
            int_tmp_result_k60Ocge = 1;
        }
        else
        {
            if( (delta3 < -eps) )
            {
                int_tmp_result_k60Ocge = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return ((((int_tmp_result == 0) && (int_tmp_result_FFWKCAA == 0)) && (int_tmp_result_k60Ocge == 0)) ? 0 : 1);
}

/******* extracted from predicates.cpp *******/


// This makes sure the compiler will not optimize y = a*x+b
// with fused multiply-add, this would break the exact
// predicates.
#ifdef GEO_COMPILER_MSVC
#pragma fp_contract(off)
#endif

#include <algorithm>

#define FPG_UNCERTAIN_VALUE 0


#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace {

    using namespace GEO;

    GEO::PCK::SOSMode SOS_mode_ = GEO::PCK::SOS_ADDRESS;

    class LexicoCompare {
    public:

	LexicoCompare(index_t dim) : dim_(dim) {
	}

	bool operator()(const double* x, const double* y) const {
	    for(index_t i=0; i<dim_-1; ++i) {
		if(x[i] < y[i]) {
		    return true;
		}
		if(x[i] > y[i]) {
		    return false;
		}
	    }
	    return (x[dim_-1] < y[dim_-1]);
	}
    private:
	index_t dim_;
    };

    bool lexico_compare_3d(const double* x, const double* y) {
	if(x[0] < y[0]) {
	    return true;
	}
	if(x[0] > y[0]) {
	    return false;
	}
	if(x[1] < y[1]) {
	    return true;
	}
	if(x[1] > y[1]) {
	    return false;
	}
	return x[2] < y[2];
    }

    void SOS_sort(const double** begin, const double** end, index_t dim) {
	if(SOS_mode_ == PCK::SOS_ADDRESS) {
	    std::sort(begin, end);
	} else {
	    if(dim == 3) {
		std::sort(begin, end, lexico_compare_3d);
	    } else {
		std::sort(begin, end, LexicoCompare(dim));
	    }
	}
    }

    inline double max4(double x1, double x2, double x3, double x4) {
#ifdef __SSE2__
	double result;
	__m128d X1 =_mm_load_sd(&x1);
	__m128d X2 =_mm_load_sd(&x2);
	__m128d X3 =_mm_load_sd(&x3);
	__m128d X4 =_mm_load_sd(&x4);
	X1 = _mm_max_sd(X1,X2);
	X3 = _mm_max_sd(X3,X4);
	X1 = _mm_max_sd(X1,X3);
	_mm_store_sd(&result, X1);
	return result;
#else
	return std::max(std::max(x1,x2),std::max(x3,x4));
#endif
    }


    inline void get_minmax3(
	double& m, double& M, double x1, double x2, double x3
    ) {
#ifdef __SSE2__
	__m128d X1 =_mm_load_sd(&x1);
	__m128d X2 =_mm_load_sd(&x2);
	__m128d X3 =_mm_load_sd(&x3);
	__m128d MIN12 = _mm_min_sd(X1,X2);
	__m128d MAX12 = _mm_max_sd(X1,X2);
	X1 = _mm_min_sd(MIN12, X3);
	X3 = _mm_max_sd(MAX12, X3);
	_mm_store_sd(&m, X1);
	_mm_store_sd(&M, X3);
#else
	m = std::min(std::min(x1,x2), x3);
	M = std::max(std::max(x1,x2), x3);
#endif
    }

#ifdef __AVX2__

    inline __m256d avx2_vecdet(__m256d A, __m256d B, __m256d C, __m256d D) {
	__m256d AD = _mm256_mul_pd(A,D);
	__m256d BC = _mm256_mul_pd(B,C);
	return _mm256_sub_pd(AD,BC);
    }

    inline double avx2_det4x4(
	__m256d C11,
	__m256d C12,
	__m256d C13,
	__m256d C14
    ) {
	// We develop w.r.t. the first column and
	// compute the 4 minors simultaneously.

	__m256d C41 = _mm256_permute4x64_pd(C11, _MM_SHUFFLE(2,1,0,3));

	__m256d C22 = _mm256_permute4x64_pd(C12, _MM_SHUFFLE(0,3,2,1));
	__m256d C32 = _mm256_permute4x64_pd(C12, _MM_SHUFFLE(1,0,3,2));

	__m256d C23 = _mm256_permute4x64_pd(C13, _MM_SHUFFLE(0,3,2,1));
	__m256d C33 = _mm256_permute4x64_pd(C13, _MM_SHUFFLE(1,0,3,2));

        __m256d C24 = _mm256_permute4x64_pd(C14, _MM_SHUFFLE(0,3,2,1));
        __m256d C34 = _mm256_permute4x64_pd(C14, _MM_SHUFFLE(1,0,3,2));

	__m256d M1 = _mm256_mul_pd(C12,avx2_vecdet(C23,C24,C33,C34));
	__m256d M2 = _mm256_mul_pd(C22,avx2_vecdet(C13,C14,C33,C34));
	__m256d M3 = _mm256_mul_pd(C32,avx2_vecdet(C13,C14,C23,C24));

	// compute the 4 3x3 minors simulateously by
	// assembling the 3*4 2x2 minors
	__m256d M = _mm256_add_pd(_mm256_sub_pd(M1,M2),M3);

	// multiply the 4 3x3 minors by the 4 coefficients of the
	// first column (permutted so that a41 comes first).
	M = _mm256_mul_pd(M, C41);

	// Compute -m0 +m1 -m2 +m3
	M = _mm256_permute4x64_pd(M, _MM_SHUFFLE(2,0,3,1));
	__m128d M_a = _mm256_extractf128_pd(M, 0);
	__m128d M_b = _mm256_extractf128_pd(M, 1);
	__m128d Mab = _mm_sub_pd(M_a,M_b);
	double m[2];
	_mm_store_pd(m, Mab);
	return m[0]+m[1];
    }

    inline int in_sphere_3d_filter_avx2(
        const double* p, const double* q,
        const double* r, const double* s, const double* t
    ) {

	// Mask to load just three doubles from the points.
	__m256i XYZonly = _mm256_set_epi64x(
	    0,~__int64_t(0),~__int64_t(0),~__int64_t(0)
	);
	__m256d P = _mm256_maskload_pd(p, XYZonly);
	__m256d Q = _mm256_maskload_pd(q, XYZonly);
	__m256d R = _mm256_maskload_pd(r, XYZonly);
	__m256d S = _mm256_maskload_pd(s, XYZonly);
	__m256d T = _mm256_maskload_pd(t, XYZonly);

	__m256d PT = _mm256_sub_pd(P,T);
	__m256d QT = _mm256_sub_pd(Q,T);
	__m256d RT = _mm256_sub_pd(R,T);
	__m256d ST = _mm256_sub_pd(S,T);

	// Absolute values by masking sign bit.
	__m256d sign_mask = _mm256_set1_pd(-0.);
	__m256d absPT     = _mm256_andnot_pd(PT, sign_mask);
	__m256d absQT     = _mm256_andnot_pd(QT, sign_mask);
	__m256d absRT     = _mm256_andnot_pd(RT, sign_mask);
	__m256d absST     = _mm256_andnot_pd(ST, sign_mask);
	__m256d maxXYZ    = _mm256_max_pd(
	    _mm256_max_pd(absPT, absQT), _mm256_max_pd(absRT, absST)
	);

	// Separating maxX, maxY, maxZ in three different registers
	__m128d maxX = _mm256_extractf128_pd(maxXYZ,0);
	__m128d maxZ = _mm256_extractf128_pd(maxXYZ,1);
	__m128d maxY = _mm_shuffle_pd(maxX, maxX, 1);
	__m128d max_max = _mm_max_pd(maxX, _mm_max_pd(maxY, maxZ));

	// Computing dynamic filter
	__m128d eps     = _mm_set1_pd(1.2466136531027298e-13);
	        eps     = _mm_mul_pd(eps, _mm_mul_pd(maxX, _mm_mul_pd(maxY, maxZ)));
		eps     = _mm_mul_pd(eps, _mm_mul_pd(max_max, max_max));

	// Transpose [PT, QT, RT, ST] -> X,Y,Z (last column is 0, ignored)
	__m256d tmp0 = _mm256_shuffle_pd(PT, QT, 0x0);
	__m256d tmp2 = _mm256_shuffle_pd(PT, QT, 0xF);
	__m256d tmp1 = _mm256_shuffle_pd(RT, ST, 0x0);
	__m256d tmp3 = _mm256_shuffle_pd(RT, ST, 0xF);
	__m256d X  = _mm256_permute2f128_pd(tmp0, tmp1, 0x20);
	__m256d Y  = _mm256_permute2f128_pd(tmp2, tmp3, 0x20);
	__m256d Z  = _mm256_permute2f128_pd(tmp0, tmp1, 0x31);

	// Compute first column with squared lengths of vectors
	__m256d X2 = _mm256_mul_pd(X,X);
	__m256d Y2 = _mm256_mul_pd(Y,Y);
	__m256d Z2 = _mm256_mul_pd(Z,Z);
	__m256d L2 = _mm256_add_pd(_mm256_add_pd(X2,Y2),Z2);

	double det = avx2_det4x4(L2,X,Y,Z);

	double epsval;
	_mm_store_pd1(&epsval, eps);

	// Note: inverted as compared to CGAL
	//   CGAL: in_sphere_3d (called side_of_oriented_sphere())
	//      positive side is outside the sphere.
	//   PCK: in_sphere_3d : positive side is inside the sphere

	return (det > epsval) * -1 + (det < -epsval);
    }

#endif

    inline int in_sphere_3d_filter_optim(
        const double* p, const double* q,
        const double* r, const double* s, const double* t
    ) {
        double ptx = p[0] - t[0];
        double pty = p[1] - t[1];
        double ptz = p[2] - t[2];
        double pt2 = geo_sqr(ptx) + geo_sqr(pty) + geo_sqr(ptz);

        double qtx = q[0] - t[0];
        double qty = q[1] - t[1];
        double qtz = q[2] - t[2];
        double qt2 = geo_sqr(qtx) + geo_sqr(qty) + geo_sqr(qtz);

        double rtx = r[0] - t[0];
        double rty = r[1] - t[1];
        double rtz = r[2] - t[2];
        double rt2 = geo_sqr(rtx) + geo_sqr(rty) + geo_sqr(rtz);

        double stx = s[0] - t[0];
        double sty = s[1] - t[1];
        double stz = s[2] - t[2];
        double st2 = geo_sqr(stx) + geo_sqr(sty) + geo_sqr(stz);

        // Compute the semi-static bound.
        double maxx = ::fabs(ptx);
        double maxy = ::fabs(pty);
        double maxz = ::fabs(ptz);

        double aqtx = ::fabs(qtx);
        double artx = ::fabs(rtx);
        double astx = ::fabs(stx);

        double aqty = ::fabs(qty);
        double arty = ::fabs(rty);
        double asty = ::fabs(sty);

        double aqtz = ::fabs(qtz);
        double artz = ::fabs(rtz);
        double astz = ::fabs(stz);

	maxx = max4(maxx, aqtx, artx, astx);
	maxy = max4(maxy, aqty, arty, asty);
	maxz = max4(maxz, aqtz, artz, astz);

        double eps = 1.2466136531027298e-13 * maxx * maxy * maxz;

	double min_max;
	double max_max;
	get_minmax3(min_max, max_max, maxx, maxy, maxz);

        double det = det4x4(
                        ptx,pty,ptz,pt2,
                        rtx,rty,rtz,rt2,
                        qtx,qty,qtz,qt2,
                        stx,sty,stz,st2
                     );

        if (min_max < 1e-58)  { /* sqrt^5(min_double/eps) */
            // Protect against underflow in the computation of eps.
            return FPG_UNCERTAIN_VALUE;
        } else if (max_max < 1e61)  { /* sqrt^5(max_double/4 [hadamard]) */
            // Protect against overflow in the computation of det.
            eps *= (max_max * max_max);
            // Note: inverted as compared to CGAL
            //   CGAL: in_sphere_3d (called side_of_oriented_sphere())
            //      positive side is outside the sphere.
            //   PCK: in_sphere_3d : positive side is inside the sphere
            if (det > eps)  return -1;
            if (det < -eps) return  1;
        }

        return FPG_UNCERTAIN_VALUE;
    }


    using namespace GEO;

    index_t cnt_side1_total = 0;
    index_t cnt_side1_exact = 0;
    index_t cnt_side1_SOS = 0;
    index_t len_side1 = 0;

    index_t cnt_side2_total = 0;
    index_t cnt_side2_exact = 0;
    index_t cnt_side2_SOS = 0;
    index_t len_side2_num = 0;
    index_t len_side2_denom = 0;
    index_t len_side2_SOS = 0;

    index_t cnt_side3_total = 0;
    index_t cnt_side3_exact = 0;
    index_t cnt_side3_SOS = 0;
    index_t len_side3_num = 0;
    index_t len_side3_denom = 0;
    index_t len_side3_SOS = 0;

    index_t cnt_side3h_total = 0;
    index_t cnt_side3h_exact = 0;
    index_t cnt_side3h_SOS = 0;
    index_t len_side3h_num = 0;
    index_t len_side3h_denom = 0;
    index_t len_side3h_SOS = 0;

    index_t cnt_side4_total = 0;
    index_t cnt_side4_exact = 0;
    index_t cnt_side4_SOS = 0;
    index_t len_side4_num = 0;
    index_t len_side4_denom = 0;
    index_t len_side4_SOS = 0;

    index_t cnt_orient2d_total = 0;
    index_t cnt_orient2d_exact = 0;
    index_t len_orient2d = 0;

    index_t cnt_orient3d_total = 0;
    index_t cnt_orient3d_exact = 0;
    index_t len_orient3d = 0;

    index_t cnt_orient3dh_total = 0;
    index_t cnt_orient3dh_exact = 0;
    index_t cnt_orient3dh_SOS = 0;
    index_t len_orient3dh_num = 0;
    index_t len_orient3dh_denom = 0;
    index_t len_orient3dh_SOS = 0;

    // ================= side1 =========================================

    Sign side1_exact_SOS(
        const double* p0, const double* p1,
        const double* q0,
        coord_index_t dim
    ) {
        cnt_side1_exact++;
        expansion& l = expansion_sq_dist(p0, p1, dim);
        expansion& a = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
        expansion& r = expansion_diff(l, a);
        Sign r_sign = r.sign();
        // Symbolic perturbation, Simulation of Simplicity
        if(r_sign == ZERO) {
            cnt_side1_SOS++;
            return (p0 < p1) ? POSITIVE : NEGATIVE;
        }
        len_side1 = std::max(len_side1, r.length());
        return r_sign;
    }

    Sign side1_3d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_3d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 3);
        }
        return result;
    }

    Sign side1_4d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_4d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 4);
        }
        return result;
    }

    Sign side1_6d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_6d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 6);
        }
        return result;
    }

    Sign side1_7d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_7d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 7);
        }
        return result;
    }

    Sign side1_8d_SOS(
        const double* p0, const double* p1, const double* q0
    ) {
        Sign result = Sign(side1_8d_filter(p0, p1, q0));
        if(result == ZERO) {
            result = side1_exact_SOS(p0, p1, q0, 8);
        }
        return result;
    }

    // ================= side2 =========================================

    Sign side2_exact_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1,
        coord_index_t dim
    ) {
        cnt_side2_exact++;

        const expansion& l1 = expansion_sq_dist(p1, p0, dim);
        const expansion& l2 = expansion_sq_dist(p2, p0, dim);

        const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
        const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
        const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
        const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);

        const expansion& Delta = expansion_diff(a11, a10);

        Sign Delta_sign = Delta.sign();
        // Should not occur with symbolic
        // perturbation done at previous steps.
        geo_assert(Delta_sign != ZERO);

        //       [ Lambda0 ]   [ -1 ]        [  a11 ]
        // Delta [         ] = [    ] * l1 + [      ]
        //       [ Lambda1 ]   [  1 ]        [ -a10 ]

        const expansion& DeltaLambda0 = expansion_diff(a11, l1);
        const expansion& DeltaLambda1 = expansion_diff(l1, a10);

        // r = Delta*l2 - ( a20*DeltaLambda0 + a21*DeltaLambda1 )

        const expansion& r0 = expansion_product(Delta, l2);
        const expansion& r1 = expansion_product(a20, DeltaLambda0).negate();
        const expansion& r2 = expansion_product(a21, DeltaLambda1).negate();
        const expansion& r = expansion_sum3(r0, r1, r2);

        Sign r_sign = r.sign();

        // Statistics
        len_side2_num = std::max(len_side2_num, r.length());
        len_side2_denom = std::max(len_side2_denom, Delta.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(r_sign == ZERO) {
            cnt_side2_SOS++;
            const double* p_sort[3];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;

            SOS_sort(p_sort, p_sort + 3, dim);

            for(index_t i = 0; i < 3; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1 = expansion_diff(Delta, a21);
                    const expansion& z = expansion_sum(z1, a20);
                    Sign z_sign = z.sign();
                    len_side2_SOS = std::max(len_side2_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                }
                if(p_sort[i] == p1) {
                    const expansion& z = expansion_diff(a21, a20);
                    Sign z_sign = z.sign();
                    len_side2_SOS = std::max(len_side2_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                }
                if(p_sort[i] == p2) {
                    return NEGATIVE;
                }
            }
            geo_assert_not_reached;
        }

        return Sign(Delta_sign * r_sign);
    }

    Sign side2_3d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_3d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 3);
        }
        return result;
    }

    Sign side2_4d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_4d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 4);
        }
        return result;
    }

    Sign side2_6d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_6d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 6);
        }
        return result;
    }

    Sign side2_7d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_7d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 7);
        }
        return result;
    }

    Sign side2_8d_SOS(
        const double* p0, const double* p1, const double* p2,
        const double* q0, const double* q1
    ) {
        Sign result = Sign(side2_8d_filter(p0, p1, p2, q0, q1));
        if(result == ZERO) {
            result = side2_exact_SOS(p0, p1, p2, q0, q1, 8);
        }
        return result;
    }

    // ================= side3 =========================================

    Sign side3_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2,
        coord_index_t dim
    ) {
        cnt_side3_exact++;

        const expansion& l1 = expansion_sq_dist(p1, p0, dim);
        const expansion& l2 = expansion_sq_dist(p2, p0, dim);
        const expansion& l3 = expansion_sq_dist(p3, p0, dim);

        const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
        const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
        const expansion& a12 = expansion_dot_at(p1, q2, p0, dim).scale_fast(2.0);
        const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
        const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);
        const expansion& a22 = expansion_dot_at(p2, q2, p0, dim).scale_fast(2.0);

        const expansion& a30 = expansion_dot_at(p3, q0, p0, dim).scale_fast(2.0);
        const expansion& a31 = expansion_dot_at(p3, q1, p0, dim).scale_fast(2.0);
        const expansion& a32 = expansion_dot_at(p3, q2, p0, dim).scale_fast(2.0);

        // [ b00 b01 b02 ]           [  1   1   1  ]-1
        // [ b10 b11 b12 ] = Delta * [ a10 a11 a12 ]
        // [ b20 b21 b22 ]           [ a20 a21 a22 ]

        const expansion& b00 = expansion_det2x2(a11, a12, a21, a22);
        const expansion& b01 = expansion_diff(a21, a22);
        const expansion& b02 = expansion_diff(a12, a11);
        const expansion& b10 = expansion_det2x2(a12, a10, a22, a20);
        const expansion& b11 = expansion_diff(a22, a20);
        const expansion& b12 = expansion_diff(a10, a12);
        const expansion& b20 = expansion_det2x2(a10, a11, a20, a21);
        const expansion& b21 = expansion_diff(a20, a21);
        const expansion& b22 = expansion_diff(a11, a10);

        const expansion& Delta = expansion_sum3(b00, b10, b20);
        Sign Delta_sign = Delta.sign();
        // Should not occur with symbolic
        // perturbation done at previous steps.
        geo_assert(Delta_sign != ZERO);

        //       [ Lambda0 ]   [ b01 b02 ]   [ l1 ]   [ b00 ]
        // Delta [ Lambda1 ] = [ b11 b12 ] * [    ] + [ b10 ]
        //       [ Lambda2 ]   [ b21 b22 ]   [ l2 ]   [ b20 ]

        const expansion& b01_l1 = expansion_product(b01, l1);
        const expansion& b02_l2 = expansion_product(b02, l2);
        const expansion& DeltaLambda0 = expansion_sum3(b01_l1, b02_l2, b00);

        const expansion& b11_l1 = expansion_product(b11, l1);
        const expansion& b12_l2 = expansion_product(b12, l2);
        const expansion& DeltaLambda1 = expansion_sum3(b11_l1, b12_l2, b10);

        const expansion& b21_l1 = expansion_product(b21, l1);
        const expansion& b22_l2 = expansion_product(b22, l2);
        const expansion& DeltaLambda2 = expansion_sum3(b21_l1, b22_l2, b20);

        // r = Delta*l3-(a30*DeltaLambda0+a31*DeltaLambda1+a32*DeltaLambda2)

        const expansion& r0 = expansion_product(Delta, l3);
        const expansion& r1 = expansion_product(a30, DeltaLambda0).negate();
        const expansion& r2 = expansion_product(a31, DeltaLambda1).negate();
        const expansion& r3 = expansion_product(a32, DeltaLambda2).negate();
        const expansion& r = expansion_sum4(r0, r1, r2, r3);
        Sign r_sign = r.sign();

        // Statistics
        len_side3_num = std::max(len_side3_num, r.length());
        len_side3_denom = std::max(len_side3_denom, Delta.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(r_sign == ZERO) {
            cnt_side3_SOS++;
            const double* p_sort[4];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            SOS_sort(p_sort, p_sort + 4, dim);
            for(index_t i = 0; i < 4; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1_0 = expansion_sum(b01, b02);
                    const expansion& z1 = expansion_product(a30, z1_0).negate();
                    const expansion& z2_0 = expansion_sum(b11, b12);
                    const expansion& z2 = expansion_product(a31, z2_0).negate();
                    const expansion& z3_0 = expansion_sum(b21, b22);
                    const expansion& z3 = expansion_product(a32, z3_0).negate();
                    const expansion& z = expansion_sum4(Delta, z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3_SOS = std::max(len_side3_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    const expansion& z1 = expansion_product(a30, b01);
                    const expansion& z2 = expansion_product(a31, b11);
                    const expansion& z3 = expansion_product(a32, b21);
                    const expansion& z = expansion_sum3(z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3_SOS = std::max(len_side3_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p2) {
                    const expansion& z1 = expansion_product(a30, b02);
                    const expansion& z2 = expansion_product(a31, b12);
                    const expansion& z3 = expansion_product(a32, b22);
                    const expansion& z = expansion_sum3(z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3_SOS = std::max(len_side3_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p3) {
                    return NEGATIVE;
                }
            }
            geo_assert_not_reached;
        }
        return Sign(Delta_sign * r_sign);
    }


    Sign side3h_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        double h0, double h1, double h2, double h3,
        const double* q0, const double* q1, const double* q2
    ) {
        cnt_side3h_exact++;

        const expansion& l1 = expansion_diff(h1,h0);
        const expansion& l2 = expansion_diff(h2,h0);
        const expansion& l3 = expansion_diff(h3,h0);

        const expansion& a10 = expansion_dot_at(p1, q0, p0, 3).scale_fast(2.0);
        const expansion& a11 = expansion_dot_at(p1, q1, p0, 3).scale_fast(2.0);
        const expansion& a12 = expansion_dot_at(p1, q2, p0, 3).scale_fast(2.0);
        const expansion& a20 = expansion_dot_at(p2, q0, p0, 3).scale_fast(2.0);
        const expansion& a21 = expansion_dot_at(p2, q1, p0, 3).scale_fast(2.0);
        const expansion& a22 = expansion_dot_at(p2, q2, p0, 3).scale_fast(2.0);

        const expansion& a30 = expansion_dot_at(p3, q0, p0, 3).scale_fast(2.0);
        const expansion& a31 = expansion_dot_at(p3, q1, p0, 3).scale_fast(2.0);
        const expansion& a32 = expansion_dot_at(p3, q2, p0, 3).scale_fast(2.0);

        // [ b00 b01 b02 ]           [  1   1   1  ]-1
        // [ b10 b11 b12 ] = Delta * [ a10 a11 a12 ]
        // [ b20 b21 b22 ]           [ a20 a21 a22 ]

        const expansion& b00 = expansion_det2x2(a11, a12, a21, a22);
        const expansion& b01 = expansion_diff(a21, a22);
        const expansion& b02 = expansion_diff(a12, a11);
        const expansion& b10 = expansion_det2x2(a12, a10, a22, a20);
        const expansion& b11 = expansion_diff(a22, a20);
        const expansion& b12 = expansion_diff(a10, a12);
        const expansion& b20 = expansion_det2x2(a10, a11, a20, a21);
        const expansion& b21 = expansion_diff(a20, a21);
        const expansion& b22 = expansion_diff(a11, a10);

        const expansion& Delta = expansion_sum3(b00, b10, b20);
        Sign Delta_sign = Delta.sign();
        // Should not occur with symbolic
        // perturbation done at previous steps.
        geo_assert(Delta_sign != ZERO);

        //       [ Lambda0 ]   [ b01 b02 ]   [ l1 ]   [ b00 ]
        // Delta [ Lambda1 ] = [ b11 b12 ] * [    ] + [ b10 ]
        //       [ Lambda2 ]   [ b21 b22 ]   [ l2 ]   [ b20 ]

        const expansion& b01_l1 = expansion_product(b01, l1);
        const expansion& b02_l2 = expansion_product(b02, l2);
        const expansion& DeltaLambda0 = expansion_sum3(b01_l1, b02_l2, b00);

        const expansion& b11_l1 = expansion_product(b11, l1);
        const expansion& b12_l2 = expansion_product(b12, l2);
        const expansion& DeltaLambda1 = expansion_sum3(b11_l1, b12_l2, b10);

        const expansion& b21_l1 = expansion_product(b21, l1);
        const expansion& b22_l2 = expansion_product(b22, l2);
        const expansion& DeltaLambda2 = expansion_sum3(b21_l1, b22_l2, b20);

        // r = Delta*l3-(a30*DeltaLambda0+a31*DeltaLambda1+a32*DeltaLambda2)

        const expansion& r0 = expansion_product(Delta, l3);
        const expansion& r1 = expansion_product(a30, DeltaLambda0).negate();
        const expansion& r2 = expansion_product(a31, DeltaLambda1).negate();
        const expansion& r3 = expansion_product(a32, DeltaLambda2).negate();
        const expansion& r = expansion_sum4(r0, r1, r2, r3);
        Sign r_sign = r.sign();

        // Statistics
        len_side3h_num = std::max(len_side3h_num, r.length());
        len_side3h_denom = std::max(len_side3h_denom, Delta.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(r_sign == ZERO) {
            cnt_side3h_SOS++;
            const double* p_sort[4];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;

	    SOS_sort(p_sort, p_sort + 4, 3);
            for(index_t i = 0; i < 4; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1_0 = expansion_sum(b01, b02);
                    const expansion& z1 = expansion_product(a30, z1_0).negate();
                    const expansion& z2_0 = expansion_sum(b11, b12);
                    const expansion& z2 = expansion_product(a31, z2_0).negate();
                    const expansion& z3_0 = expansion_sum(b21, b22);
                    const expansion& z3 = expansion_product(a32, z3_0).negate();
                    const expansion& z = expansion_sum4(Delta, z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3h_SOS = std::max(len_side3h_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    const expansion& z1 = expansion_product(a30, b01);
                    const expansion& z2 = expansion_product(a31, b11);
                    const expansion& z3 = expansion_product(a32, b21);
                    const expansion& z = expansion_sum3(z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3h_SOS = std::max(len_side3h_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p2) {
                    const expansion& z1 = expansion_product(a30, b02);
                    const expansion& z2 = expansion_product(a31, b12);
                    const expansion& z3 = expansion_product(a32, b22);
                    const expansion& z = expansion_sum3(z1, z2, z3);
                    Sign z_sign = z.sign();
                    len_side3h_SOS = std::max(len_side3h_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p3) {
                    return NEGATIVE;
                }
            }
            geo_assert_not_reached;
        }
        return Sign(Delta_sign * r_sign);
    }


    Sign side3_3d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_3d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 3);
        }
        return result;
    }


    Sign side3_4d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_4d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 4);
        }
        return result;
    }

    Sign side3_6d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_6d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 6);
        }
        return result;
    }

    Sign side3_7d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_7d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 7);
        }
        return result;
    }

    Sign side3_8d_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* q0, const double* q1, const double* q2
    ) {
        Sign result = Sign(side3_8d_filter(p0, p1, p2, p3, q0, q1, q2));
        if(result == ZERO) {
            result = side3_exact_SOS(p0, p1, p2, p3, q0, q1, q2, 8);
        }
        return result;
    }

    // ================= side4 =========================================

    Sign side4_3d_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* p4, bool sos = true
    ) {
        cnt_side4_exact++;

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);
        const expansion& a14 = expansion_sq_dist(p1, p0, 3).negate();

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);
        const expansion& a24 = expansion_sq_dist(p2, p0, 3).negate();

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);
        const expansion& a34 = expansion_sq_dist(p3, p0, 3).negate();

        const expansion& a41 = expansion_diff(p4[0], p0[0]);
        const expansion& a42 = expansion_diff(p4[1], p0[1]);
        const expansion& a43 = expansion_diff(p4[2], p0[2]);
        const expansion& a44 = expansion_sq_dist(p4, p0, 3).negate();

        // This commented-out version does not reuse
        // the 2x2 minors.
/*
        const expansion& Delta1 = expansion_det3x3(
            a21, a22, a23,
            a31, a32, a33,
            a41, a42, a43
        );
        const expansion& Delta2 = expansion_det3x3(
            a11, a12, a13,
            a31, a32, a33,
            a41, a42, a43
        );
        const expansion& Delta3 = expansion_det3x3(
            a11, a12, a13,
            a21, a22, a23,
            a41, a42, a43
        );
        const expansion& Delta4 = expansion_det3x3(
            a11, a12, a13,
            a21, a22, a23,
            a31, a32, a33
        );
*/

        // Optimized version that reuses the 2x2 minors

        const expansion& m12 = expansion_det2x2(a12,a13,a22,a23);
        const expansion& m13 = expansion_det2x2(a12,a13,a32,a33);
        const expansion& m14 = expansion_det2x2(a12,a13,a42,a43);
        const expansion& m23 = expansion_det2x2(a22,a23,a32,a33);
        const expansion& m24 = expansion_det2x2(a22,a23,a42,a43);
        const expansion& m34 = expansion_det2x2(a32,a33,a42,a43);


        const expansion& z11 = expansion_product(a21,m34);
        const expansion& z12 = expansion_product(a31,m24).negate();
        const expansion& z13 = expansion_product(a41,m23);
        const expansion& Delta1 = expansion_sum3(z11,z12,z13);

        const expansion& z21 = expansion_product(a11,m34);
        const expansion& z22 = expansion_product(a31,m14).negate();
        const expansion& z23 = expansion_product(a41,m13);
        const expansion& Delta2 = expansion_sum3(z21,z22,z23);

        const expansion& z31 = expansion_product(a11,m24);
        const expansion& z32 = expansion_product(a21,m14).negate();
        const expansion& z33 = expansion_product(a41,m12);
        const expansion& Delta3 = expansion_sum3(z31,z32,z33);

        const expansion& z41 = expansion_product(a11,m23);
        const expansion& z42 = expansion_product(a21,m13).negate();
        const expansion& z43 = expansion_product(a31,m12);
        const expansion& Delta4 = expansion_sum3(z41,z42,z43);


        Sign Delta4_sign = Delta4.sign();
        geo_assert(Delta4_sign != ZERO);

        const expansion& r_1 = expansion_product(Delta1, a14);
        const expansion& r_2 = expansion_product(Delta2, a24).negate();
        const expansion& r_3 = expansion_product(Delta3, a34);
        const expansion& r_4 = expansion_product(Delta4, a44).negate();
        const expansion& r = expansion_sum4(r_1, r_2, r_3, r_4);
        Sign r_sign = r.sign();

        // Statistics
        len_side4_num = std::max(len_side4_num, r.length());
        len_side4_denom = std::max(len_side4_denom, Delta1.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(sos && r_sign == ZERO) {
            cnt_side4_SOS++;
            const double* p_sort[5];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            p_sort[4] = p4;
            SOS_sort(p_sort, p_sort + 5, 3);
            for(index_t i = 0; i < 5; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1 = expansion_diff(Delta2, Delta1);
                    const expansion& z2 = expansion_diff(Delta4, Delta3);
                    const expansion& z = expansion_sum(z1, z2);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta4_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    Sign Delta1_sign = Delta1.sign();
                    if(Delta1_sign != ZERO) {
                        len_side4_SOS = std::max(len_side4_SOS, Delta1.length());
                        return Sign(Delta4_sign * Delta1_sign);
                    }
                } else if(p_sort[i] == p2) {
                    Sign Delta2_sign = Delta2.sign();
                    if(Delta2_sign != ZERO) {
                        len_side4_SOS = std::max(len_side4_SOS, Delta2.length());
                        return Sign(-Delta4_sign * Delta2_sign);
                    }
                } else if(p_sort[i] == p3) {
                    Sign Delta3_sign = Delta3.sign();
                    if(Delta3_sign != ZERO) {
                        len_side4_SOS = std::max(len_side4_SOS, Delta3.length());
                        return Sign(Delta4_sign * Delta3_sign);
                    }
                } else if(p_sort[i] == p4) {
                    return NEGATIVE;
                }
            }
        }
        return Sign(Delta4_sign * r_sign);
    }

    Sign side4_exact_SOS(
        const double* p0, const double* p1, const double* p2, const double* p3,
        const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3,
        coord_index_t dim
    ) {
        cnt_side4_exact++;

        const expansion& l1 = expansion_sq_dist(p1, p0, dim);
        const expansion& l2 = expansion_sq_dist(p2, p0, dim);
        const expansion& l3 = expansion_sq_dist(p3, p0, dim);
        const expansion& l4 = expansion_sq_dist(p4, p0, dim);

        const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
        const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
        const expansion& a12 = expansion_dot_at(p1, q2, p0, dim).scale_fast(2.0);
        const expansion& a13 = expansion_dot_at(p1, q3, p0, dim).scale_fast(2.0);

        const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
        const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);
        const expansion& a22 = expansion_dot_at(p2, q2, p0, dim).scale_fast(2.0);
        const expansion& a23 = expansion_dot_at(p2, q3, p0, dim).scale_fast(2.0);

        const expansion& a30 = expansion_dot_at(p3, q0, p0, dim).scale_fast(2.0);
        const expansion& a31 = expansion_dot_at(p3, q1, p0, dim).scale_fast(2.0);
        const expansion& a32 = expansion_dot_at(p3, q2, p0, dim).scale_fast(2.0);
        const expansion& a33 = expansion_dot_at(p3, q3, p0, dim).scale_fast(2.0);

        const expansion& a40 = expansion_dot_at(p4, q0, p0, dim).scale_fast(2.0);
        const expansion& a41 = expansion_dot_at(p4, q1, p0, dim).scale_fast(2.0);
        const expansion& a42 = expansion_dot_at(p4, q2, p0, dim).scale_fast(2.0);
        const expansion& a43 = expansion_dot_at(p4, q3, p0, dim).scale_fast(2.0);

        // [ b00 b01 b02 b03 ]           [  1   1   1   1  ]-1
        // [ b10 b11 b12 b13 ]           [ a10 a11 a12 a13 ]
        // [ b20 b21 b22 b23 ] = Delta * [ a20 a21 a22 a23 ]
        // [ b30 b31 b32 b33 ]           [ a30 a31 a32 a33 ]

        // Note: we could probably reuse some of the co-factors
        // (but for now I'd rather keep this form that is easier to
        //  read ... and to debug if need be !)

        const expansion& b00 = expansion_det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33);
        const expansion& b01 = expansion_det_111_2x3(a21, a22, a23, a31, a32, a33).negate();
        const expansion& b02 = expansion_det_111_2x3(a11, a12, a13, a31, a32, a33);
        const expansion& b03 = expansion_det_111_2x3(a11, a12, a13, a21, a22, a23).negate();

        const expansion& b10 = expansion_det3x3(a10, a12, a13, a20, a22, a23, a30, a32, a33).negate();
        const expansion& b11 = expansion_det_111_2x3(a20, a22, a23, a30, a32, a33);
        const expansion& b12 = expansion_det_111_2x3(a10, a12, a13, a30, a32, a33).negate();
        const expansion& b13 = expansion_det_111_2x3(a10, a12, a13, a20, a22, a23);

        const expansion& b20 = expansion_det3x3(a10, a11, a13, a20, a21, a23, a30, a31, a33);
        const expansion& b21 = expansion_det_111_2x3(a20, a21, a23, a30, a31, a33).negate();
        const expansion& b22 = expansion_det_111_2x3(a10, a11, a13, a30, a31, a33);
        const expansion& b23 = expansion_det_111_2x3(a10, a11, a13, a20, a21, a23).negate();

        const expansion& b30 = expansion_det3x3(a10, a11, a12, a20, a21, a22, a30, a31, a32).negate();
        const expansion& b31 = expansion_det_111_2x3(a20, a21, a22, a30, a31, a32);
        const expansion& b32 = expansion_det_111_2x3(a10, a11, a12, a30, a31, a32).negate();
        const expansion& b33 = expansion_det_111_2x3(a10, a11, a12, a20, a21, a22);

        const expansion& Delta = expansion_sum4(b00, b10, b20, b30);
        Sign Delta_sign = Delta.sign();
        geo_assert(Delta_sign != ZERO);

        //       [ Lambda0 ]   [ b01 b02 b03 ]   [ l1 ]   [ b00 ]
        //       [ Lambda1 ]   [ b11 b12 b13 ]   [ l2 ]   [ b10 ]
        // Delta [ Lambda2 ] = [ b21 b22 b23 ] * [ l3 ] + [ b20 ]
        //       [ Lambda3 ]   [ b31 b32 b33 ]   [ l4 ]   [ b30 ]

        const expansion& b01_l1 = expansion_product(b01, l1);
        const expansion& b02_l2 = expansion_product(b02, l2);
        const expansion& b03_l3 = expansion_product(b03, l3);
        const expansion& DeltaLambda0 = expansion_sum4(b01_l1, b02_l2, b03_l3, b00);

        const expansion& b11_l1 = expansion_product(b11, l1);
        const expansion& b12_l2 = expansion_product(b12, l2);
        const expansion& b13_l3 = expansion_product(b13, l3);
        const expansion& DeltaLambda1 = expansion_sum4(b11_l1, b12_l2, b13_l3, b10);

        const expansion& b21_l1 = expansion_product(b21, l1);
        const expansion& b22_l2 = expansion_product(b22, l2);
        const expansion& b23_l3 = expansion_product(b23, l3);
        const expansion& DeltaLambda2 = expansion_sum4(b21_l1, b22_l2, b23_l3, b20);

        const expansion& b31_l1 = expansion_product(b31, l1);
        const expansion& b32_l2 = expansion_product(b32, l2);
        const expansion& b33_l3 = expansion_product(b33, l3);
        const expansion& DeltaLambda3 = expansion_sum4(b31_l1, b32_l2, b33_l3, b30);

        // r = Delta*l4 - (
        //    a40*DeltaLambda0+
        //    a41*DeltaLambda1+
        //    a42*DeltaLambda2+
        //    a43*DeltaLambda3
        // )

        const expansion& r0 = expansion_product(Delta, l4);
        const expansion& r1 = expansion_product(a40, DeltaLambda0);
        const expansion& r2 = expansion_product(a41, DeltaLambda1);
        const expansion& r3 = expansion_product(a42, DeltaLambda2);
        const expansion& r4 = expansion_product(a43, DeltaLambda3);
        const expansion& r1234 = expansion_sum4(r1, r2, r3, r4);
        const expansion& r = expansion_diff(r0, r1234);
        Sign r_sign = r.sign();

        // Simulation of Simplicity (symbolic perturbation)
        if(r_sign == ZERO) {
            cnt_side4_SOS++;
            const double* p_sort[5];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            p_sort[4] = p4;
            SOS_sort(p_sort, p_sort + 5, dim);
            for(index_t i = 0; i < 5; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1_0 = expansion_sum3(b01, b02, b03);
                    const expansion& z1 = expansion_product(a30, z1_0);
                    const expansion& z2_0 = expansion_sum3(b11, b12, b13);
                    const expansion& z2 = expansion_product(a31, z2_0);
                    const expansion& z3_0 = expansion_sum3(b21, b22, b23);
                    const expansion& z3 = expansion_product(a32, z3_0);
                    const expansion& z4_0 = expansion_sum3(b31, b32, b33);
                    const expansion& z4 = expansion_product(a33, z4_0);
                    const expansion& z1234 = expansion_sum4(z1, z2, z3, z4);
                    const expansion& z = expansion_diff(Delta, z1234);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    const expansion& z1 = expansion_product(a30, b01);
                    const expansion& z2 = expansion_product(a31, b11);
                    const expansion& z3 = expansion_product(a32, b21);
                    const expansion& z4 = expansion_product(a33, b31);
                    const expansion& z = expansion_sum4(z1, z2, z3, z4);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p2) {
                    const expansion& z1 = expansion_product(a30, b02);
                    const expansion& z2 = expansion_product(a31, b12);
                    const expansion& z3 = expansion_product(a32, b22);
                    const expansion& z4 = expansion_product(a33, b32);
                    const expansion& z = expansion_sum4(z1, z2, z3, z4);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p3) {
                    const expansion& z1 = expansion_product(a30, b03);
                    const expansion& z2 = expansion_product(a31, b13);
                    const expansion& z3 = expansion_product(a32, b23);
                    const expansion& z4 = expansion_product(a33, b33);
                    const expansion& z = expansion_sum4(z1, z2, z3, z4);
                    Sign z_sign = z.sign();
                    len_side4_SOS = std::max(len_side4_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta_sign * z_sign);
                    }
                } else if(p_sort[i] == p4) {
                    return NEGATIVE;
                }
            }
            geo_assert_not_reached;
        }
        return Sign(r_sign * Delta_sign);
    }

    Sign side4_4d_SOS(
        const double* p0,
        const double* p1, const double* p2, const double* p3, const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3
    ) {
        Sign result = Sign(side4_4d_filter(p0, p1, p2, p3, p4, q0, q1, q2, q3));
        if(result == ZERO) {
            result = side4_exact_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3, 4);
        }
        return result;
    }

    Sign side4_6d_SOS(
        const double* p0,
        const double* p1, const double* p2, const double* p3, const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3
    ) {
        Sign result = Sign(side4_6d_filter(p0, p1, p2, p3, p4, q0, q1, q2, q3));
        if(result == ZERO) {
            result = side4_exact_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3, 6);
        }
        return result;
    }

    Sign side4_7d_SOS(
        const double* p0,
        const double* p1, const double* p2, const double* p3, const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3
    ) {
        Sign result = Sign(side4_7d_filter(p0, p1, p2, p3, p4, q0, q1, q2, q3));
        if(result == ZERO) {
            result = side4_exact_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3, 7);
        }
        return result;
    }

    Sign side4_8d_SOS(
        const double* p0,
        const double* p1, const double* p2, const double* p3, const double* p4,
        const double* q0, const double* q1, const double* q2, const double* q3
    ) {
        Sign result = Sign(side4_8d_filter(p0, p1, p2, p3, p4, q0, q1, q2, q3));
        if(result == ZERO) {
            result = side4_exact_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3, 8);
        }
        return result;
    }

    // ============ orient2d ==============================================

    Sign orient_2d_exact(
        const double* p0, const double* p1, const double* p2
    ) {
        cnt_orient2d_exact++;

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);

        const expansion& Delta = expansion_det2x2(
            a11, a12, a21, a22
        );

        len_orient2d = std::max(len_orient2d, Delta.length());

        return Delta.sign();
    }


    // ============ orient3d ==============================================

    Sign orient_3d_exact(
        const double* p0, const double* p1,
        const double* p2, const double* p3
    ) {
        cnt_orient3d_exact++;

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);

        const expansion& Delta = expansion_det3x3(
            a11, a12, a13, a21, a22, a23, a31, a32, a33
        );

        len_orient3d = std::max(len_orient3d, Delta.length());

        return Delta.sign();
    }

    Sign side4h_3d_exact_SOS(
        const double* p0, const double* p1,
        const double* p2, const double* p3, const double* p4,
        double h0, double h1, double h2, double h3, double h4,
        bool sos = true
    ) {
        cnt_orient3dh_exact++;

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(p1[2], p0[2]);
        const expansion& a14 = expansion_diff(h0,h1);

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(p2[2], p0[2]);
        const expansion& a24 = expansion_diff(h0,h2);

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(p3[2], p0[2]);
        const expansion& a34 = expansion_diff(h0,h3);

        const expansion& a41 = expansion_diff(p4[0], p0[0]);
        const expansion& a42 = expansion_diff(p4[1], p0[1]);
        const expansion& a43 = expansion_diff(p4[2], p0[2]);
        const expansion& a44 = expansion_diff(h0,h4);

        // Note: we could probably reuse some of the 2x2 co-factors
        // (but for now I'd rather keep this form that is easier to
        //  read ... and to debug if need be !)
        const expansion& Delta1 = expansion_det3x3(
            a21, a22, a23,
            a31, a32, a33,
            a41, a42, a43
        );
        const expansion& Delta2 = expansion_det3x3(
            a11, a12, a13,
            a31, a32, a33,
            a41, a42, a43
        );
        const expansion& Delta3 = expansion_det3x3(
            a11, a12, a13,
            a21, a22, a23,
            a41, a42, a43
        );
        const expansion& Delta4 = expansion_det3x3(
            a11, a12, a13,
            a21, a22, a23,
            a31, a32, a33
        );

        Sign Delta4_sign = Delta4.sign();
        geo_assert(Delta4_sign != ZERO);

        const expansion& r_1 = expansion_product(Delta1, a14);
        const expansion& r_2 = expansion_product(Delta2, a24).negate();
        const expansion& r_3 = expansion_product(Delta3, a34);
        const expansion& r_4 = expansion_product(Delta4, a44).negate();
        const expansion& r = expansion_sum4(r_1, r_2, r_3, r_4);

        Sign r_sign = r.sign();

        // Statistics
        len_orient3dh_num = std::max(len_orient3dh_num, r.length());
        len_orient3dh_denom = std::max(len_orient3dh_denom, Delta1.length());

        // Simulation of Simplicity (symbolic perturbation)
        if(sos && r_sign == ZERO) {
            cnt_orient3dh_SOS++;
            const double* p_sort[5];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            p_sort[4] = p4;

	    SOS_sort(p_sort, p_sort + 5, 3);
            for(index_t i = 0; i < 5; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1 = expansion_diff(Delta2, Delta1);
                    const expansion& z2 = expansion_diff(Delta4, Delta3);
                    const expansion& z = expansion_sum(z1, z2);
                    Sign z_sign = z.sign();
                    len_orient3dh_SOS = std::max(len_orient3dh_SOS, z.length());
                    if(z_sign != ZERO) {
                        return Sign(Delta4_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    Sign Delta1_sign = Delta1.sign();
                    if(Delta1_sign != ZERO) {
                        len_orient3dh_SOS = std::max(len_orient3dh_SOS, Delta1.length());
                        return Sign(Delta4_sign * Delta1_sign);
                    }
                } else if(p_sort[i] == p2) {
                    Sign Delta2_sign = Delta2.sign();
                    if(Delta2_sign != ZERO) {
                        len_orient3dh_SOS = std::max(len_orient3dh_SOS, Delta2.length());
                        return Sign(-Delta4_sign * Delta2_sign);
                    }
                } else if(p_sort[i] == p3) {
                    Sign Delta3_sign = Delta3.sign();
                    if(Delta3_sign != ZERO) {
                        len_orient3dh_SOS = std::max(len_orient3dh_SOS, Delta3.length());
                        return Sign(Delta4_sign * Delta3_sign);
                    }
                } else if(p_sort[i] == p4) {
                    return NEGATIVE;
                }
            }
        }
        return Sign(Delta4_sign * r_sign);
    }


    Sign side3h_2d_exact_SOS(
        const double* p0, const double* p1,
	const double* p2, const double* p3,
        double h0, double h1, double h2, double h3,
        bool sos = true
    ) {

        const expansion& a11 = expansion_diff(p1[0], p0[0]);
        const expansion& a12 = expansion_diff(p1[1], p0[1]);
        const expansion& a13 = expansion_diff(h0,h1);

        const expansion& a21 = expansion_diff(p2[0], p0[0]);
        const expansion& a22 = expansion_diff(p2[1], p0[1]);
        const expansion& a23 = expansion_diff(h0,h2);

        const expansion& a31 = expansion_diff(p3[0], p0[0]);
        const expansion& a32 = expansion_diff(p3[1], p0[1]);
        const expansion& a33 = expansion_diff(h0,h3);

        const expansion& Delta1 = expansion_det2x2(
            a21, a22,
            a31, a32
        );
        const expansion& Delta2 = expansion_det2x2(
            a11, a12,
            a31, a32
        );
        const expansion& Delta3 = expansion_det2x2(
            a11, a12,
            a21, a22
        );

        Sign Delta3_sign = Delta3.sign();
        geo_assert(Delta3_sign != ZERO);

        const expansion& r_1 = expansion_product(Delta1, a13);
        const expansion& r_2 = expansion_product(Delta2, a23).negate();
        const expansion& r_3 = expansion_product(Delta3, a33);
        const expansion& r = expansion_sum3(r_1, r_2, r_3);

        Sign r_sign = r.sign();

        // Simulation of Simplicity (symbolic perturbation)
        if(sos && r_sign == ZERO) {
            const double* p_sort[4];
            p_sort[0] = p0;
            p_sort[1] = p1;
            p_sort[2] = p2;
            p_sort[3] = p3;
            SOS_sort(p_sort, p_sort + 4, 3);
            for(index_t i = 0; i < 4; ++i) {
                if(p_sort[i] == p0) {
                    const expansion& z1 = expansion_diff(Delta2, Delta1);
                    const expansion& z = expansion_sum(z1, Delta3);
                    Sign z_sign = z.sign();
                    if(z_sign != ZERO) {
                        return Sign(Delta3_sign * z_sign);
                    }
                } else if(p_sort[i] == p1) {
                    Sign Delta1_sign = Delta1.sign();
                    if(Delta1_sign != ZERO) {
                        return Sign(Delta3_sign * Delta1_sign);
                    }
                } else if(p_sort[i] == p2) {
                    Sign Delta2_sign = Delta2.sign();
                    if(Delta2_sign != ZERO) {
                        return Sign(-Delta3_sign * Delta2_sign);
                    }
                } else if(p_sort[i] == p3) {
		    return NEGATIVE;
                }
            }
        }
        return Sign(Delta3_sign * r_sign);
    }


    // ================================ det and dot =======================

    Sign det_3d_exact(
	const double* p0, const double* p1, const double* p2
    ) {
	const expansion& p0_0 = expansion_create(p0[0]);
	const expansion& p0_1 = expansion_create(p0[1]);
	const expansion& p0_2 = expansion_create(p0[2]);

	const expansion& p1_0 = expansion_create(p1[0]);
	const expansion& p1_1 = expansion_create(p1[1]);
	const expansion& p1_2 = expansion_create(p1[2]);

	const expansion& p2_0 = expansion_create(p2[0]);
	const expansion& p2_1 = expansion_create(p2[1]);
	const expansion& p2_2 = expansion_create(p2[2]);

	const expansion& Delta = expansion_det3x3(
	    p0_0, p0_1, p0_2,
	    p1_0, p1_1, p1_2,
	    p2_0, p2_1, p2_2
	);
	return Delta.sign();
    }


    bool aligned_3d_exact(
	const double* p0, const double* p1, const double* p2
    ) {
	const expansion& U_0 = expansion_diff(p1[0],p0[0]);
	const expansion& U_1 = expansion_diff(p1[1],p0[1]);
	const expansion& U_2 = expansion_diff(p1[2],p0[2]);

	const expansion& V_0 = expansion_diff(p2[0],p0[0]);
	const expansion& V_1 = expansion_diff(p2[1],p0[1]);
	const expansion& V_2 = expansion_diff(p2[2],p0[2]);

	const expansion& N_0 = expansion_det2x2(U_1, V_1, U_2, V_2);
	const expansion& N_1 = expansion_det2x2(U_2, V_2, U_0, V_0);
	const expansion& N_2 = expansion_det2x2(U_0, V_0, U_1, V_1);

	return(
	    N_0.sign() == 0 &&
	    N_1.sign() == 0 &&
	    N_2.sign() == 0
	);
    }

    Sign dot_3d_exact(
	const double* p0, const double* p1, const double* p2
    ) {
	const expansion& U_0 = expansion_diff(p1[0],p0[0]);
	const expansion& U_1 = expansion_diff(p1[1],p0[1]);
	const expansion& U_2 = expansion_diff(p1[2],p0[2]);

	const expansion& V_0 = expansion_diff(p2[0],p0[0]);
	const expansion& V_1 = expansion_diff(p2[1],p0[1]);
	const expansion& V_2 = expansion_diff(p2[2],p0[2]);

	const expansion& UV_0 = expansion_product(U_0, V_0);
	const expansion& UV_1 = expansion_product(U_1, V_1);
	const expansion& UV_2 = expansion_product(U_2, V_2);

	const expansion& Delta = expansion_sum3(UV_0, UV_1, UV_2);

	return Delta.sign();
    }

    Sign dot_compare_3d_exact(
	const double* v0, const double* v1, const double* v2
    ) {
	const expansion& d01_0 = expansion_product(v0[0], v1[0]);
	const expansion& d01_1 = expansion_product(v0[1], v1[1]);
	const expansion& d01_2 = expansion_product(v0[2], v1[2]);
	const expansion& d01_12 = expansion_sum(d01_1, d01_2);
	const expansion& d01 = expansion_sum(d01_0, d01_12);

	const expansion& d02_0 = expansion_product(v0[0], v2[0]);
	const expansion& d02_1 = expansion_product(v0[1], v2[1]);
	const expansion& d02_2 = expansion_product(v0[2], v2[2]);
	const expansion& d02_12 = expansion_sum(d02_1, d02_2);
	const expansion& d02 = expansion_sum(d02_0, d02_12);

	const expansion& result = expansion_diff(d01, d02);

	return result.sign();
    }

    // ================================ statistics ========================

    inline double percent(index_t a, index_t b) {
        if(a == 0 && b == 0) {
            return 0;
        }
        return double(a) * 100.0 / double(b);
    }

    void show_stats_plain(
        const std::string& name, index_t cnt1, index_t cnt2
    ) {
        Logger::out(name)
            << "Tot:" << cnt1
            << " Exact:" << cnt2
            << std::endl;
        Logger::out(name)
            << " Exact: " << percent(cnt2, cnt1) << "% "
            << std::endl;
    }

    void show_stats_sos(
        const std::string& name, index_t cnt1, index_t cnt2, index_t cnt3
    ) {
        Logger::out(name)
            << "Tot:" << cnt1
            << " Exact:" << cnt2
            << " SOS:" << cnt3 << std::endl;
        Logger::out(name)
            << " Exact: " << percent(cnt2, cnt1) << "% "
            << " SOS: " << percent(cnt3, cnt1) << "% "
            << std::endl;
    }

    void show_stats_sos(
        const std::string& name, index_t cnt1, index_t cnt2, index_t cnt3,
        index_t len
    ) {
        show_stats_sos(name, cnt1, cnt2, cnt3);
        Logger::out(name) << " Len: " << len << std::endl;
    }

    void show_stats_plain(
        const std::string& name, index_t cnt1, index_t cnt2,
        index_t len
    ) {
        show_stats_plain(name, cnt1, cnt2);
        Logger::out(name) << " Len: " << len << std::endl;
    }

    void show_stats_sos(
        const std::string& name, index_t cnt1, index_t cnt2, index_t cnt3,
        index_t num_len, index_t denom_len, index_t SOS_len
    ) {
        show_stats_sos(name, cnt1, cnt2, cnt3);
        Logger::out(name)
            << " Num len: " << num_len
            << " Denom len: " << denom_len
            << " SOS len: " << SOS_len
            << std::endl;
    }
}



namespace GEO {

    namespace PCK {

	void set_SOS_mode(SOSMode m) {
	    SOS_mode_ = m;
	}

	SOSMode get_SOS_mode() {
	    return SOS_mode_;
	}


        Sign side1_SOS(
            const double* p0, const double* p1,
            const double* q0,
            coord_index_t DIM
        ) {
            cnt_side1_total++;
            switch(DIM) {
            case 3:
                return side1_3d_SOS(p0, p1, q0);
            case 4:
                return side1_4d_SOS(p0, p1, q0);
            case 6:
                return side1_6d_SOS(p0, p1, q0);
            case 7:
                return side1_7d_SOS(p0, p1, q0);
            case 8:
                return side1_8d_SOS(p0, p1, q0);
            }
            geo_assert_not_reached;
        }

        Sign side2_SOS(
            const double* p0, const double* p1, const double* p2,
            const double* q0, const double* q1,
            coord_index_t DIM
        ) {
            cnt_side2_total++;
            switch(DIM) {
            case 3:
                return side2_3d_SOS(p0, p1, p2, q0, q1);
            case 4:
                return side2_4d_SOS(p0, p1, p2, q0, q1);
            case 6:
                return side2_6d_SOS(p0, p1, p2, q0, q1);
            case 7:
                return side2_7d_SOS(p0, p1, p2, q0, q1);
            case 8:
                return side2_8d_SOS(p0, p1, p2, q0, q1);
            }
            geo_assert_not_reached;
        }

        Sign side3_SOS(
            const double* p0, const double* p1, const double* p2, const double* p3,
            const double* q0, const double* q1, const double* q2,
            coord_index_t DIM
        ) {
            cnt_side3_total++;
            switch(DIM) {
            case 3:
                return side3_3d_SOS(p0, p1, p2, p3, q0, q1, q2);
            case 4:
                return side3_4d_SOS(p0, p1, p2, p3, q0, q1, q2);
            case 6:
                return side3_6d_SOS(p0, p1, p2, p3, q0, q1, q2);
            case 7:
                return side3_7d_SOS(p0, p1, p2, p3, q0, q1, q2);
            case 8:
                return side3_8d_SOS(p0, p1, p2, p3, q0, q1, q2);
            }
            geo_assert_not_reached;
        }


        Sign side3_3dlifted_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3,
            double h0, double h1, double h2, double h3,
            const double* q0, const double* q1, const double* q2,
	    bool SOS
        ) {
            Sign result = Sign(side3h_3d_filter(p0, p1, p2, p3, h0, h1, h2, h3, q0, q1, q2));
            if(SOS && result == ZERO) {
                result = side3h_exact_SOS(p0, p1, p2, p3, h0, h1, h2, h3, q0, q1, q2);
            }
            return result;
        }

        Sign side4_SOS(
            const double* p0,
            const double* p1, const double* p2, const double* p3, const double* p4,
            const double* q0, const double* q1, const double* q2, const double* q3,
            coord_index_t DIM
        ) {
            switch(DIM) {
            case 3:
                // 3d is a special case for side4()
                //   (intrinsic dim == ambient dim)
                // therefore embedding tet q0,q1,q2,q3 is not needed.
                // WARNING: cnt_side4_total is not incremented here,
                // because it is
                // incremented in side4_3d_SOS().
                return side4_3d_SOS(p0, p1, p2, p3, p4);
            case 4:
                cnt_side4_total++;
                return side4_4d_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3);
            case 6:
                cnt_side4_total++;
                return side4_6d_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3);
            case 7:
                cnt_side4_total++;
                return side4_7d_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3);
            case 8:
                cnt_side4_total++;
                return side4_8d_SOS(p0, p1, p2, p3, p4, q0, q1, q2, q3);
            }
            geo_assert_not_reached;
        }


        Sign side4_3d(
            const double* p0, const double* p1, const double* p2, const double* p3,
            const double* p4
        ) {
            cnt_side4_total++;
            Sign result = Sign(side4_3d_filter(p0, p1, p2, p3, p4));
            if(result == 0) {
                // last argument is false: do not apply symbolic perturbation
                result = side4_3d_exact_SOS(p0, p1, p2, p3, p4, false);
            }
            return result;
        }

        Sign side4_3d_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3,
            const double* p4
        ) {
            cnt_side4_total++;
            Sign result = Sign(side4_3d_filter(p0, p1, p2, p3, p4));
            if(result == 0) {
                result = side4_3d_exact_SOS(p0, p1, p2, p3, p4);
            }
            return result;
        }


        Sign in_sphere_3d_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3,
            const double* p4
        ) {
            // in_sphere_3d is simply implemented using side4_3d.
            // Both predicates are equivalent through duality as can
            // be easily seen:
            // side4_3d(p0,p1,p2,p3,p4) returns POSITIVE if
            //    d(q,p0) < d(q,p4)
            //    where q denotes the circumcenter of (p0,p1,p2,p3)
            // Note that d(q,p0) = R  (radius of circumscribed sphere)
            // In other words, side4_3d(p0,p1,p2,p3,p4) returns POSITIVE if
            //   d(q,p4) > R which means whenever p4 is not in the
            //   circumscribed sphere of (p0,p1,p2,p3).
            // Therefore:
            // in_sphere_3d(p0,p1,p2,p3,p4) = -side4_3d(p0,p1,p2,p3,p4)

            cnt_side4_total++;

            // This specialized filter supposes that orient_3d(p0,p1,p2,p3) > 0

#ifdef __AVX2__
	    Sign result = Sign(in_sphere_3d_filter_avx2(p0, p1, p2, p3, p4));
#else
            Sign result = Sign(in_sphere_3d_filter_optim(p0, p1, p2, p3, p4));
#endif
            if(result == 0) {
                result = side4_3d_exact_SOS(p0, p1, p2, p3, p4);
            }
            return Sign(-result);
        }

        Sign GEOGRAM_API in_circle_2d_SOS(
            const double* p0, const double* p1, const double* p2,
            const double* p3
        ) {
            // in_circle_2d is simply implemented using side3_2d.
            // Both predicates are equivalent through duality as can
            // be easily seen:
            // side3_2d(p0,p1,p2,p3,p0,p1,p2) returns POSITIVE if
            //    d(q,p0) < d(q,p3)
            //    where q denotes the circumcenter of (p0,p1,p2)
            // Note that d(q,p0) = R  (radius of circumscribed circle)
            // In other words, side3_2d(p0,p1,p2,p3,p4) returns POSITIVE if
            //   d(q,p3) > R which means whenever p3 is not in the
            //   circumscribed circle of (p0,p1,p2).
            // Therefore:
            // in_circle_2d(p0,p1,p2,p3) = -side3_2d(p0,p1,p2,p3)

	    // TODO: implement specialized filter like the one used
	    // by "in-sphere".
	    Sign s = Sign(-side3_2d_filter(p0, p1, p2, p3, p0, p1, p2));
	    if(s != ZERO) {
		return s;
	    }
	    return Sign(-side3_exact_SOS(p0, p1, p2, p3, p0, p1, p2, 2));
        }

        Sign GEOGRAM_API in_circle_3d_SOS(
            const double* p0, const double* p1, const double* p2,
            const double* p3
        ) {
            // in_circle_3d is simply implemented using side3_3d.
            // Both predicates are equivalent through duality as can
            // be easily seen:
            // side3_3d(p0,p1,p2,p3,p0,p1,p2) returns POSITIVE if
            //    d(q,p0) < d(q,p3)
            //    where q denotes the circumcenter of (p0,p1,p2)
            // Note that d(q,p0) = R  (radius of circumscribed circle)
            // In other words, side3_3d(p0,p1,p2,p3,p4) returns POSITIVE if
            //   d(q,p3) > R which means whenever p3 is not in the
            //   circumscribed circle of (p0,p1,p2).
            // Therefore:
            // in_circle_3d(p0,p1,p2,p3) = -side3_3d(p0,p1,p2,p3)
            return Sign(-side3_3d_SOS(p0,p1,p2,p3,p0,p1,p2));
        }

        Sign GEOGRAM_API in_circle_3dlifted_SOS(
            const double* p0, const double* p1, const double* p2,
            const double* p3,
            double h0, double h1, double h2, double h3,
	    bool SOS
        ) {
            // in_circle_3dlifted is simply implemented using side3_3dlifted.
            // Both predicates are equivalent through duality
            // (see comment in in_circle_3d_SOS(), the same
            //  remark applies).
            return Sign(-side3_3dlifted_SOS(p0,p1,p2,p3,h0,h1,h2,h3,p0,p1,p2,SOS));
        }


        Sign orient_2d(
            const double* p0, const double* p1, const double* p2
        ) {
            cnt_orient2d_total++;
            Sign result = Sign(orient_2d_filter(p0, p1, p2));
            if(result == 0) {
                result = orient_2d_exact(p0, p1, p2);
            }
            return result;
        }


        Sign orient_2dlifted_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3,
            double h0, double h1, double h2, double h3
	) {
            Sign result = Sign(
                side3_2dlifted_2d_filter(
                    p0, p1, p2, p3, h0, h1, h2, h3
                    )
                );
            if(result == 0) {
                result = side3h_2d_exact_SOS(
                    p0, p1, p2, p3, h0, h1, h2, h3
                );
            }
            // orient_3d() is opposite to side3h()
            // (like in_sphere() that is opposite to side3())
	    return result;
	}


        Sign orient_3d(
            const double* p0, const double* p1,
            const double* p2, const double* p3
            ) {
            cnt_orient3d_total++;
            Sign result = Sign(orient_3d_filter(p0, p1, p2, p3));
            if(result == 0) {
                result = orient_3d_exact(p0, p1, p2, p3);
            }
            return result;
        }


        Sign orient_3dlifted(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4,
            double h0, double h1, double h2, double h3, double h4
        ) {
            cnt_orient3dh_total++;
            Sign result = Sign(
                side4h_3d_filter(
                    p0, p1, p2, p3, p4, h0, h1, h2, h3, h4
                    )
                );
            if(result == 0) {
                // last argument is false -> do not perturb.
                result = side4h_3d_exact_SOS(
                    p0, p1, p2, p3, p4, h0, h1, h2, h3, h4, false
                );
            }
            // orient_4d() is opposite to side4h()
            // (like in_sphere() that is opposite to side4())
            return Sign(-result);
        }

        Sign orient_3dlifted_SOS(
            const double* p0, const double* p1,
            const double* p2, const double* p3, const double* p4,
            double h0, double h1, double h2, double h3, double h4
        ) {
            cnt_orient3dh_total++;
            Sign result = Sign(
                side4h_3d_filter(
                    p0, p1, p2, p3, p4, h0, h1, h2, h3, h4
                    )
                );
            if(result == 0) {
                result = side4h_3d_exact_SOS(
                    p0, p1, p2, p3, p4, h0, h1, h2, h3, h4
                );
            }
            // orient_4d() is opposite to side4h()
            // (like in_sphere() that is opposite to side4())
            return Sign(-result);
        }

	Sign det_3d(
	    const double* p0, const double* p1, const double* p2
	) {
	    Sign result = Sign(
		det_3d_filter(p0, p1, p2)
	    );
	    if(result == 0) {
		result = det_3d_exact(p0, p1, p2);
	    }
	    return result;
	}


	Sign det_4d(
	    const double* p0, const double* p1, const double* p2, const double* p3
	) {
	    Sign result = Sign(
		det_4d_filter(p0, p1, p2, p3)
	    );

	    if(result == 0) {
		const expansion& p0_0 = expansion_create(p0[0]);
		const expansion& p0_1 = expansion_create(p0[1]);
		const expansion& p0_2 = expansion_create(p0[2]);
		const expansion& p0_3 = expansion_create(p0[3]);

		const expansion& p1_0 = expansion_create(p1[0]);
		const expansion& p1_1 = expansion_create(p1[1]);
		const expansion& p1_2 = expansion_create(p1[2]);
		const expansion& p1_3 = expansion_create(p1[3]);

		const expansion& p2_0 = expansion_create(p2[0]);
		const expansion& p2_1 = expansion_create(p2[1]);
		const expansion& p2_2 = expansion_create(p2[2]);
		const expansion& p2_3 = expansion_create(p2[3]);

		const expansion& p3_0 = expansion_create(p3[0]);
		const expansion& p3_1 = expansion_create(p3[1]);
		const expansion& p3_2 = expansion_create(p3[2]);
		const expansion& p3_3 = expansion_create(p3[3]);

		result = sign_of_expansion_determinant(
		    p0_0, p0_1, p0_2, p0_3,
		    p1_0, p1_1, p1_2, p1_3,
		    p2_0, p2_1, p2_2, p2_3,
		    p3_0, p3_1, p3_2, p3_3
		);
	    }
	    return result;
	}


	Sign det_compare_4d(
	    const double* p0, const double* p1,
	    const double* p2, const double* p3,
	    const double* p4
	) {
	    Sign result = Sign(
		det_compare_4d_filter(p0, p1, p2, p3, p4)
	    );
	    if(result == 0) {
		const expansion& p0_0 = expansion_create(p0[0]);
		const expansion& p0_1 = expansion_create(p0[1]);
		const expansion& p0_2 = expansion_create(p0[2]);
		const expansion& p0_3 = expansion_create(p0[3]);

		const expansion& p1_0 = expansion_create(p1[0]);
		const expansion& p1_1 = expansion_create(p1[1]);
		const expansion& p1_2 = expansion_create(p1[2]);
		const expansion& p1_3 = expansion_create(p1[3]);

		const expansion& p2_0 = expansion_create(p2[0]);
		const expansion& p2_1 = expansion_create(p2[1]);
		const expansion& p2_2 = expansion_create(p2[2]);
		const expansion& p2_3 = expansion_create(p2[3]);

		const expansion& a3_0 = expansion_diff(p4[0],p3[0]);
		const expansion& a3_1 = expansion_diff(p4[1],p3[1]);
		const expansion& a3_2 = expansion_diff(p4[2],p3[2]);
		const expansion& a3_3 = expansion_diff(p4[3],p3[3]);

		result = sign_of_expansion_determinant(
		    p0_0, p0_1, p0_2, p0_3,
		    p1_0, p1_1, p1_2, p1_3,
		    p2_0, p2_1, p2_2, p2_3,
		    a3_0, a3_1, a3_2, a3_3
		);
	    }
	    return result;
	}


	bool aligned_3d(
	    const double* p0, const double* p1, const double* p2
	) {
	    /*
	    Sign result = Sign(
		aligned_3d_filter(p0,p1,p2)
	    );
	    if(result != 0) {
		return false;
	    }
	    */
	    return aligned_3d_exact(p0, p1, p2);
	}

	Sign dot_3d(
	    const double* p0, const double* p1, const double* p2
	) {
	    Sign result = Sign(det_3d_filter(p0, p1, p2));
	    if(result == 0) {
		result = dot_3d_exact(p0, p1, p2);
	    }
	    return result;
	}

	Sign dot_compare_3d(
	    const double* v0, const double* v1, const double* v2
	) {
	    Sign result = Sign(dot_compare_3d_filter(v0, v1, v2));
	    if(result == 0) {
		result = dot_compare_3d_exact(v0, v1, v2);
	    }
	    return result;
	}


	bool points_are_identical_2d(
	    const double* p1,
	    const double* p2
	) {
	    return
		(p1[0] == p2[0]) &&
		(p1[1] == p2[1])
	    ;
	}

	bool points_are_identical_3d(
	    const double* p1,
	    const double* p2
	) {
	    return
		(p1[0] == p2[0]) &&
		(p1[1] == p2[1]) &&
		(p1[2] == p2[2])
	    ;
	}

	bool points_are_colinear_3d(
	    const double* p1,
	    const double* p2,
	    const double* p3
	) {
	    // Colinearity is tested by using four coplanarity
	    // tests with four points that are not coplanar.
	    // TODO: use PCK::aligned_3d() instead (to be tested)
	    static const double q000[3] = {0.0, 0.0, 0.0};
	    static const double q001[3] = {0.0, 0.0, 1.0};
	    static const double q010[3] = {0.0, 1.0, 0.0};
	    static const double q100[3] = {1.0, 0.0, 0.0};
	    return
		PCK::orient_3d(p1, p2, p3, q000) == ZERO &&
		PCK::orient_3d(p1, p2, p3, q001) == ZERO &&
		PCK::orient_3d(p1, p2, p3, q010) == ZERO &&
		PCK::orient_3d(p1, p2, p3, q100) == ZERO
	    ;
	}

    inline Sign orient_3d_inexact(
	    const double* p0, const double* p1,
	    const double* p2, const double* p3
	) {
	    double a11 = p1[0] - p0[0] ;
	    double a12 = p1[1] - p0[1] ;
	    double a13 = p1[2] - p0[2] ;

	    double a21 = p2[0] - p0[0] ;
	    double a22 = p2[1] - p0[1] ;
	    double a23 = p2[2] - p0[2] ;

	    double a31 = p3[0] - p0[0] ;
	    double a32 = p3[1] - p0[1] ;
	    double a33 = p3[2] - p0[2] ;

	    double Delta = det3x3(
		a11,a12,a13,
		a21,a22,a23,
		a31,a32,a33
	    );

	    return geo_sgn(Delta);
	}

        void initialize() {
            expansion::initialize();
        }

        void terminate() {
            // Nothing to do.
        }

        void show_stats() {
            show_stats_plain(
                "orient2d",
                cnt_orient2d_total, cnt_orient2d_exact,
                len_orient2d
            );
            show_stats_plain(
                "orient3d",
                cnt_orient3d_total, cnt_orient3d_exact,
                len_orient3d
            );
            show_stats_sos(
                "orient3dh",
                cnt_orient3dh_total, cnt_orient3dh_exact, cnt_orient3dh_SOS,
                len_orient3dh_num, len_orient3dh_denom, len_orient3dh_SOS
            );
            show_stats_sos(
                "side1",
                cnt_side1_total, cnt_side1_exact, cnt_side1_SOS,
                len_side1
            );
            show_stats_sos(
                "side2",
                cnt_side2_total, cnt_side2_exact, cnt_side2_SOS,
                len_side2_num, len_side2_denom, len_side2_SOS
            );
            show_stats_sos(
                "side3",
                cnt_side3_total, cnt_side3_exact, cnt_side3_SOS,
                len_side3_num, len_side3_denom, len_side3_SOS
            );
            show_stats_sos(
                "side3h",
                cnt_side3h_total, cnt_side3h_exact, cnt_side3h_SOS,
                len_side3h_num, len_side3h_denom, len_side3h_SOS
            );
            show_stats_sos(
                "side4/insph.",
                cnt_side4_total, cnt_side4_exact, cnt_side4_SOS,
                len_side4_num, len_side4_denom, len_side4_SOS
            );
        }
    }
}

