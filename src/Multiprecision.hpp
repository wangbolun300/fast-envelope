#pragma once

#include <gmp.h>
#include <iostream>

namespace fastEnvelope {

	class Multiprecision
	{
	public:
		static void set_precision(const int digits)
		{
			mpf_set_default_prec(digits);
		}

		static int get_precision()
		{
			return mpf_get_default_prec();
		}

		mpf_t value;

		int get_sign()
		{
			return mpf_sgn(value);
		}
		int get_prec_bits() {
			return mpf_get_prec(value);
		}
		Multiprecision()
		{
			mpf_init(value);
			mpf_set_d(value, 0);
		}

		Multiprecision(double d)
		{
			mpf_init(value);
			mpf_set_d(value, d);
		}

		Multiprecision(double d, int precision)
		{
			mpf_init2(value, precision);
			mpf_set_d(value, d);
		}

		Multiprecision(const mpf_t &v_)
		{
			mpf_init2(value, mpf_get_prec(v_));
			mpf_set(value, v_);
		}

		Multiprecision(const Multiprecision &other)
		{
			mpf_init2(value, mpf_get_prec(other.value));
			mpf_set(value, other.value);
		}

		~Multiprecision()
		{
			mpf_clear(value);
		}

		friend Multiprecision operator+(const Multiprecision &x, const Multiprecision &y)
		{
			static Multiprecision r_out;
			mpf_add(r_out.value, x.value, y.value);
			return r_out;
		}

		friend Multiprecision operator-(const Multiprecision &x, const Multiprecision &y)
		{
			static Multiprecision r_out;
			mpf_sub(r_out.value, x.value, y.value);
			return r_out;
		}

		friend Multiprecision operator*(const Multiprecision &x, const Multiprecision &y)
		{
			static Multiprecision r_out;
			mpf_mul(r_out.value, x.value, y.value);
			return r_out;
		}

		friend Multiprecision operator/(const Multiprecision &x, const Multiprecision &y)
		{
			static Multiprecision r_out;
			mpf_div(r_out.value, x.value, y.value);
			return r_out;
		}

		Multiprecision &operator=(const Multiprecision &x)
		{
			if (this == &x)
				return *this;
			//mpf_init2(value, prec);
			mpf_set(value, x.value);
			return *this;
		}

		Multiprecision &operator=(const double x)
		{
			//mpf_init2(value, prec);
			mpf_set_d(value, x);
			return *this;
		}

		//> < ==
		friend bool operator<(const Multiprecision &r, const Multiprecision &r1)
		{
			return mpf_cmp(r.value, r1.value) < 0;
		}

		friend bool operator>(const Multiprecision &r, const Multiprecision &r1)
		{
			return mpf_cmp(r.value, r1.value) > 0;
		}

		friend bool operator<=(const Multiprecision &r, const Multiprecision &r1)
		{
			return mpf_cmp(r.value, r1.value) <= 0;
		}

		friend bool operator>=(const Multiprecision &r, const Multiprecision &r1)
		{
			return mpf_cmp(r.value, r1.value) >= 0;
		}

		friend bool operator==(const Multiprecision &r, const Multiprecision &r1)
		{
			return mpf_cmp(r.value, r1.value) == 0;
		}

		friend bool operator!=(const Multiprecision &r, const Multiprecision &r1)
		{
			return mpf_cmp(r.value, r1.value) != 0;
		}

		//to double
		double to_double()
		{
			return mpf_get_d(value);
		}

		//<<
		friend std::ostream &operator<<(std::ostream &os, const Multiprecision &r)
		{
			os << mpf_get_d(r.value);
			return os;
		}

		Multiprecision sqrt(const Multiprecision &mp)
		{
			Multiprecision res;
			mpf_sqrt(res.value, mp.value);

			return res;
		}
	};
}
