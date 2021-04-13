#ifndef BOOST_MP_QD_BACKEND_HPP
#define BOOST_MP_QD_BACKEND_HPP
#define _QD_SPLIT 134217729.0                // = 2^27 + 1
#define _QD_SPLIT_MAX 6.69692879491417e+299  // = 2^996

#include <iostream>
#include <limits>
#include <cmath>

namespace boost {
	namespace multiprecision {
		namespace backends {



		} // namespace backends ends
		namespace qd {
			static const double _d_nan = std::numeric_limits<double>::quiet_NaN();
			static const double _d_inf = std::numeric_limits<double>::infinity();

			// Basic Operations

			// func(a + b) and err(a + b)
			// Assumes |a| >= |b|
			inline double quick_two_sum(double a, double b, double &err)
			{
				double s = a + b;
				err = b - (s - a);
				return s;
			}

			// func(a - b) and err(a - b)
			// Assumes |a| >= |b|
			inline double quick_two_diff(double a, double b, double &err)
			{
				double s = a - b;
				err = (a - s) - b;
				return s;
			}

			// func(a + b) and err(a + b)
			inline double two_sum(double a, double b, double &err)
			{
				double s = a + b;
				double bb = s - a;
				err = (a - (s - bb)) + (b - bb);
				return s;
			}

			// func(a - b) and err(a - b)
			inline double two_diff(double a, double b, double &err)
			{
				double s = a - b;
				double bb = s - a;
				err = (a - (s - bb)) - (b + bb);
				return s;
			}


#ifndef QD_FMS

			inline void split(double a, double &high, double &low)
			{
				double temp;
				if (a > _QD_SPLIT_MAX || a < (-1 * _QD_SPLIT_MAX))
				{
					a *= 3.7252902984619140625e-09;  // 2^-28
					temp = _QD_SPLIT * a;
					high = temp - (temp - a);
					low = a - high;
					high *= 268435456.0;            // 2^28
					low *= 268435456.0;             // 2^28 
				}
				else
				{
					temp = _QD_SPLIT * a;
					high = temp - (temp - a);
					low = a - high;
				}
			}

#endif // !QD_FMS


			// func(a * b) and err(a * b)
			inline double two_product(double a, double b, double &err)
			{

#ifdef QD_FMS
				double product = a * b;
				err = split(a, b, p);
				return product;
#else
				double a_high, a_low, b_high, b_low;
				double product = a * b;
				split(a, a_high, a_low);
				split(b, b_high, b_low);
				err = ((a_high * b_high) - product) + (a_high * b_low) + (a_low * b_low);
				return product;
#endif // QD_FMS


				// func(a * a) and err(a * a), Faster but may lose some bits
				QD_API double two_square(double a, double &err)
				{
#ifdef QD_FMS
					double square = a * a;
					err = split(a, a, square);
					return square;
#else
					double high, low;
					double square = a * a;
					split(a, high, low);
					err = ((high * high) - square + (2.0 * high * low) + (low * low));
					return square;
#endif // QD_FMS

				}


				// finds nearest integer for double d
				inline double nearest_int(double d)
				{
					if (d == std::floor(d))
						return d;
					return std::floor(d + 0.5);
				}

				// finds the truncated integer
				inline double trunc_int(double d)
				{
					return (d >= 0.0) ? std::floor(d) : std::ceil(d);
				}

				inline double sqr(double t)
				{
					return t * t;
				}

				inline double to_double(double a)
				{
					return a;
				}

				inline int to_int(double a)
				{
					return static_cast<int>(a);
				}

#ifndef _QD_QD_COLLECTION
#define _QD_QD_COLLECTION

#ifndef QD_API
#undef QD_API
#endif

				//
				// If fused multiply-add is available, define to correct macro for using it.
				// It is invoked as QD_FMA(a, b, c) to compute func(a * b + c).
				// If correctly rounded multiply-add is not available (or if unsure), keep it undefined.
				//
#ifndef QD_FMA
#undef QD_FMA
#endif

//
// If fused multiply-subtract is available, define to correct macro for using it.
// It is invoked as QD_FMS(a, b, c) to compute func(a * b - c).
// If correctly rounded multiply-add is not available (or if unsure), keep it undefined.
//
#ifndef QD_FMS
#undef QD_FMS
#endif

//
// Set the following to 1 to define commonly used function to be inlined. 
// This should be set to 1 unless the compiler does not support the "inline" keyword,
// or if building for debugging purposes.
//
#ifndef QD_INLINE
#undef QD_INLINE
#endif

//
// Set the following to 1 to use ANSI C++ standard header files
// such as cmath, iostream, etc.  If set to zero, it will try to
// include math.h, iostream.h, etc, instead.
//
#ifndef QD_HAVE_STD
#undef QD_HAVE_STD 
#endif

//
// Set the following to 1 to make the addition and subtraction to satisfy the IEEE-style error bound
// func(a + b) = (1 + d) * (a + b)
// where |d| <= eps.  If set to 0, the addition and subtraction
// will satisfy the weaker Cray-style error bound
// func(a + b) = (1 + d1) * a + (1 + d2) * b
// where |d1| <= eps and |d2| eps.
//
#ifndef QD_IEEE_ADD
#undef QD_IEEE_ADD
#endif

// Set the following to 1 to use slightly inaccurate but faster version of multiplication.
#ifndef QD_SLOPPY_MUL
#undef QD_SLOPPY_MUL
#endif

// Set the following to 1 to use slightly inaccurate but faster version of division.
#ifndef QD_SLOPPY_DIV
#undef QD_SLOPPY_DIV
#endif

// Define this macro to be the isfinite(x) function.
#ifndef QD_ISFINITE
#undef QD_ISFINITE
#endif

// Define this macro to be the isinf(x) function.
#ifndef QD_ISINF
#undef QD_ISINF
#endif

// Define this macro to be the isnan(x) function.
#ifndef QD_ISNAN
#undef QD_ISNAN
#endif

#endif // !_QD_QD_COLLECTION



#ifndef _QD_QD_REAL
#define _QD_QD_REAL

				struct QD_API qd_real {
					double x[4];    // The components of quad double

					inline qd_real::qd_real(double x0, double x1, double x2, double x3) {
						x[0] = x0;
						x[1] = x1;
						x[2] = x2;
						x[3] = x3;
					}

					inline qd_real::qd_real(const double *xx) {
						x[0] = xx[0];
						x[1] = xx[1];
						x[2] = xx[2];
						x[3] = xx[3];
					}

					inline qd_real::qd_real(double x0) {
						x[0] = x0;
						x[1] = x[2] = x[3] = 0.0;
					}

					inline qd_real::qd_real() {
						x[0] = 0.0;
						x[1] = 0.0;
						x[2] = 0.0;
						x[3] = 0.0;
					}

					inline qd_real::qd_real(int i) {
						x[0] = static_cast<double>(i);
						x[1] = x[2] = x[3] = 0.0;
					}

					// Eliminates any zeros in the middle component(s).
					void zero_elim();
					void zero_elim(double &e);

					//
					// Renormalisation
					//
					inline void qd_real::renorm() {
						qd::renorm(x[0], x[1], x[2], x[3]);
					}

					inline void qd_real::renorm(double &e) {
						qd::renorm(x[0], x[1], x[2], x[3], e);
					}

					inline void quick_renorm(double &c0, double &c1,
						double &c2, double &c3, double &c4) {
						double t0, t1, t2, t3;
						double s;
						s = qd::quick_two_sum(c3, c4, t3);
						s = qd::quick_two_sum(c2, s, t2);
						s = qd::quick_two_sum(c1, s, t1);
						c0 = qd::quick_two_sum(c0, s, t0);

						s = qd::quick_two_sum(t2, t3, t2);
						s = qd::quick_two_sum(t1, s, t1);
						c1 = qd::quick_two_sum(t0, s, t0);

						s = qd::quick_two_sum(t1, t2, t1);
						c2 = qd::quick_two_sum(t0, s, t0);

						c3 = t0 + t1;
					}

					inline void renorm(double &c0, double &c1,
						double &c2, double &c3) {
						double s0, s1, s2 = 0.0, s3 = 0.0;

						if (QD_ISINF(c0)) return;

						s0 = qd::quick_two_sum(c2, c3, c3);
						s0 = qd::quick_two_sum(c1, s0, c2);
						c0 = qd::quick_two_sum(c0, s0, c1);

						s0 = c0;
						s1 = c1;
						if (s1 != 0.0) {
							s1 = qd::quick_two_sum(s1, c2, s2);
							if (s2 != 0.0)
								s2 = qd::quick_two_sum(s2, c3, s3);
							else
								s1 = qd::quick_two_sum(s1, c3, s2);
						}
						else {
							s0 = qd::quick_two_sum(s0, c2, s1);
							if (s1 != 0.0)
								s1 = qd::quick_two_sum(s1, c3, s2);
							else
								s0 = qd::quick_two_sum(s0, c3, s1);
						}

						c0 = s0;
						c1 = s1;
						c2 = s2;
						c3 = s3;
					}

					inline void renorm(double &c0, double &c1,
						double &c2, double &c3, double &c4) {
						double s0, s1, s2 = 0.0, s3 = 0.0;

						if (QD_ISINF(c0)) return;

						s0 = qd::quick_two_sum(c3, c4, c4);
						s0 = qd::quick_two_sum(c2, s0, c3);
						s0 = qd::quick_two_sum(c1, s0, c2);
						c0 = qd::quick_two_sum(c0, s0, c1);

						s0 = c0;
						s1 = c1;

						s0 = qd::quick_two_sum(c0, c1, s1);
						if (s1 != 0.0) {
							s1 = qd::quick_two_sum(s1, c2, s2);
							if (s2 != 0.0) {
								s2 = qd::quick_two_sum(s2, c3, s3);
								if (s3 != 0.0)
									s3 += c4;
								else
									s2 += c4;
							}
							else {
								s1 = qd::quick_two_sum(s1, c3, s2);
								if (s2 != 0.0)
									s2 = qd::quick_two_sum(s2, c4, s3);
								else
									s1 = qd::quick_two_sum(s1, c4, s2);
							}
						}
						else {
							s0 = qd::quick_two_sum(s0, c2, s1);
							if (s1 != 0.0) {
								s1 = qd::quick_two_sum(s1, c3, s2);
								if (s2 != 0.0)
									s2 = qd::quick_two_sum(s2, c4, s3);
								else
									s1 = qd::quick_two_sum(s1, c4, s2);
							}
							else {
								s0 = qd::quick_two_sum(s0, c3, s1);
								if (s1 != 0.0)
									s1 = qd::quick_two_sum(s1, c4, s2);
								else
									s0 = qd::quick_two_sum(s0, c4, s1);
							}
						}

						c0 = s0;
						c1 = s1;
						c2 = s2;
						c3 = s3;
					}
				}


				void quick_accum(double d, double &e);
				void quick_prod_accum(double a, double b, double &e);

				qd_real::qd_real(double x0, double x1, double x2, double x3);
				explicit qd_real(const double *xx);

				static const qd_real _2pi;
				static const qd_real _pi;
				static const qd_real _3pi4;
				static const qd_real _pi2;
				static const qd_real _pi4;
				static const qd_real _e;
				static const qd_real _log2;
				static const qd_real _log10;
				static const qd_real _nan;
				static const qd_real _inf;

				static const double _eps;
				static const double _min_normalized;
				static const qd_real _max;
				static const qd_real _safe_max;
				static const int _ndigits;

				qd_real();
				qd_real(const char *s);
				qd_real(double d);
				qd_real(int i);

				double operator[](int i) const;
				double &operator[](int i);

				static void error(const char *msg);

				bool isnan() const;
				bool isfinite() const { return QD_ISFINITE(x[0]); }
				bool isinf() const { return QD_ISINF(x[0]); }

				static qd_real ieee_add(const qd_real &a, const qd_real &b);
				static qd_real sloppy_add(const qd_real &a, const qd_real &b);

				qd_real &operator+=(double a);
				qd_real &operator+=(const qd_real &a);

				qd_real &operator-=(double a);
				qd_real &operator-=(const qd_real &a);

				static qd_real sloppy_mul(const qd_real &a, const qd_real &b);
				static qd_real accurate_mul(const qd_real &a, const qd_real &b);

				qd_real &operator*=(double a);
				qd_real &operator*=(const qd_real &a);

				static qd_real sloppy_div(const qd_real &a, const qd_real &b);
				static qd_real accurate_div(const qd_real &a, const qd_real &b);

				qd_real &operator/=(double a);
				qd_real &operator/=(const qd_real &a);

				qd_real operator^(int n) const;

				qd_real operator-() const;

				qd_real &operator=(double a);
				qd_real &operator=(const char *s);

				bool is_zero() const;
				bool is_one() const;
				bool is_positive() const;
				bool is_negative() const;

				static qd_real rand(void);

				void to_digits(char *s, int &expn, int precision = _ndigits) const;
				void write(char *s, int len, int precision = _ndigits, bool showpos = false, bool uppercase = false) const;
				std::string to_string(int precision = _ndigits, int width = 0, std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), bool showpos = false, bool uppercase = false, char fill = ' ') const;
				static int read(const char *s, qd_real &a);

				/* Debugging methods */
				void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
				void dump_bits(const std::string &name = "",
					std::ostream &os = std::cerr) const;

				static qd_real debug_rand();

				//
				// Multiplications
				//
				qd_real n_int(const qd_real &a) {
					double x0, x1, x2, x3;

					x0 = n_int(a[0]);
					x1 = x2 = x3 = 0.0;

					if (x0 == a[0]) {
						/* First double is already an integer. */
						x1 = n_int(a[1]);

						if (x1 == a[1]) {
							/* Second double is already an integer. */
							x2 = n_int(a[2]);

							if (x2 == a[2]) {
								/* Third double is already an integer. */
								x3 = n_int(a[3]);
							}
							else {
								if (std::abs(x2 - a[2]) == 0.5 && a[3] < 0.0) {
									x2 -= 1.0;
								}
							}

						}
						else {
							if (std::abs(x1 - a[1]) == 0.5 && a[2] < 0.0) {
								x1 -= 1.0;
							}
						}

					}
					else {
						/* First double is not an integer. */
						if (std::abs(x0 - a[0]) == 0.5 && a[1] < 0.0) {
							x0 -= 1.0;
						}
					}

					renorm(x0, x1, x2, x3);
					return qd_real(x0, x1, x2, x3);
				}

				qd_real floor(const qd_real &a) {
					double x0, x1, x2, x3;
					x1 = x2 = x3 = 0.0;
					x0 = std::floor(a[0]);

					if (x0 == a[0]) {
						x1 = std::floor(a[1]);

						if (x1 == a[1]) {
							x2 = std::floor(a[2]);

							if (x2 == a[2]) {
								x3 = std::floor(a[3]);
							}
						}

						renorm(x0, x1, x2, x3);
						return qd_real(x0, x1, x2, x3);
					}

					return qd_real(x0, x1, x2, x3);
				}

				qd_real ceil(const qd_real &a) {
					double x0, x1, x2, x3;
					x1 = x2 = x3 = 0.0;
					x0 = std::ceil(a[0]);

					if (x0 == a[0]) {
						x1 = std::ceil(a[1]);

						if (x1 == a[1]) {
							x2 = std::ceil(a[2]);

							if (x2 == a[2]) {
								x3 = std::ceil(a[3]);
							}
						}

						renorm(x0, x1, x2, x3);
						return qd_real(x0, x1, x2, x3);
					}

					return qd_real(x0, x1, x2, x3);
				}

				int qd_real::read(const char *s, qd_real &qd) {
					const char *p = s;
					char ch;
					int sign = 0;
					int point = -1;  /* location of decimal point */
					int nd = 0;      /* number of digits read */
					int e = 0;       /* exponent. */
					bool done = false;
					qd_real r = 0.0;  /* number being read */

					/* Skip any leading spaces */
					while (*p == ' ') p++;

					while (!done && (ch = *p) != '\0') {
						if (ch >= '0' && ch <= '9') {
							/* It's a digit */
							int d = ch - '0';
							r *= 10.0;
							r += static_cast<double>(d);
							nd++;
						}
						else {
							/* Non-digit */
							switch (ch) {
							case '.':
								if (point >= 0)
									return -1;   /* we've already encountered a decimal point. */
								point = nd;
								break;
							case '-':
							case '+':
								if (sign != 0 || nd > 0)
									return -1;  /* we've already encountered a sign, or if its
													  not at first position. */
								sign = (ch == '-') ? -1 : 1;
								break;
							case 'E':
							case 'e':
								int nread;
								nread = std::sscanf(p + 1, "%d", &e);
								done = true;
								if (nread != 1)
									return -1;  /* read of exponent failed. */
								break;
							case ' ':
								done = true;
								break;
							default:
								return -1;

							}
						}

						p++;
					}



					/* Adjust exponent to account for decimal point */
					if (point >= 0) {
						e -= (nd - point);
					}

					/* Multiply the the exponent */
					if (e != 0) {
						r *= (qd_real(10.0) ^ e);
					}

					qd = (sign < 0) ? -r : r;
					return 0;
				}


			};

			namespace std {
				template <>
				class numeric_limits<qd_real> : public numeric_limits<double> {
				public:
					inline static double epsilon()
					{
						return qd_real::_eps;
					}
					inline static double min()
					{
						return qd_real::_min_normalized;
					}
					inline static qd_real max()
					{
						return qd_real::_max;
					}
					inline static qd_real safe_max()
					{
						return qd_real::_safe_max;
					}
					static const int digits = 209;
					static const int digits10 = 62;
				};
			}

			QD_API qd_real polyeval(const qd_real *c, int n, const qd_real &x);
			QD_API qd_real polyroot(const qd_real *c, int n, const qd_real &x0, int max_iter = 64, double thresh = 0.0);

			QD_API qd_real qdrand(void);
			QD_API qd_real sqrt(const qd_real &a);

			QD_API inline bool isnan(const qd_real &a) { return a.isnan(); }
			QD_API inline bool isfinite(const qd_real &a) { return a.isfinite(); }
			QD_API inline bool isinf(const qd_real &a) { return a.isinf(); }

			//
			// Computes  qd * d  where d is known to be a power of 2.
			// This can be done component wise.
			//
			QD_API qd_real mul_pwr2(const qd_real &qd, double d);

			QD_API qd_real operator+(const qd_real &a, const qd_real &b);
			QD_API qd_real operator+(const qd_real &a, double b);
			QD_API qd_real operator+(double a, const qd_real &b);

			QD_API qd_real operator-(const qd_real &a, const qd_real &b);
			QD_API qd_real operator-(const qd_real &a, double b);
			QD_API qd_real operator-(double a, const qd_real &b);

			QD_API qd_real operator*(const qd_real &a, const qd_real &b);
			QD_API qd_real operator*(const qd_real &a, double b);
			QD_API qd_real operator*(double a, const qd_real &b);

			QD_API qd_real operator/(const qd_real &a, const qd_real &b);
			QD_API qd_real operator/(const qd_real &a, double b);
			QD_API qd_real operator/(double a, const qd_real &b);

			QD_API qd_real sqr(const qd_real &a);
			QD_API qd_real sqrt(const qd_real &a);
			QD_API qd_real pow(const qd_real &a, int n);
			QD_API qd_real pow(const qd_real &a, const qd_real &b);
			QD_API qd_real npwr(const qd_real &a, int n);

			QD_API qd_real nroot(const qd_real &a, int n);

			QD_API qd_real rem(const qd_real &a, const qd_real &b);
			QD_API qd_real drem(const qd_real &a, const qd_real &b);
			QD_API qd_real divrem(const qd_real &a, const qd_real &b, qd_real &r);

			double  to_double(const qd_real &a);
			int     to_int(const qd_real &a);


			//
			// Equality comparison
			//
			inline bool operator==(const qd_real &a, double b) {
				return (a[0] == b && a[1] == 0.0 && a[2] == 0.0 && a[3] == 0.0);
			}

			inline bool operator==(double a, const qd_real &b) {
				return (b == a);
			}

			inline bool operator==(const qd_real &a, const qd_real &b) {
				return (a[0] == b[0] && a[1] == b[1] &&
					a[2] == b[2] && a[3] == b[3]);
			}


			//
			// Less than comparison
			//
			inline bool operator<(const qd_real &a, double b) {
				return (a[0] < b || (a[0] == b && a[1] < 0.0));
			}

			inline bool operator<(double a, const qd_real &b) {
				return (b > a);
			}

			inline bool operator<(const qd_real &a, const qd_real &b) {
				return (a[0] < b[0] ||
					(a[0] == b[0] && (a[1] < b[1] ||
					(a[1] == b[1] && (a[2] < b[2] ||
						(a[2] == b[2] && a[3] < b[3]))))));
			}


			//
			// Greater than comparison
			//
			inline bool operator>(const qd_real &a, double b) {
				return (a[0] > b || (a[0] == b && a[1] > 0.0));
			}

			inline bool operator>(double a, const qd_real &b) {
				return (b < a);
			}

			inline bool operator>(const qd_real &a, const qd_real &b) {
				return (a[0] > b[0] ||
					(a[0] == b[0] && (a[1] > b[1] ||
					(a[1] == b[1] && (a[2] > b[2] ||
						(a[2] == b[2] && a[3] > b[3]))))));
			}


			//
			// Less than or Equal to comparison
			//
			inline bool operator<=(const qd_real &a, double b) {
				return (a[0] < b || (a[0] == b && a[1] <= 0.0));
			}

			inline bool operator<=(double a, const qd_real &b) {
				return (b >= a);
			}

			inline bool operator<=(const qd_real &a, const qd_real &b) {
				return (a[0] < b[0] ||
					(a[0] == b[0] && (a[1] < b[1] ||
					(a[1] == b[1] && (a[2] < b[2] ||
						(a[2] == b[2] && a[3] <= b[3]))))));
			}


			//
			// Greater than or Equal to comparison
			//
			inline bool operator>=(const qd_real &a, double b) {
				return (a[0] > b || (a[0] == b && a[1] >= 0.0));
			}

			inline bool operator>=(double a, const qd_real &b) {
				return (b <= a);
			}

			inline bool operator>=(const qd_real &a, const qd_real &b) {
				return (a[0] > b[0] ||
					(a[0] == b[0] && (a[1] > b[1] ||
					(a[1] == b[1] && (a[2] > b[2] ||
						(a[2] == b[2] && a[3] >= b[3]))))));
			}

			QD_API bool operator!=(const qd_real &a, const qd_real &b);
			QD_API bool operator!=(double a, const qd_real &b);
			QD_API bool operator!=(const qd_real &a, double b);
			// Not equal to comparison
			inline bool operator!=(const qd_real &a, double b) {
				return !(a == b);
			}

			inline bool operator!=(double a, const qd_real &b) {
				return !(a == b);
			}

			inline bool operator!=(const qd_real &a, const dd_real &b) {
				return !(a == b);
			}

			inline bool operator!=(const dd_real &a, const qd_real &b) {
				return !(a == b);
			}

			inline bool operator!=(const qd_real &a, const qd_real &b) {
				return !(a == b);
			}


			QD_API qd_real fabs(const qd_real &a);
			QD_API qd_real abs(const qd_real &a);    /* same as fabs */

			QD_API qd_real ldexp(const qd_real &a, int n);

			QD_API qd_real n_int(const qd_real &a);
			QD_API qd_real quick_nint(const qd_real &a);
			QD_API qd_real floor(const qd_real &a);
			QD_API qd_real ceil(const qd_real &a);
			QD_API qd_real trunc_int(const qd_real &a);

			QD_API qd_real sin(const qd_real &a);
			QD_API qd_real cos(const qd_real &a);
			QD_API qd_real tan(const qd_real &a);
			QD_API void sincos(const qd_real &a, qd_real &s, qd_real &c);

			QD_API qd_real asin(const qd_real &a);
			QD_API qd_real acos(const qd_real &a);
			QD_API qd_real atan(const qd_real &a);
			QD_API qd_real atan2(const qd_real &y, const qd_real &x);

			QD_API qd_real exp(const qd_real &a);
			QD_API qd_real log(const qd_real &a);
			QD_API qd_real log10(const qd_real &a);

			QD_API qd_real sinh(const qd_real &a);
			QD_API qd_real cosh(const qd_real &a);
			QD_API qd_real tanh(const qd_real &a);
			QD_API void sincosh(const qd_real &a, qd_real &sin_qd, qd_real &cos_qd);

			QD_API qd_real asinh(const qd_real &a);
			QD_API qd_real acosh(const qd_real &a);
			QD_API qd_real atanh(const qd_real &a);

			QD_API qd_real qdrand(void);

			QD_API qd_real max(const qd_real &a, const qd_real &b);
			QD_API qd_real max(const qd_real &a, const qd_real &b, const qd_real &c);
			QD_API qd_real min(const qd_real &a, const qd_real &b);
			QD_API qd_real min(const qd_real &a, const qd_real &b, const qd_real &c);

			QD_API qd_real fmod(const qd_real &a, const qd_real &b);

			QD_API std::ostream &operator<<(std::ostream &s, const qd_real &a);
			QD_API std::istream &operator>>(std::istream &s, qd_real &a);

#endif // !_QD_QD_REAL

		} // namespace qd ends

	} // namespace multiprecision ends

} // namespace boost ends

#endif // !BOOST_MP_QD_BACKEND_HPP
}