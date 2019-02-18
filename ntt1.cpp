#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <utility>

#define EXTERN_INLINE 1

#if EXTERN_INLINE
extern uint32_t mcn_mulm(uint32_t, uint32_t, uint32_t);
extern uint32_t mcn_addm(uint32_t, uint32_t, uint32_t);
extern uint32_t mcn_subm(uint32_t, uint32_t, uint32_t);
extern uint32_t mcn_invm(uint32_t, uint32_t, uint32_t);
extern uint64_t mcn_mulm(uint64_t, uint64_t, uint64_t);
extern uint64_t mcn_addm(uint64_t, uint64_t, uint64_t);
extern uint64_t mcn_subm(uint64_t, uint64_t, uint64_t);
extern uint64_t mcn_invm(uint64_t, uint64_t);
extern uint32_t mpn_add_n(uint32_t *, const uint32_t *, const uint32_t *, int);
extern uint32_t mpn_add_1(uint32_t *, const uint32_t *, int, uint32_t);
extern uint32_t mpn_mul_1(uint32_t *, const uint32_t *, int, uint32_t);
extern uint64_t mpn_add_n(uint64_t *, const uint64_t *, const uint64_t *, long);
extern uint64_t mpn_add_1(uint64_t *, const uint64_t *, long, uint64_t);
extern uint64_t mpn_mul_1(uint64_t *, const uint64_t *, long, uint64_t);
#endif

template <typename T> void dump(const char *, T *, int);

using std::swap;

namespace FMT {

template <typename T>
void umul_11(T& w1, T& w0, T u, T v)
{
static const T half_bits = sizeof(T) * 4;
static const T x100 = ( 1L << half_bits);
static const T x0ff = (x100 - 1);

	T ul = u & x0ff;
	T uh = u >> half_bits;
	T vl = v & x0ff;
	T vh = v >> half_bits;

	T x0 = ul * vl;
	T x1 = ul * vh;
	T x2 = uh * vl;
	T x3 = uh * vh;

	x1 += x0 >> half_bits;
	x1 += x2;
	if (x1 < x2) x3 += x100;

	w1 = x3 + (x1 >> half_bits);
	w0 = (x1 << half_bits) + (x0 & x0ff);
}

template <typename T>
T mpx_mul_1(T *rp, const T *up, long n, T vl)
{
	T ul, cl, ph, pl;

	assert (n >= 1);

	cl = 0;
	do {
		ul = *up++;
		umul_11(ph, pl, ul, vl);
		pl += cl;
		cl = ph + (pl < cl);
		*rp++ = pl;
	} while (--n != 0);

	return cl;
}

#if EXTERN_INLINE

template <> inline
uint64_t mpx_mul_1(uint64_t *f, const uint64_t *g, long n, uint64_t h)
{
	return mpn_mul_1(f, g, n, h);
}

template <> inline
uint32_t mpx_mul_1(uint32_t *f, const uint32_t *g, long n, uint32_t h)
{
	return mpn_mul_1(f, g, n, h);
}

#endif

template <class T>
int mul_1(T *rp, const T *ap, int an, T h)
{
	if (an == 0 || h == 0) {
		rp[0] = 0;
		return 0;
	}
	rp[an] = mpx_mul_1<T>(rp, ap, an, h);
	an += (rp[an] != 0);
	return an;
}

template <typename T>
T mpx_add_1(T *rp, const T *ap, long n, T b)
{
	long i = 0;

	do {
		T r = ap[i] + b;
		/* Carry out */
		b = (r < b);
		rp[i] = r;
	} while (++i < n);

	return b;
}

#if EXTERN_INLINE

template <> inline
uint32_t mpx_add_1(uint32_t *rp, const uint32_t *ap, long n, uint32_t b)
{
	return mpn_add_1(rp, ap, n, b);
}

template <> inline
uint64_t mpx_add_1(uint64_t *rp, const uint64_t *ap, long n, uint64_t b)
{
	return mpn_add_1(rp, ap, n, b);
}

#endif

template <typename T>
T mpx_add_n(T *rp, const T *ap, const T *bp, long n)
{
	T cy; long i;

	for (i = 0, cy = 0; i < n; i++)
	{
		T a, b, r;
		a = ap[i]; b = bp[i];
		r = a + cy;
		cy = (r < cy);
		r += b;
		cy += (r < b);
		rp[i] = r;
	}

	return cy;
}

#if EXTERN_INLINE

template <> inline
uint32_t mpx_add_n(uint32_t *rp, const uint32_t *ap, const uint32_t *bp, long n)
{
	return mpn_add_n(rp, ap, bp, n);
}

template <> inline
uint64_t mpx_add_n(uint64_t *rp, const uint64_t *ap, const uint64_t *bp, long n)
{
	return mpn_add_n(rp, ap, bp, n);
}

#endif

template <class T>
T mpx_add(T *rp, const T *ap, long an, const T *bp, long bn)
{
	T cy = mpx_add_n(rp, ap, bp, bn);
	if (an > bn)
		cy = mpx_add_1(rp + bn, ap + bn, an - bn, cy);
	return cy;
}

template <class T>
int add(T *rp, const T *ap, int an, const T *bp, int bn)
{
	if (an < bn) {
		swap(ap, bp);
		swap(an, bn);
	}
	if (an == 0) {
		if (bn == 0) rp[0] = 0; else
		for (int i = 0; i < bn; ++i) rp[i] = bp[i];
		return bn;
	} else if (bn == 0) {
		for (int i = 0; i < an; ++i) rp[i] = ap[i];
		return an;
	}
	rp[an] = mpx_add<T>(rp, ap, an, bp, bn);
	return an + (rp[an] != 0);
}

template <typename T>
T invm(T a, T mod)
{
	int s = 0, t = 0;
	T u = 0, v = 1, b = mod;
	while (a) {
		T q = b / a, w = q * v;
		if (u < w) {
			swap(u, w);
			s = ~t;
		}
		if (s ^ t) u += w;
		else	   u -= w;
		swap(u, v);
		swap(s, t);
		swap(b %= a, a);
	}
	if (s) u = mod - u;
	return u;
}

#if EXTERN_INLINE

template <> inline
uint32_t invm(uint32_t a, uint32_t mod)
{
	return mcn_invm(a, mod);
}

template <> inline
uint64_t invm(uint64_t a, uint64_t mod)
{
	return mcn_invm(a, mod);
}

#endif

template <class T>
T addm(T a, T b, T mod)
{
	T s = a + b;
	if (s < a) s -= mod;
	return s % mod;
}

#if EXTERN_INLINE

template <> inline
uint32_t addm(uint32_t a, uint32_t b, uint32_t mod)
{
	return mcn_addm(a, b, mod);
}

template <> inline
uint64_t addm(uint64_t a, uint64_t b, uint64_t mod)
{
	return mcn_addm(a, b, mod);
}

#endif

template <class T>
T subm(T a, T b, T mod)
{
	T s = a - b;
	if (s > a) s += mod;
	return s % mod;
}

#if EXTERN_INLINE

template <> inline
uint32_t subm(uint32_t a, uint32_t b, uint32_t mod)
{
	return mcn_subm(a, b, mod);
}

template <> inline
uint64_t subm(uint64_t a, uint64_t b, uint64_t mod)
{
	return mcn_subm(a, b, mod);
}

#endif

template <class T>
T mulm(T a, T b, T mod)
{
static const int n = 8 * sizeof(T);
static const T m = 1L << (n - 1);

	umul_11(a, b, a, b);
	a %= mod;
	if (a != 0) {
		for (int s = 1; s <= n; ++s) {
			T c = a & m;
			a = (a << 1) | ((b & m) != 0); b <<= 1;
			(a >= mod || c != 0) && (a -= mod);
		}
		return a;
	} else
		return b % mod;
}

#if EXTERN_INLINE

template <> inline
uint32_t mulm(uint32_t a, uint32_t b, uint32_t mod)
{
	return mcn_mulm(a, b, mod);
}

template <> inline
uint64_t mulm(uint64_t a, uint64_t b, uint64_t mod)
{
	return mcn_mulm(a, b, mod);
}

#endif

template <class T>
T divm(T a, T b, T mod)
{
	return mulm<T>(a, invm<T>(b, mod), mod);
}

template <class T>
T powm(T a, T n, T mod)
{
	T ret = 1;
	T p = a % mod;
	while (n) {
		if (n & 1)
			ret = mulm<T>(ret, p, mod);
		p = mulm<T>(p, p, mod);
		n >>= 1;
	}
	return ret;
}

constexpr int MAX_H = 20;

/* NTT_PRIMES */

constexpr unsigned long LTT_PRIMES[][2] = {
	{     286354059558913UL, 5}, // 2^36 * 3^2 * 463 + 1
	{     285838663483393UL, 5}, // 2^35 * 3 * 47 * 59 + 1
	{     279782759596033UL, 5}, // 2^33 * 3^2 * 7 * 11 * 47 + 1
	{ 6195545712378249217UL, 5}, // 2^48 * 3 * 11 * 23 * 29 + 1
	{ 1139410705724735489UL, 3}, // 2^52 * 11 * 23 + 1
	{     300956948365313UL, 3}, // 2^35 * 19 * 461 + 1
	{     293741403308033UL, 3}, // 2^35 * 83 * 103 + 1
	{     288570262683649UL,13}, // 2^34 * 3 * 11 * 509 + 1
	{     278442729799681UL, 7}, // 2^33 * 3 * 5 * 2161 + 1
#if EXTERN_INLINE
	{15564440312192434177UL,87}, // 2^59 * 3^3 + 1
	{18446744069414584321UL, 7}, // 2^32 * 3 * 5 * 17 * 257 * 65537 + 1
#endif
//	{        0x70000000 0000000000000001,11}, // 2^92 * 7 + 1
//	{0xA200000000000000 0000000000000001,57}, // 2^121 * 81 + 1
//	{0x6000000000000000 0000000000000000 0000000000000001, 3}, // 2^189 * 3 + 1
//	{0xD000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000001, 7}, // 2^316 * 13 + 1
};

constexpr int LTT_SIZE = sizeof(LTT_PRIMES) / sizeof(*LTT_PRIMES);

constexpr unsigned int MTT_PRIMES[][2] = {
	{3221225473, 5}, // 2^30 * 3 + 1
	{3489660929,11}, // 2^28 * 13 + 1
	{1541406721,17}, // 2^21 * 3 * 5 * 7^2 + 1,
	{1224736769, 3}, // 2^24 * 73 + 1,
	{1053818881, 7}, // 2^20 * 3 * 5 * 67 + 1
	{1051721729, 6}, // 2^20 * 17 * 59 + 1
	{1045430273, 3}, // 2^20 * 997 + 1
	{1012924417, 5}, // 2^21 * 3 * 7 * 23 + 1
	{1007681537, 3}, // 2^20 * 31^2 + 1
	{1004535809, 3}, // 2^21 * 479 + 1
	{ 998244353, 3}, // 2^23 * 7 * 17 + 1
	{ 985661441, 3}, // 2^22 * 5 * 47 + 1
	{ 976224257, 3}, // 2^20 * 7^2 * 19 + 1
	{ 975175681,17}, // 2^21 * 3 * 5 * 31 + 1
	{ 962592769, 7}, // 2^21 * 3^3 * 17 + 1
	{ 950009857, 7}, // 2^21 * 4 * 151 + 1
	{ 943718401, 7}, // 2^22 * 3^2 * 5^2 + 1
	{ 935329793, 3}, // 2^22 * 223 + 1
	{ 924844033, 5}, // 2^21 * 3^2 * 7^2 + 1
	{ 897581057, 3}, // 2^23 * 107 + 1
	{ 645922817, 3}, // 2^23 * 7 * 11 + 1
	{ 595591169, 3}, // 2^23 * 71 + 1
	{ 469762049, 3}, // 2^26 * 7 + 1
	{ 167772161, 3}, // 2^25 * 5 + 1

	{ 21626881, 67}, // 2^17 * 3 * 5 * 11 + 1
	{ 21495809,  3}, // 2^19 * 41 + 1
	{ 20316161,  3}, // 2^17 * 5 * 31 + 1
	{ 20054017,  5}, // 2^17 * 3^2 * 17 + 1
	{ 19529729,  3}, // 2^17 * 149 + 1
	{ 17006593,  5}, // 2^15 * 3 * 173 + 1
	{ 16957441, 11}, // 2^14 * 5 * 207 + 1
	{  1376257,  5}, // 2^16 * 3 * 7 + 1
};

constexpr int MTT_SIZE = sizeof(MTT_PRIMES) / sizeof(*MTT_PRIMES);

/* end of NTT_PRIMES */

template <int H, typename T> class NTT
{
	T mod, primitive_root;

	void _ntt(T *a, int n, int sign)
	{
		assert((n ^ (n & -n)) == 0); // n = 2^k
		assert(((mod - 1) % n) == 0);

		T base = powm<T>(primitive_root, (mod - 1) / n, mod);
		if (sign <= 0) base = invm<T>(base, mod);
		int Pj = n;
		while (Pj > 1) {
			int Nj = Pj >> 1;
			for (int k = 0; k < n; k += Pj) {
				T u = a[k];
				a[k] = addm<T>(u, a[k+Nj], mod);
				a[k+Nj] = subm<T>(u, a[k+Nj], mod);
			}
			T w = base;
			for (int L = 1; L < Nj; ++L) {
				for (int k = L; k < n; k += Pj) {
					T u = a[k];
					a[k] = addm<T>(u, a[k+Nj], mod);
					a[k+Nj] = mulm<T>(w, subm<T>(u, a[k+Nj], mod), mod);
				}
				w = mulm<T>(w, base, mod);
			}
			base = mulm<T>(base, base, mod);
			Pj = Nj;
		}
		// bit reverse
		for (int i = 0, j = 1; j < n - 1; ++j) {
			for (int k = n >> 1; k >(i ^= k); k >>= 1);
			if (j < i) swap(a[i], a[j]);
		}
		if (sign < 0) {
			T h = invm<T>(n, mod);
			for (int i = 0; i < n; ++i) a[i] = mulm<T>(a[i], h, mod);
		}
	}

  public:

	T get_mod() { return mod; }

	NTT(T m, T p) : mod(m), primitive_root(p)
	{
		assert((mod & ((1L << H) - 1)) == 1);
	}

	void forward(T *input, int n)
	{
		_ntt(input, n, 1);
	}

	void backward(T *input, int n)
	{
		_ntt(input, n, 0);
	}

	void inverse(T *input, int n)
	{
		_ntt(input, n, -1);
	}

	// operate convolution.
	int convolution(T *rp, const T *ap, int an, const T *bp, int bn, int n = 0)
	{
		int size = 0;
		if (n == 0) {
			size = an + bn - 1;
			n = 1;
			while (n < size) n <<= 1;
		}
		if (size == 0) size = n;

		T *_a = new T[n];
		T *_b = new T[n];
		for (int i = 0; i < n; ++i) {
			_a[i] = i < an ? ap[i] % mod : 0;
			_b[i] = i < bn ? bp[i] % mod : 0;
		}

		forward(_a, n); forward(_b, n);
		for (int i = 0; i < n; ++i)
			_a[i] = mulm<T>(_a[i], _b[i], mod);
		inverse(_a, n);

		memcpy(rp, _a, size * sizeof(T));
		delete[] _a;
		delete[] _b;
		return size;
	}
};

#define NTT_PRIMES LTT_PRIMES
#define NTT_SIZE LTT_SIZE

template <int N, typename T>
int garner_convolution(T *rp, const T *ap, int an, const T *bp, int bn)
{
	static_assert(N >= 1, "N must be positive");
	static_assert(N <= NTT_SIZE, "N is too big");

	int rn = an + bn, n = 1;
	while (n < rn) n <<= 1;

	T *workspace = new T[(N+1)*n], m[N], coeffs[N];
	T *retvalues = new T[N+n], coefficient[N] = {1};
	T *b = workspace, *constants = workspace;
	for (int I = 0; I < N; ++I) {
		m[I] = NTT_PRIMES[I][0];
		coeffs[I] = 1;
		constants += n;
		for (int j = 0; j < n; ++j) constants[j] = 0;
	}	for (int j = 0; j < n; ++j) retvalues[j] = 0;
	int cn = 1, un = 0;
	for (int I = 0; I < N; ++I)
	{
		NTT<MAX_H, T>& ntt = *((NTT<MAX_H, T>*)NTT_PRIMES[I]);
		int bc = ntt.convolution(b, ap, an, bp, bn, n); assert(bc == n);
		// Garner Algorithm
		// for each step, we solve "coeffs[I] * t[I] + constants[I] = b[I] (mod. m[I])"
		//      coeffs[I] = m[0]m[1]...m[I-1]
		//   constants[I] = t[0] + t[1]m[0] + ... + t[I-1]m[0]m[1]...m[I-2]
		for (int i = 0; i < rn; ++i) {
			constants = b + n;
			// coeffs[I] * t + constants[i] == b[i] (mod m[I]) を解く
			T t = divm<T>(subm<T>(b[i], constants[i], m[I]), coeffs[I], m[I]);
			for (int k = I + 1; k < N; k++) {
				constants += n;
			// (constants[i] += coeffs[k] * t) %= m[k];
				constants[i] = addm<T>(constants[i], mulm<T>(coeffs[k], t, m[k]), m[k]);
			}
			if (t != 0) {
				T temp[cn+1]; int vn = un > i ? un - i : 0;
				un = i + add<T>(&retvalues[i], &retvalues[i], vn, temp,
						 mul_1<T>(temp, coefficient, cn, t));
			}
		}
		if ((b += n) < constants) {
			for (int k = I + 1; k < N; k++) {
			// (coeffs[k] *= m[I]) %= m[k];
				coeffs[k] = mulm<T>(coeffs[k], m[I], m[k]);
			}
			// coefficient *= m[I];
			cn = mul_1<T>(coefficient, coefficient, cn, m[I]);
		}
	}
	while (rn > 0 && retvalues[rn-1] == 0) --rn;
	for (int i = 0; i < rn; ++i) rp[i] = retvalues[i];
	delete[] retvalues;
	delete[] workspace;
	return rn;
}

}	// end of FMT

#include <cstdio>
#include <iostream>
#include <iomanip>

using namespace std;

using FMT::add;
using FMT::mul_1;
using FMT::NTT;
using FMT::garner_convolution;

template <typename T>
void dump(const char *fmt, T *c, int n, const char *end)
{
	for (int i = 0; i < n; ++i) {
		printf(fmt, c[i]);
		if (i < n-1) printf(", ");
	}
	if (end) printf("%s", end);
}

void ntt_test()
{
	typedef unsigned int Type;
	Type v[8];
	for (int i = 0; i < 8; ++i) v[i] = (10 + i);

	Type v2[8];
	memcpy(v2, v, 8 * sizeof(Type));
	NTT<12, Type>& ntt = *((NTT<12, Type>*)FMT::MTT_PRIMES[0]);

	ntt.forward(v2, 8);
	ntt.inverse(v2, 8);

	int j;
	for (j = 0; j < 8; ++j) if (v[j] != v2[j]) break;
	assert(j == 8);

	int n = j;
	dump("%d", v2, n, "\n");
}

void convolution_test1()
{
	typedef unsigned int Type;

	Type v[4] = { 82, 91, 58, 31 };
	Type u[4] = { 13, 79, 17, 54 };
	const int N = 4, M = 4;
	Type vu[N+M]; int n, s = FMT::MTT_SIZE;
	for (int k = 0; k < s; ++k) {
		NTT<20, Type>& ntt = *((NTT<20, Type>*)FMT::MTT_PRIMES[k]);
	
		n = ntt.convolution(vu, v, N, u, M);
	//	vu = [1066, 7661, 9337, 10960, 8349, 3659, 1674];
		assert(n == N + M - 1);
		int i;
		for (i = 0; i < n; ++i) {
			Type w = 0;
			for (int j = 0; j <= i; ++j)
				if (j < N && i-j < M) w += v[j] * u[i-j];
			if (vu[i] != w) break;
		}
		assert(i == n);
	}
	printf("["); dump("%u", vu, n, "] = ");
	long su = vu[n-1];
	for (int i = n-1; i > 0; --i) {
		su = su * 100 + vu[i-1];
	}	printf("%ld\n", su);
//	su = 17 11 43 59 54 13 71 66;
}

void convolution_test2()
{
	typedef unsigned long Type;

	Type v[4] = { 87, 26, 73, 35 };
	Type u[4] = { 81, 90, 43, 56 };
	const int N = 4, M = 4;
	Type vu[N+M]; int n, s = FMT::LTT_SIZE;
	for (int k = 0; k < s; ++k) {
		NTT<32, Type>& ntt = *((NTT<32, Type>*)FMT::LTT_PRIMES[k]);
	
		n = ntt.convolution(vu, v, N, u, M);
	//	vu = [1960, 5593, 7745, 15395, 11994, 9936, 7047];
		assert(n == N + M - 1);
		int i;
		for (i = 0; i < n; ++i) {
			Type w = 0;
			for (int j = 0; j <= i; ++j)
				if (j < N && i-j < M) w += v[j] * u[i-j];
			if (vu[i] != w) break;
		}
		assert(i == n);
	}
	printf("["); dump("%lu", vu, n, "]");
	long su = vu[n-1];
	for (int i = n-1; i > 0; --i) {
		su = su * 100 + vu[i-1];
	}	printf(" = %ld\n", su);
//	su = 20 16 72 00 15 94 06 47;
}

void garner_convolution_test()
{
	typedef unsigned long Type;

	const int N = 7, M = 2;
	Type x[N] = {0xf000000000000000,0xffffffffffffffff,0xffffffffffff,0xfe00000000000000,0xffffffffffffffff,0xffffffffffffffff,0x3fff};
	Type y[M] = {0x3ffffffffffc0000,0xffffffffffffffff};
	Type z[N+M];
	int i, wn = 0, L = garner_convolution<3, Type>(z, x, N, y, M);
	Type w[L+1];
	for (i = 0; i < L; ++i) {
		for (int j = 0; j <= i; ++j) {
			Type temp[2]; int vn = wn > i ? wn - i : 0;
			if (j < N && i-j < M)
				wn = i + add<Type>(&w[i], &w[i], vn, temp,
						 mul_1<Type>(temp, &x[j], 1, y[i-j]));
		}
		if (z[i] != w[i]) break;
	}
	assert(i == L);
	printf("0x%lx", z[L-1]);
	for (i = L-1; i > 0; --i) {
		printf("%016lx", z[i-1]);
	}	printf("\n");
}

int main()
{
	ntt_test();
	convolution_test1();
	convolution_test2();
	garner_convolution_test();
	return 0;
}
