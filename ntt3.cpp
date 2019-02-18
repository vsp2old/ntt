#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <utility>
#include <cstdio>

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
int mul_1(T *f, const T *g, int n, T h)
{
	if (n == 0 || h == 0) {
		f[0] = 0;
		return 0;
	}
	f[n] = mpx_mul_1<T>(f, g, n, h);
	n += (f[n] != 0);
	return n;
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
int add(T *f, const T *g, int n, const T *h, int m)
{
	if (n < m) {
		swap(g, h);
		swap(n, m);
	}
	if (n == 0) {
		if (m == 0) f[0] = 0; else
		for (int i = 0; i < m; ++i) f[i] = h[i];
		return m;
	} else if (m == 0) {
		for (int i = 0; i < n; ++i) f[i] = g[i];
		return n;
	}
	f[n] = mpx_add<T>(f, g, n, h, m);
	return n + (f[n] != 0);
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
T powm(T a, long n, T mod)
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

template <class T>
void addz(T& zr, T& zi, T xr, T xi, T mod)
{
	zr = addm<T>(zr, xr, mod);
	zi = addm<T>(zi, xi, mod);
}

template <class T>
void subz(T& zr, T& zi, T xr, T xi, T mod)
{
	zr = subm<T>(zr, xr, mod);
	zi = subm<T>(zi, xi, mod);
}

template <class T>
void mulz(T& zr, T& zi, T xr, T xi, T mod)
{
	T u = zr;
	zr = subm<T>(mulm<T>(u, xr, mod), mulm<T>(zi, xi, mod), mod);
	zi = addm<T>(mulm<T>(u, xi, mod), mulm<T>(zi, xr, mod), mod);
}

template <class T>
void divz(T& zr, T& zi, T xr, T xi, T mod)
{
	T u = zr, v = invm<T>(addm<T>(mulm<T>(xr, xr, mod), mulm<T>(xi, xi, mod), mod), mod);
	zr = mulm<T>(addm<T>(mulm<T>(u, xr, mod), mulm<T>(zi, xi, mod), mod), v, mod);
	zi = mulm<T>(subm<T>(mulm<T>(u, xi, mod), mulm<T>(zi, xr, mod), mod), v, mod);
}

template <class T>
void invz(T& zr, T& zi, T xr, T xi, T mod)
{
	T u = invm<T>(addm<T>(mulm<T>(xr, xr, mod), mulm<T>(xi, xi, mod), mod), mod);
	zr = mulm<T>(xr, u, mod);
	zi = mulm<T>(subm<T>(mod, xi, mod), u, mod);
}

template <class T>
void powz(T& zr, T& zi, T xr, T xi, long n, T mod)
{
	zr = 1, zi = 0;
	xr %= mod, xi %= mod;
	while (n) {
		if (n & 1)
			mulz<T>(zr, zi, xr, xi, mod);
		mulz<T>(xr, xi, xr, xi, mod);
		n >>= 1;
	}
}

template <typename T> class NTT_BASE
{
  public:

	virtual ~NTT_BASE(){};
	virtual T get_mod() = 0;
	virtual int convolution(T*, const T*, int, const T*, int, int) = 0;
};

template <int H, typename T> class NTT : public NTT_BASE<T>
{
	T mod, primitive_root; long t;

	void _ntt(T *a, int n, int sign)
	{
		int h = __builtin_ctz(n);
		assert((n ^ (n & -n)) == 0); // n = 2^k
		assert( H >= h );

		// primitive n-th root of unity.
		T base = powm<T>(primitive_root, t << (H - h), mod);
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
		t = (mod - 1) >> H;
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

template <int H, typename T> class ZNTT : public NTT_BASE<T>
{
	T mod, pr, pi; long t;

	void _ntt(T *a, int n, int sign)
	{
		int N = n >> 1, h = __builtin_ctz(N);
		assert((N ^ (N & -N)) == 0); // N = 2^k
		assert( H >= h );

		T wc = pr, ws = pi; // primitive N-th root of unity.
		for (int i = 0; i < H - h; ++i)
			mulz<T>(wc, ws, wc, ws, mod);
		powz<T>(wc, ws, wc, ws, t, mod);
		if (sign <= 0) invz<T>(wc, ws, wc, ws, mod);

		int Pj = N;
		while (Pj > 1) {
			int Nj = Pj >> 1;
			for (int k = 0; k < N; k += Pj) {
				int i = k << 1, j = (k + Nj) << 1;
				T ur = a[i], ui = a[i+1];
				addz<T>(a[i], a[i+1], a[j], a[j+1], mod);
				swap(a[j], ur), swap(a[j+1], ui);
				subz<T>(a[j], a[j+1], ur, ui, mod);
			}
			T wr = wc, wi = ws;
			for (int L = 1; L < Nj; ++L) {
				for (int k = L; k < N; k += Pj) {
					int i = k << 1, j = (k + Nj) << 1;
					T ur = a[i], ui = a[i+1];
					addz<T>(a[i], a[i+1], a[j], a[j+1], mod);
					swap(a[j], ur), swap(a[j+1], ui);
					subz<T>(a[j], a[j+1], ur, ui, mod);
					mulz<T>(a[j], a[j+1], wr, wi, mod);
				}
				mulz<T>(wr, wi, wc, ws, mod);
			}
			mulz<T>(wc, ws, wc, ws, mod);
			Pj = Nj;
		}
		// bit reverse
		for (int i = 0, j = 1; j < N - 1; ++j) {
			for (int k = N >> 1; k >(i ^= k); k >>= 1);
			if (j < i) {
				swap(a[2*i],   a[2*j]  );
				swap(a[2*i+1], a[2*j+1]);
			}
		}
		if (sign < 0) {
			T h = invm<T>(N, mod);
			for (int i = 0; i < n; ++i) a[i] = mulm<T>(a[i], h, mod);
		}
	}

  public:

	T get_mod() { return mod; }

	ZNTT(T m, T qr, T qi) : mod(m), pr(qr), pi(qi), t((m - 1) >> 1) {}

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
	int convolution(T *rp, const T *ap, int an, const T *bp, int bn, int N = 0)
	{
		static const T _c = invm<T>(2, mod);
		int size = 0;
		if (N == 0) {
			size = an + bn - 1;
			N = 1;
			while (N < size) N <<= 1;
		}
		if (size == 0) size = N;

		int n = N << 1;
		T *w = new T[n];
		for (int k = 0; k < N; ++k) {
			int i = k << 1, j = i + 1;
			w[i] = k < an ? ap[k] % mod : 0;
			w[j] = k < bn ? bp[k] % mod : 0;
		}

		forward(w, n);

		w[0] = mulm<T>(w[0], w[1], mod); w[1] = 0;
		for (int i = 2; i < N; i += 2) {
			int j = n - i;
			T rep = mulm<T>(addm<T>(w[ i ], w[ j ], mod), _c, mod);
			T rem = mulm<T>(subm<T>(w[ i ], w[ j ], mod), _c, mod);
			T aip = mulm<T>(addm<T>(w[i+1], w[j+1], mod), _c, mod);
			T aim = mulm<T>(subm<T>(w[i+1], w[j+1], mod), _c, mod);
			w[i] = rep; w[i+1] = aim;
			w[j] = rep; w[j+1] = mod - aim;
			mulz<T>(w[i], w[i+1], aip, mod - rem, mod);
			mulz<T>(w[j], w[j+1], aip, rem, mod);
		}
		w[N] = mulm<T>(w[N], w[N+1], mod); w[N+1] = 0;

		inverse(w, n);

		for (int i = 0; i < size; ++i) rp[i] = w[i << 1];

		delete[] w;
		return size;
	}
};

template <int N, typename T>
int garner_convolution(NTT_BASE<T> *NTT_TABLE[N], T *rp, const T *ap, int an, const T *bp, int bn)
{
	static_assert(N > 0, "N must be positive");

	int rn = an + bn, n = 1;
	while (n < rn) n <<= 1;

	T *workspace = new T[(N+1)*n], m[N], coeffs[N];
	T *retvalues = new T[N+n], coefficient[N] = {1};
	T *b = workspace, *constants = workspace;
	for (int I = 0; I < N; ++I) {
		m[I] = NTT_TABLE[I]->get_mod();
		coeffs[I] = 1;
		constants += n;
		for (int j = 0; j < n; ++j) constants[j] = 0;
	}	for (int j = 0; j < n; ++j) retvalues[j] = 0;
	int cn = 1, un = 0;
	for (int I = 0; I < N; ++I)
	{
		int bc = NTT_TABLE[I]->convolution(b, ap, an, bp, bn, n); assert(bc == n);
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

using FMT::addm;
using FMT::mulm;
using FMT::mulz;
using FMT::powz;
using FMT::add;
using FMT::mul_1;
using FMT::NTT;
using FMT::ZNTT;
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

typedef unsigned long Type;

const int NTT_SIZE = 3;

FMT::NTT_BASE<Type> *NTT_TABLE[NTT_SIZE];

void ntt_test()
{
	static const int N = 8, n = N << 1;
	Type v[n];
	for (int i = 0; i < N; ++i) {
		v[2*i] = (10 + i);
		v[2*i+1] = (20 + i);
	}
	Type v2[n];
	memcpy(v2, v, n * sizeof(Type));

	Type M61 = (1L << 61) - 1;
	// primitive root = 5 + 2i;
	// primitive n-th root = (5 + 2i)^t; t = 2^(62 - n)*(2^60 - 1);
	ZNTT<62, Type> ntt(M61, 5, 2);

	ntt.forward(v2, n);
	ntt.inverse(v2, n);

	int j;
	for (j = 0; j < n; ++j) if (v[j] != v2[j]) break;
	assert(j == n);

	dump("%ld", v2, n, "\n");
}

void convolution_test1()
{
	const int N = 4, M = 4;
	Type v[N] = { 82, 91, 58, 31 };
	Type u[M] = { 13, 79, 17, 54 };
	Type vu[N+M]; int n;
	Type M61 = (1L << 61) - 1;
	ZNTT<62, Type>* ntt = new ZNTT<62, Type>(M61, 5, 2);
	NTT_TABLE[2] = ntt;

	n = ntt->convolution(vu, v, N, u, M);
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
	printf("["); dump("%lu", vu, n, "] = ");
	long su = vu[n-1];
	for (int i = n-1; i > 0; --i) {
		su = su * 100 + vu[i-1];
	}	printf("%ld\n", su);
//	su = 17 11 43 59 54 13 71 66;
}

void convolution_test2()
{
	Type v[4] = { 87, 26, 73, 35 };
	Type u[4] = { 81, 90, 43, 56 };
	const int N = 4, M = 4;
	Type vu[N+M]; int n;
	Type T5 = (((1L << 32) - 1) << 32) + 1, p = 7;
	// primitive root = 7;
	// primitivr n-th root = 7^t; t = 2^(32 - n)*(2^32 - 1) = (T5 - 1) / 2^n;
	NTT<32, Type>* ntt = new NTT<32, Type>(T5, p);
	NTT_TABLE[1] = ntt;
	
	n = ntt->convolution(vu, v, N, u, M);
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
	printf("["); dump("%lu", vu, n, "]");
	long su = vu[n-1];
	for (int i = n-1; i > 0; --i) {
		su = su * 100 + vu[i-1];
	}	printf(" = %ld\n", su);
//	su = 20 16 72 00 15 94 06 47;
}

void garner_convolution_test()
{
	NTT_TABLE[0] = new NTT<20, Type>((1L << 34) * 3 * 11 * 509 + 1, 13);

	const int N = 7, M = 2;
	Type x[N] = {0xf000000000000000,0xffffffffffffffff,0xffffffffffff,0xfe00000000000000,0xffffffffffffffff,0xffffffffffffffff,0x3fff};
	Type y[M] = {0x3ffffffffffc0000,0xffffffffffffffff};
	Type z[N+M];

	int i, wn = 0, L = garner_convolution<NTT_SIZE, Type>(NTT_TABLE, z, x, N, y, M);

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

void M61_test()
{
	const int N = 62, M = 12;
	Type M61 = (1L << 61) - 1;
	Type zr[N], zi[N], mr[N], mi[N];
	// N = 50 : (1000, 643)
	// N = 51 : (81, 46)
	// N = 52 : (171, 124)
	// N = 53 : (57, 32)
	// N = 54 : (195, 118)
	// N = 55 : (140, 109)
	// N = 56 : (43, 12)
	// N = 57 : (35, 8)
	// N = 58 : (17, 2)
	// N = 59 : (18, 5)
	// N = 60 : (9, 4)
	// N = 61 : (2, 1)
	// N = 62 : (5, 2)
	int L = 0;
	for (Type wc = 2; wc <= M; ++wc) {
		for (Type ws = (wc % 2) ? 2 : 1; ws < wc; ws += 2) {
			powz<Type>(zr[0], zi[0], wc, ws, (M61 - 1), M61);
			for (int i = 1; i < N; ++i) {
				zr[i] = zr[i-1], zi[i] = zi[i-1];
				mulz<Type>(zr[i], zi[i], zr[i-1], zi[i-1], M61);
			}
			if (zr[N-1] == 1 && zi[N-1] == 0)
				if (zr[N-2] == (M61 - 1) && zi[N-2] == 0)
					if (zr[N-3] == 0 && zi[N-3] == (M61 - 1)) {
						for (int i = 0; i < N; ++i)
							mr[i] = zr[i], mi[i] = zi[i];
						printf("(%lu, %lu)", wc, ws);
						if (ws < M-1) printf(",");
						++L;
					}
		}
	}	printf("\n");
	if (L > 0) for (int i = 0; i < N; ++i) {
		Type norm = addm<Type>(mulm<Type>(mr[i], mr[i], M61), mulm<Type>(mi[i], mi[i], M61), M61);
		printf("Z(%d) = |(%lu, %lu)| = %lu\n", N - 1 - i, mr[i], mi[i], norm);
	}
}

int main()
{
	ntt_test();
	convolution_test1();
	convolution_test2();
	garner_convolution_test();
//	M61_test();
	for (int i = 0; i < NTT_SIZE; ++i)
		delete NTT_TABLE[i];
	return 0;
}
