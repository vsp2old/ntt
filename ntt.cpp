#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <utility>

using std::swap;

namespace FMT {

typedef unsigned int Type;

template <typename T>
T invm(T a, T mod)
{
	int s = 0, t = 0;
	T u = 0, v = 1, b = mod;
	while (a) {
		T q = b / a, w = q * v;
		if (u < w) {
			swap(u, w);
			s ^= (s == t);
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

template <class T>
T addm(T a, T b, T mod)
{
	T s = a + b;
	if (s < a) s -= mod;
	return s % mod;
}

template <class T>
T subm(T a, T b, T mod)
{
	T s = a - b;
	if (s > a) s += mod;
	return s % mod;
}

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

template <int H, typename T> class NTT
{
	T mod, List[H], InvList[H];

	void _ntt(T *a, int n, int sign)
	{
		assert((n ^ (n & -n)) == 0); // n = 2^k
		assert((mod - 1) % n == 0);
		// bit reverse
		for (int i = 0, j = 1; j < n - 1; ++j) {
			for (int k = n >> 1; k >(i ^= k); k >>= 1);
			if (j < i) swap(a[i], a[j]);
		}
		int Nj = 1, q = 0;
		while (Nj < n) {
			int Pj = Nj << 1;
			T base = (sign <= 0) ? InvList[q] : List[q];
			for (int k = 0; k < n; k += Pj) {
				T u = a[k];
				T d = a[k+Nj];
				a[k] = addm<T>(u, d, mod);
				a[k+Nj] = subm<T>(u, d, mod);
			}
			T w = base;
			for (int L = 1; L < Nj; ++L) {
				for (int k = L; k < n; k += Pj) {
					T u = a[k];
					T d = mulm<T>(a[k+Nj], w, mod);
					a[k] = addm<T>(u, d, mod);
					a[k+Nj] = subm<T>(u, d, mod);
				}
				w = mulm<T>(w, base, mod);
			}
			++q;
			Nj = Pj;
		}
		if (sign < 0) {
			T h = invm<T>(n, mod);
			for (int i = 0; i < n; ++i)
				a[i] = mulm<T>(a[i], h, mod);
		}
	}

  public:

	T get_mod() { return mod; }

	NTT(T m, T primitive_root) : mod(m)
	{
		assert((mod & ((1 << H) - 1)) == 1);
		   List[H - 1] = powm<T>(primitive_root, (mod - 1) / (1 << H), mod);
		InvList[H - 1] = invm<T>(List[H - 1], mod);
		for(int i = H - 1; i > 0; --i) {
			   List[i - 1] = mulm<T>(   List[i],    List[i], mod);
			InvList[i - 1] = mulm<T>(InvList[i], InvList[i], mod);
		}
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

	// operate convolution
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
			_a[i] = i < an ? ap[i] : 0;
			_b[i] = i < bn ? bp[i] : 0;
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

//	NTT_OBJECTS

NTT<MAX_H, Type> NTT_TABLE[23] = {
	{3221225473, 5}, // 2^30 * 3 + 1
	{3489660929, 3}, // 2^28 * 13 + 1
	{ 469762049, 3}, // 2^26 * 7 + 1
	{ 167772161, 3}, // 2^25 * 5 + 1
	{1224736769, 3}, // 2^24 * 73 + 1,
	{ 998244353, 3}, // 2^23 * 7 * 17 + 1
	{ 897581057, 3}, // 2^23 * 107 + 1
	{ 645922817, 3}, // 2^23 * 7 * 11 + 1
	{ 595591169, 3}, // 2^23 * 71 + 1
	{ 985661441, 3}, // 2^22 * 5 * 47 + 1
	{ 935329793, 3}, // 2^22 * 223 + 1
	{ 943718401, 7}, // 2^22 * 3^2 * 5^2 + 1
	{1004535809, 3}, // 2^21 * 479 + 1
	{1012924417, 5}, // 2^21 * 3 * 7 * 23 + 1
	{ 924844033, 5}, // 2^21 * 3^2 * 7^2 + 1
	{ 962592769, 7}, // 2^21 * 3^3 * 17 + 1
	{ 950009857, 7}, // 2^21 * 4 * 151 + 1
	{ 975175681,17}, // 2^21 * 3 * 5 * 31 + 1
	{1045430273, 3}, // 2^20 * 997 + 1
	{1007681537, 3}, // 2^20 * 31^2 + 1
	{ 976224257, 3}, // 2^20 * 7^2 * 19 + 1
	{1051721729, 6}, // 2^20 * 17 * 59 + 1
	{1053818881, 7}  // 2^20 * 3 * 5 * 67 + 1
};

template <int N, typename T>
//	convolution in any mod.
int garner_convolution(T *rp, const T *ap, int an, const T *bp, int bn, T mod)
{
	static_assert(N >= 1, "N must be positive");
	static_assert(N <= sizeof(NTT_TABLE) / sizeof(*NTT_TABLE), "N is too big");
	int rn = an + bn - 1, n = 1;
	while (n < rn) n <<= 1;

	T m[N+1], coeffs[N+1], *workspace = new T[(N+2)*n];
	T *b = workspace, *constants = workspace;
	for (int I = 0; I <= N; ++I) {
		m[I] = (I < N) ? NTT_TABLE[I].get_mod() : mod;
		coeffs[I] = 1;
		constants += n;
		for (int j = 0; j < n; ++j) constants[j] = 0;
	}
	for (int I = 0; I < N; ++I)
	{
#ifndef NDEBUG
		int bc = 
#endif
		NTT_TABLE[I].convolution(b, ap, an, bp, bn, n); assert(bc == n);
		// Garner Algorithm
		// for each step, we solve "coeffs[I] * t[I] + constants[I] = b[I] (mod. m[I])"
		//      coeffs[I] = m[0]m[1]...m[I-1]
		//   constants[I] = t[0] + t[1]m[0] + ... + t[I-1]m[0]m[1]...m[I-2]
		for (int i = 0; i < n; ++i) {
			constants = b + n;
			// solve 'coeffs[I] * t + constants[i] == b[i] (mod. m[I])'
			T t = divm<T>(subm<T>(b[i], constants[i], m[I]), coeffs[I], m[I]);
			for (int k = I + 1; k <= N; k++) {
				constants += n;
			// (constants[i] += coeffs[k] * t) %= m[k];
				constants[i] = addm<T>(constants[i], mulm<T>(coeffs[k], t, m[k]), m[k]);
			}
		}
		if ((b += n) < constants) {
			for (int k = I + 1; k <= N; k++) {
			// (coeffs[k] *= m[I]) %= m[k];
				coeffs[k] = mulm<T>(coeffs[k], m[I], m[k]);
			}
		}
	}
	for (int i = 0; i < rn; ++i) rp[i] = constants[i];
	delete[] workspace;
	return rn;
}

}	// end of FMT

#include <cstdio>
#include <iostream>
#include <iomanip>

using namespace std;

using Type = FMT::Type;
using FMT::addm;
using FMT::mulm;

void ntt_test()
{
	const int N = 16;
	Type v[N];
	for (int i = 0; i < N; ++i) v[i] = (10 + i);

	Type v2[N];
	memcpy(v2, v, N * sizeof(Type));
	FMT::NTT<N, Type> ntt(3221225473, 5);

	ntt.forward(v2, N);
	ntt.inverse(v2, N);

	int j;
	for (j = 0; j < N; ++j) if (v[j] != v2[j]) break;
	assert(j == N);

	for (int i = 0; i < N; ++i) {
		printf("%u", v2[i]);
		if (i < N-1) printf(", ");
	}	printf("\n");
}

void convolution_test()
{
	Type v[4] = { 82, 91, 58, 31 };
	Type u[4] = { 13, 79, 17, 54 };
	const int N = 4, M = 4;
	Type vu[N+M]; int n;

	for (int I = 0; I < 22; ++I) {
		n = FMT::NTT_TABLE[I].convolution(vu, v, N, u, M);
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
	printf("[");
	for (int i = 0; i < n; ++i) {
		printf("%u", vu[i]);
		if (i < n-1) printf(", ");
	}	printf("]");
	long su = vu[n-1];
	for (int i = n-1; i > 0; --i) {
		su = su * 100 + vu[i-1];
	}	printf(" = %ld\n", su);
//	su = 17 11 43 59 54 13 71 66;
}

void garner_convolution_test()
{
	const int N = 10, M = 10;
	Type x[20], y[20], z[40];
	for (int i = 0; i < 20; ++i) {
		x[i] = i < N ? (1e8 + i) : 0;
		y[i] = i < M ? (1e8 + i) : 0;
	}
	Type mod = 1e9+7;
	int i, zn = FMT::garner_convolution<3, Type>(z, x, N, y, M, mod);
	for (i = 0; i < zn; ++i) {
		Type w = 0;
		for (int j = 0; j <= i; ++j)
			if (j < N && i-j < M) w = addm<Type>(w, mulm<Type>(x[j], y[i-j], mod), mod);
		if (z[i] != w) break;
	}
	assert(i == zn);
	int n = i;
	for (int i = 0; i < n; ++i) {
		printf("%u", z[i]);
		if (i < n-1) printf(", ");
	}	printf("\n");
}

int main()
{
	ntt_test();
	convolution_test();
	garner_convolution_test();
	return 0;
}
