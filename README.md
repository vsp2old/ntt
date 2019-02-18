# 数論変換による多倍長整数乗算の実装
FFT乗算は畳み込み定理を応用したものですが、実数演算が必要なため実行速度や演算精度に難点があります。しかし、整数のみでおこなうNTT（Number-Theoretical Transform）があることを知りました。理屈は分かっても具体的な実装方法がわからず、探してみたところ次の記事を見つけました。  
・[任意modでの畳み込み演算をO(n log(n))でmath314のブログ](https://math314.hateblo.jp/entry/2015/05/07/014908)  
・[数論変換(Number Theoretic Transform; NTT)|るまライブラリ](https://lumakernel.github.io/ecasdqina/math/FFT/NTT)  
これらを参考に多倍長整数の乗算を実装してみましたが、そのままではうまくいきません。原因は、ここにあるプログラムはlong long型を基本に作られているため、unsigned型を使うと破綻してしまいます。多くの多倍長整数では符号無し64ビット配列を使用しています。GMP（GNU Multi-Precision Library）では63ビット配列などを使うオプションが用意されているようですが、できれば64ビットのままで使いたいものです。  
## 1. モジュロ演算関数
まず、四則演算を剰余演算として実行するための次の関数を作ります。  
~~~
  加算  	T addm<T>(T a, T b, T mod);  
  減算  	T subm<T>(T a, T b, T mod);  
  乗算  	T mulm<T>(T a, T b, T mod);  
  除算  	T divm<T>(T a, T b, T mod);  
  逆数  	T invm<T>(T a, T mod);  
  巾乗  	T powm<T>(T a, T n, T mod);  
~~~
`T`は`unsigned int` または`unsigned long`です。  
## 2. 数論変換クラス(NTT)
数論変換を使った畳み込み演算を実行するためのclassを作ります。
~~~
class NTT<H,T>
{
    T mod, List[H], InvList[H];
  public:
    NTT(T m, T p);
    T get_mod();
    void forward(T *a, int n);
    void inverse(T *a, int n);
    int convolution(T *rp, T *ap, int an, T *bp, int bn, int n = 0);
};
~~~
法として素数`m`と原始根`p`を引数にクラスインスタンスを作ります。変換はサイズ`n`の配列`(a,n)`を引数にして`forward`を呼びます。復元は`inverse`です。内部では回転因子をあらかじめ計算しておいた`List[H]`,`InvList[H]`を参照して高速剰余変換をおこないます。ここで`H`は配列の長さですが`forward`,`inverse`の引数`n`に対して`n`≤`2^H` でなければなりません。`convolution`は配列`(ap,an)`と`(bp,bn)`の畳み込み剰余演算を行い、配列`(rp,an+bn)`に書き込みます。戻り値は引数`n`を省略すると`an+bn-1`を、`n`を指定すると`n`を返します。
NTTクラスインスタンスはNTT_TABLEに保存しておき、ガーナーアルゴリズムの適用ではこれを参照します。
## 3. ガーナーアルゴリズムの実装
NTTクラスでは法として素数を指定しなければなりませんが、次に任意の法で畳み込み剰余演算を行なう関数を作ります。
~~~
int garner_convolution<N,T>(T *rp, const T *ap, int an, const T *bp, int bn, T mod);
~~~
引数はNTTクラスの`convolution`とほぼ同じですが、最後の引数に法として任意の数値`mod`を指定します。NはNTT_TABLEからN個の異なるNTTクラスインスタンスを選ぶことを意味します。これに`mod`を加えたN+1個の連立方程式を作り、ガーナーアルゴリズムを適用してこれを解き、結果を配列`(rp,an+bn)`に書き込みます。返り値は`an+bn-1`です。
以上を実現したプログラムが`ntt.cpp`です。これは、
~~~
gcc –g –O2 ntt.cpp –std=c++11 –lstdc++ -o ntt
~~~
でコンパイルできます。
## 4. 補助演算関数
いよいよ多倍長整数乗算を実装しますが、それには次の補助演算関数が必要です。
~~~
int add<T>(T *rp, T *ap, int an, T *bp, int bn);
int mul_1<T>(T *rp, T *ap, int an, T h);
~~~
`add`は配列`(ap,an)`と`(bp,bn)`を多倍長整数として加算し`(rp,max(an,bn)+1)`に書き込みます。`mul_1`は配列`(ap,an)`と数値`h`の乗算結果を`(rp,an+1)`に書き込みます。戻り値は共に配列`rp`の有効長です。
## 5. 多倍長整数乗算の実装
多倍長整数乗算を行なうには、`ntt.cpp`の中の`garner_convolution`で、任意`mod`での剰余を求める次の部分、
~~~
constants[N][i] = addm<T>(constants[N][i], mulm<T>(coeffs[N], t, m[N]), m[N]);
~~~
に対して桁上げ処理を施せばよいわけです。（`i`は桁位置を表す。）  
まず、`mulm<T>(coeffs[N], t, m[N])`を次のように変更します。
~~~
tn = mul_1<T>(temp, coefficient, cn, t);
~~~
ここで配列`(coefficient,cn)`は`m[N]`に桁上げ処理を施したものです。これと`t`との積を一時変数`(temp,tn)`に書き込みます。さらに`constants[N]`に足し込みます。
~~~
add<T>(&constants[N][i], &constants[N][i], vn, temp, tn);
~~~
これをすべての桁について（`i=0,…,n-1`）繰り返すと求める結果が得られます。
これを実現したプログラムが`ntt1.cpp`です。ここではモジュロ演算関数をアセンブル言語で実装し`(gmath.cpp)`、高速化を行なっています。
~~~
gcc –g –O2 ntt1.cpp gmath.cpp –std=c++11 –lstdc++ -o ntt1
~~~
でコンパイルできます。
## 6. 複素モジュロ演算関数
次に、数論変換クラスを複素整数に拡張します。そのため次の複素モジュロ演算関数を作ります。
~~~
加算	void addz<T>(T& zr, T& zi, T xr, T xi, T mod);
減算	void subz<T>(T& zr, T& zi, T xr, T xi, T mod);
乗算	void mulz<T>(T& zr, T& zi, T xr, T xi, T mod);
除算	void divz<T>(T& zr, T& zi, T xr, T xi, T mod);
逆数	void invz<T>(T& zr, T& zi, T xr, T xi, T mod);
巾乗	void powz<T>(T& zr, T& zi, T xr, T xi, long n, T mod);
~~~
`addz`は複素変数`[zr,zi]`に複素数値`[xr,xi]`を足し込みます。`subz,mulz,divz`も同様です。`invz`は`[xr,xi]`の逆数を`[zr,zi]`に書き込みます。`powz`は`[xr,xi]^n`を`[zr,zi]`に書き込みます。いずれも法`mod`による剰余を求めていることに注意してください。
## 7. 複素数論変換クラス(ZNTT)
複素モジュロ演算関数を使えば複素数論変換クラスを作れます。
~~~
class ZNTT<H,T>
{
    T mod, pr, pi; long t;
  public:
    ZNTT(T m, T zr, T zi);
    T get_mod();
    void forward(T *a, int n);
    void inverse(T *a, int n);
    int convolution(T *rp, T *ap, int an, T *bp, int bn, int n = 0);
};
~~~
コンストラクタの引数は素数`m`と原始根`[zr,zi]`を指定します。
## 8. 32ビットで計算可能な複素数論変換クラス
メルセンヌ素数`M31 = 2^31-1`を法とする剰余演算は32ビットで計算できます。ガウス整数 `Z[i]` のイデアル`(M31)`による剰余は
```
Z[i]⁄(M31) = {a+bi│0 ≤ a, b < M31}
```
と表現でき、32ビットで計算可能な複素数論変換クラスを作れます。ここで原始根として `10+9i` を選ぶと配列の長さを `n ≤ 2^32` にとることが可能です。すなわち、
~~~
ZNTT<32, unsigned int>(M31, 10, 9);
~~~
としてインスタンスを作り、数論変換を行なうことができます。これを実現したプログラムが`ntt2.cpp`です。これは
~~~
gcc –g –O2 ntt2.cpp gmath.cpp –std=c++11 –lstdc++ -o ntt2
~~~
でコンパイルできます。
## 9. 64ビットで計算可能な複素数論変換クラス
メルセンヌ素数`M61 = 2^61-1`を法とする剰余演算は64ビットで計算できます。ガウス整数`Z[i]`のイデアル`(M61)` による剰余は
```
Z[i]⁄(M61) = {a+bi│0 ≤ a, b < M61}
```
と表現でき、64ビットで計算可能な複素数論変換クラスを作れます。ここで原始根として `5+2i` を選ぶと配列の長さを`n ≤ 2^62` にとることが可能です。すなわち、
~~~
ZNTT<62, unsigned long>(M61, 5, 2);
~~~
としてインスタンスを作り、数論変換を行なうことができます。これを実現したプログラムが`ntt3.cpp`です。これは
~~~
gcc –g –O2 ntt3.cpp gmath.cpp –std=c++11 –lstdc++ -o ntt3
~~~
でコンパイルできます。  
# Auther 
Shinji Hayashi (vsp2old) mailto:yatcho@hotmail.co.jp


