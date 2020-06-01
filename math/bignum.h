/* 
 * Big Number Class in a single file
 * 
 * Support +, -, *, /, %
 * And +=, -=, *=, /=, %=
 * Overload operator<< and operator>> for IO
 * Use fft(Fast Fourier Transform) algorithm to do multiplications!
 * 
 * Welcome to use it !
 * Code by frezcirno
 */
#include <cmath>
#include <cstring>
#include <vector>
#include <iomanip>
#include <assert.h>

class bignum {
    using sign_type = int;
    using unit_type = long long;
    using base_type = std::vector<unit_type>;
    using digitcount_t = int;

public:
    static constexpr double PI = 3.141592653589793;
    static const digitcount_t unit_base_exp = 4;
    static const unit_type unit_base = 10000;

    struct complex {
        double re, im;
        complex(): re(0), im(0){}
        complex(double _x, double _y) : re(_x), im(_y){}
        complex operator+(const complex& b) const { return complex(re + b.re, im + b.im); }
        complex operator-(const complex& b) const { return complex(re - b.re, im - b.im); }
        complex operator*(const complex& b) const { return complex(re * b.re - im * b.im, re * b.im + im * b.re); }
        complex& operator+=(const complex& b) { re += b.re; im += b.im; return *this; }
        complex& operator-=(const complex& b) { re -= b.re; im -= b.im; return *this; }
        complex& operator*=(const complex& b) { return *this = *this * b; }
    };

private:
    base_type _base;
    sign_type _sign = 0; // 0-positive, 1-negative
    static std::vector<bignum> powTwo;

public:
    bignum(long long x = 0) : _base() {
        unsigned long long ux = x;
        if (x < 0) { ux = -x; _sign = 1; }
        while (ux) { _base.push_back(ux % bignum::unit_base); ux /= bignum::unit_base; }
    }
    bignum(const char* s) : _base() { assign(s); }
    bignum(const bignum& x) : _base(x._base), _sign(x._sign) { }
    bignum& operator=(const bignum& x)
    {
        _sign = x._sign;
        _base = x._base;
        return *this;
    }
    void clear() { _base.clear(); _sign = 0; }
    bignum& assign(const char* s){
        int i;
        long long idx = 1;
        clear();
        for (i = strlen(s) - 1; i >= 0 && isdigit(s[i]); i--, idx *= 10) {
            if (idx >= bignum::unit_base) idx = 1;
            if (idx == 1) _base.push_back(s[i] - '0');
            else _base.back() += idx * (s[i] - '0');
        }
        if (s[i] == '-') _sign = 1;
        trim();
        return *this;
    }
    const base_type& base() const { return _base; }
    bool sign() const { return _sign; }
    operator bool() const { return !_base.empty(); }
    friend std::istream &operator>>(std::istream &is, bignum &b){
        std::string s;
        is >> s;
        b.assign(s.data());
        return is;
    }
    friend std::ostream &operator<<(std::ostream &os, const bignum &b){
        if(b._sign) os << '-';
        int first = 1;
        for (int i = b._base.size() - 1; i >= 0; i--)
            if (first && first--) 
                os << b._base[i];
            else
                os << std::setfill('0') << std::setw(bignum::unit_base_exp) << b._base[i];
        if (first) os << '0';
        return os;
    }
    bignum& operator+=(long long x) { return *this = *this + x; }
    bignum& operator+=(const bignum& b) { return *this = *this + b; }
    bignum& operator-=(long long x) { return *this = *this - x; }
    bignum& operator-=(const bignum& b) { return *this = *this - b; }
    bignum &operator*=(long long x) { return *this = *this * x; }
    bignum &operator*=(const bignum& b) { return *this = *this * b; }
    bignum &operator/=(long long x) { return *this = *this / x; }
    bignum &operator/=(const bignum& b) { return *this = *this / b; }
    bignum &operator%=(long long x) { return *this = *this % x; }
    bignum &operator%=(const bignum& b) { return *this = *this % b; }
    bignum operator-() const { return bignum(_base, !_sign); }
    digitcount_t digitCount() const {
        if(_base.empty()) return 0;
        digitcount_t res = bignum::unit_base_exp * (_base.size() - 1);
        unit_type back = _base.back();
        if (back <= 10) return res + 1;
        if (back <= 100) return res + 2;
        if (back <= 1000) return res + 3;
        if (back <= 10000) return res + 4;
        if (back <= 100000) return res + 5;
        if (back <= 1000000) return res + 6;
        if (back <= 10000000) return res + 7;
        if (back <= 100000000) return res + 8;
        if (back <= 1000000000) return res + 9;
        if (back <= 10000000000) return res + 10;
        if (back <= 100000000000) return res + 11;
        if (back <= 1000000000000) return res + 12;
        if (back <= 10000000000000) return res + 13;
        if (back <= 100000000000000) return res + 14;
        if (back <= 1000000000000000) return res + 15;
        if (back <= 10000000000000000) return res + 16;
        if (back <= 100000000000000000) return res + 17;
        if (back <= 1000000000000000000) return res + 18;
        return res + 19;
    }
    unit_type getdigit(digitcount_t idx) const {
        digitcount_t seg = idx / bignum::unit_base_exp;
        if (seg >= _base.size()) return 0;
        unit_type unit = _base[seg];
        digitcount_t inseg = idx % bignum::unit_base_exp;
        while (inseg--) unit /= 10;
        return unit % 10;
    }

private:
    bignum(const base_type& base, sign_type sign) : _base(base), _sign(sign) {}
    void trim()
    {
        while(_base.size() && _base.back() == 0)
            _base.pop_back();
        if (_base.empty())
            _sign = 0;
    }
    static bignum absadd(const bignum& b1, const bignum& b2) {
        bignum res(b1._base, 0);
        int l1 = b1._base.size(), l2 = b2._base.size(), carry = 0;
        if (l2 > res._base.size()) res._base.resize(l2, 0);
        for (int i = 0; i < l2; i++) {
            res._base[i] += b2._base[i] + carry;
            carry = (res._base[i] >= bignum::unit_base);
            if (carry) res._base[i] -= bignum::unit_base;
        }
        for (int i = l2; carry && i < l1; i++) {
            res._base[i] += carry;
            carry = (res._base[i] >= bignum::unit_base);
            if (carry) res._base[i] -= bignum::unit_base;
        }
        if (carry) res._base.push_back(carry);
        res.trim();
        return res;
    }
    static bignum abssub(const bignum& b1, const bignum& b2) { // assert(|b1|>=|b2|)
        bignum res(b1._base, 0);
        int brw = 0;
        for (int i = 0, iEnd = b2._base.size(); i < iEnd; i++) {
            res._base[i] -= b2._base[i] + brw;
            brw = (res._base[i] < 0);
            if (brw) res._base[i] += bignum::unit_base;
        } assert(brw == 0);
        res.trim();
        return res;
    }
    static int abscomp(const bignum& b1, const bignum& b2) {
        int l1 = b1._base.size(), l2 = b2._base.size();
        if (l1 != l2) return l1 - l2;
        for (int i = l1 - 1; i >= 0; i--)
            if (b1._base[i] != b2._base[i])
                return b1._base[i] - b2._base[i];
        return 0;
    }
    static void fft_shuffle(bignum::complex a[], int len) {
        std::vector<int> order(len, 0);
        int l=-2;for(int i=len;i;i>>=1)l++;//l=log(len)-1
        for (int i = 0; i < len; i++)
            order[i] = (order[i >> 1] >> 1) | ((i & 1) << l);
        for (int i = 0; i < len; i++)
            if (i < order[i])
                std::swap(a[i], a[order[i]]);
    }
    static void fft(bignum::complex a[], int len, int mode) { // assert(len is 2's power)
        fft_shuffle(a, len);
        for (int n = 2, ndiv2 = 1; n <= len; n <<= 1, ndiv2 <<= 1) {
            bignum::complex wn(cos(PI / ndiv2), mode * sin(PI / ndiv2));
            for (int bst = 0; bst < len; bst += n) {
                bignum::complex w(1, 0);
                for (int k = bst, ke = bst + ndiv2; k < ke; k++, w *= wn) {
                    bignum::complex r = w * a[k + ndiv2]; // butfyswap(k, k+n/2)
                    a[k + ndiv2] = a[k] - r;
                    a[k] += r;
                }
            }
        }
        if (mode == -1)
            for (int i = 0; i < len; i++)
                a[i].re /= len, a[i].im /= len;
    }
    
public:
    friend bool operator==(const bignum& b1, long long x) { return b1 == bignum(x); }
    friend bool operator!=(const bignum& b1, long long x) { return !(b1 == bignum(x)); }
    friend bool operator<(const bignum& b1, long long x) { return b1 < bignum(x); }
    friend bool operator>=(const bignum& b1, long long x) { return !(b1 < bignum(x)); }
    friend bool operator>(const bignum& b1, long long x) { return bignum(x) < b1; }
    friend bool operator<=(const bignum& b1, long long x) { return !(bignum(x) < b1); }
    friend bool operator==(const bignum& b1, const bignum& b2) { return b1._sign == b2._sign && b1._base == b2._base; }
    friend bool operator!=(const bignum& b1, const bignum& b2) { return !(b1 == b2); }
    friend bool operator<(const bignum& b1, const bignum& b2) { return b1._sign && !b2._sign || b1._sign ^ b2._sign ^ (bignum::abscomp(b1, b2) < 0); }
    friend bool operator>=(const bignum& b1, const bignum& b2) { return !(b1 < b2); }
    friend bool operator>(const bignum& b1, const bignum& b2) { return b2 < b1; }
    friend bool operator<=(const bignum& b1, const bignum& b2) { return !(b2 < b1); }
    friend bignum operator+(const bignum& b1, long long x) { return b1 + bignum(x); }
    friend bignum operator+(const bignum& b1, const bignum& b2)
    {
        if (b1._sign == b2._sign) return (b1._sign == 0 ? bignum::absadd(b1, b2) : -bignum::absadd(b1, b2));
        else if (b1._sign == 0) return b1 - (-b2); // a >= 0, b < 0
        else return b2 - (-b1); // a < 0, b >= 0
    }
    friend bignum operator-(const bignum& b1, long long x) { return b1-bignum(x); }
    friend bignum operator-(const bignum& b1, const bignum& b2)
    {
        /* 
        0 <= b <  a     (|a|-|b|)
        a <= b <  0     (|a|-|b|)
        b <  a <  0    -(|b|-|a|)
        0 <= a <= b    -(|b|-|a|)
        */
        if (b1._sign == b2._sign) return ((b1._sign ^ (b1 >= b2)) ? bignum::abssub(b1, b2) : -bignum::abssub(b2, b1));
        else if (b1._sign == 0) return bignum::absadd(b1, -b2); // b < 0 <= a
        else return -bignum::absadd(-b1, b2); // a < 0 <= b
    }
    friend bignum operator*(const bignum &b1, long long x){
        bignum ret = b1;
        if (x < 0) {
            x = -x;
            ret._sign ^= 1;
        }
        long long carry = 0;
        for (int i = 0, iEnd=ret._base.size(); i < iEnd; i++){
            ret._base[i] = ret._base[i] * x + carry;
            carry = (ret._base[i] / bignum::unit_base);
            if (carry) ret._base[i] %= bignum::unit_base;
        }
        while (carry) {
            ret._base.push_back(carry % bignum::unit_base);
            carry /= bignum::unit_base;
        }
        ret.trim();
        return ret;
    }
    friend bignum operator*(const bignum& _b1, const bignum& b2)
    {
        bignum b1(_b1._base, _b1._sign ^ b2._sign);
        int len = 1, l1 = _b1._base.size(), l2 = b2._base.size();
        for (int s1 = 2 * l1, s2 = 2 * l2; len < s1 || len < s2;) len <<= 1;
        b1._base.resize(len);
        std::vector<bignum::complex> x1(len), x2(len);
        for (int i = 0; i < len; i++) {
            x1[i] = bignum::complex(b1._base[i], 0.0);
            if(i < l2) x2[i] = bignum::complex(b2._base[i], 0.0);
        }
        fft(x1.data(), len, 1); fft(x2.data(), len, 1);
        for (int i = 0; i < len; i++) x1[i] *= x2[i];
        fft(x1.data(), len, -1);
        for (int i = 0; i < len; i++) b1._base[i] = round(x1[i].re);
        long long carry = 0;
        for (int i = 0; i < len; i++) {
            b1._base[i] += carry;
            carry = b1._base[i] / bignum::unit_base;
            if (carry) b1._base[i] %= bignum::unit_base;
        } assert(carry==0);
        b1.trim();
        return b1;
    }
    friend bignum operator/(const bignum& b1, long long x)
    {
        bignum ret;
        if (x < 0) {
            x = -x;
            ret._sign ^= 1;
        }
        ret._base.resize( b1._base.size() );
        long long rest = 0;
        for(int i = b1._base.size() - 1 ;i >= 0; i--){
            ret._base[i] = rest * bignum::unit_base + b1._base[i];
            rest = ret._base[i] % x;
            ret._base[i] /= x;
        }
        ret.trim();
        return ret;
    }
    friend bignum operator/(const bignum& b1, const bignum& b2)
    {
        if (bignum::abscomp(b1, b2) < 0) return bignum();
        int cmp;
        bignum cur(b2._base, 0);
        std::vector<bignum> p2ntimes;
        while ((cmp = bignum::abscomp(cur, b1)) <= 0)
        {
            p2ntimes.push_back(cur); cur += cur;
            if (p2ntimes.size() >= powTwo.size())
                powTwo.push_back(powTwo.back() + powTwo.back());
        }
        bignum quo = powTwo[p2ntimes.size() - 1];
        cur = p2ntimes.back();
        for (int i = p2ntimes.size() - 2; i >= 0 && cmp != 0; --i)
            if ((cmp = bignum::abscomp(cur + p2ntimes[i], b1)) <= 0)
                cur += p2ntimes[i], quo += powTwo[i];
        quo._sign = b1._sign ^ b2._sign;
        return quo;
    }
    friend bignum operator%(const bignum& b1, long long x){
        long long d = 0;
        for (int i = b1._base.size() - 1; i >= 0; i--)
            d = (d * (long long)bignum::unit_base + (long long)b1._base[i]) % x;
        return d;
    }
    friend bignum operator%(const bignum& b1, const bignum& b2){
        if (bignum::abscomp(b1, b2) < 0) return b1;
        int cmp;
        bignum cur = b2;
        std::vector<bignum> p2ntimes;
        while ((cmp = bignum::abscomp(cur, b1)) <= 0)
        {
            p2ntimes.push_back(cur), cur += cur;
            if (p2ntimes.size() >= powTwo.size())
                powTwo.push_back(powTwo.back() + powTwo.back());
        }
        cur = p2ntimes.back();
        for (int i = p2ntimes.size() - 1; i >= 0 && cmp != 0; --i)
            if (bignum::abscomp(cur + p2ntimes[i], b1) <= 0)
                cur += p2ntimes[i];
        return b1 - cur;
    }
    friend void swap(bignum& b1, bignum& b2){
        std::swap(b1._base,b2._base);
        std::swap(b1._sign,b2._sign);
    }
};

#ifndef BIGNUM_INITIALIZED
#define BIGNUM_INITIALIZED
std::vector<bignum> bignum::powTwo = { bignum(1) };
#endif