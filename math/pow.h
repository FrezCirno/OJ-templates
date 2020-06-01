/*
 * Fast Pow Algorithm
 *
 */

int pow(int a, int b)
{
    int res = 1;
    while (b) {
        if (b & 1) res *= a;
        a *= a;
        b >>= 1;
    }
    return res;
}

long long pow(long long a, int b)
{
    long long res = 1;
    while (b) { 
        if (b & 1) res *= a;
        a *= a;
        b >>= 1;
    }
    return res;
}

template <class T>
T pow(T a, int b)
{
    T res = 1;
    while (b) {
        if (b & 1) res *= a;
        a *= a;
        b >>= 1;
    }
    return res;
}

int modpow(int a, int b, int m)
{
    a %= m;
    int res = 1;
    while (b) {
        if (b & 1) res = res * a % m;
        a = a * a % m;
        b >>= 1;
    }
    return res;
}

long long modpow(long long a, int b, long long m)
{
    a %= m;
    long long res = 1;
    while (b) {
        if (b & 1) res = res * a % m;
        a = a * a % m;
        b >>= 1;
    }
    return res;
}

template <class T>
T modpow(T a, int b, T m)
{
    a %= m;
    T res = 1;
    while (b) {
        if (b & 1) res = res * a % m;
        a = a * a % m;
        b >>= 1;
    }
    return res;
}