#include <utility>

int gcd(int a, int b) {
    while (b) { a %= b; std::swap(a, b); }
    return a;
}

long long gcd(long long a, long long b) {
    while (b) { a %= b; std::swap(a, b); }
    return a;
}

template <class T>
T gcd(T a, T b) {
    while (b) { a %= b; std::swap(a, b); }
    return a;
}