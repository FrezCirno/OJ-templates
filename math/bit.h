

// gcc buildin function (preferred)
// int __builtin_popcount(unsigned int)
// int __builtin_popcountll(unsigned long long)
// int __builtin_parity(unsigned int)
// int __builtin_parityll(unsigned long long)

int popcount(unsigned int x)
{
    int i;
    for (i = 0; x; x &= x - 1, i++) ;
    return i;
}

int popcount(unsigned long long x)
{
    int i;
    for (i = 0; x; x &= x - 1, i++) ;
    return i;
}

template <class T>
int popcount(T x)
{
    int i;
    for (i = 0; x; x &= x - 1, i++) ;
    return i;
}

bool isPowerOfTwo(int n) { return n > 0 && (n & (n - 1)) == 0; }
