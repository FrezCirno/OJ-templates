/* 
 * Online Judge Debug Macro -- very easy to use
 * 
 * Usage: 
 * 1. #define _DEBUG on your local envionment
 * 2. Use debug(anything);
 */


#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <utility>

#ifdef _DEBUG
 template <class T> void _E(T x) { std::cerr << x; } 
 void _E(const char* s) { std::cerr << '"' << s << '"'; } 
 void _E(std::string s) { _E(s.c_str()); } 
 template <class A, class B> void _E(std::pair<A, B> x) { std::cerr << '('; _E(x.first); std::cerr << ", "; _E(x.second); std::cerr << ')'; } 
 template <class T> void _E(std::vector<T> x) { std::cerr << '['; for (auto it = x.begin(); it != x.end(); ++it) { if (it != x.begin()) std::cerr << ", "; _E(*it); } std::cerr << ']'; }
 void ERR() {} 
 template <class A, class... B> void ERR(A x, B... y) { _E(x); std::cerr << (sizeof...(y) ? ", " : " "); ERR(y...); }
 #define debug(...) do { std::cerr << "Debug {" #__VA_ARGS__ "} => "; ERR(__VA_ARGS__); std::cerr<<endl; } while(false)
#else
 #define debug(...) 0
#endif
