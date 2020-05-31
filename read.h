/*
 * Read an ingeter by getchar()
 * 
 * A bit faster than scanf, cin 
 */

#include <cstdio>

int read(){
	int s=0,f=0;
	char c=getchar();
	while(c>'9'||c<'0'){if(c=='-')f=1;c=getchar();}
	while(c>='0'&&c<='9'){s=(s<<3)+(s<<1)+c-'0';c=getchar();}
	return(f?-s:s);
}