#ifndef _FPS_H
#define _FPS_H

#include<cstdio>
#include<ctime>
struct timer{
	int t,f;
	void init(){t=clock();f=0;}
	void flush(){
		++f;
		int t0=clock();
		if(t0-t>=1000){
			printf("\r%.2lffps    ",1000.0*f/(t0-t));
			t=t0;f=0;
		}
	}
};

#endif
