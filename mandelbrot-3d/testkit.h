#include"../headers/clform.h"
#include<cstdio>
#include<ctime>
bool testpresum(int n){
	int *a,*b;
	memory s;
	a=new int[n];
	b=new int[n];
	s.Buffout(n*sizeof(int));
	for(int i=0;i<n;++i)a[i]=i*(i+5);
	s.Buffwrite(a);
	for(int i=0;i<10;++i)presum(s,1);
	s.Buffread(b);
	int ti=clock();
	for(int k=0;k<10;++k)for(int i=1;i<n;++i)a[i]+=a[i-1];
	printf("cpu time:%dms\n",clock()-ti);
	for(int i=0;i<n;++i)if(a[i]!=b[i])return 0;
	return 1;
}
void testbench(int n,int step,int rep){
	memory s;
	s.Buffout(n*sizeof(int));
	for(int i=step;i<=n;i+=step){
		s.len=i*sizeof(int);
		int ti=clock();
		for(int k=0;k<rep;++k)presum(s,i);
		printf("%d:%dms\n",i,(clock()-ti)/rep);
	}
}
void test(){
	if(testpresum(140000000))puts("Yes");
	else puts("No");
	testbench(220000000,1000000,8);
}