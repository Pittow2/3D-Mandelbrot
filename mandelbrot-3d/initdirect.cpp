#include<cstdio>
#include<cmath>
#include<cstdlib>
double a[4],k,l;
int main(){
	FILE *in=fopen("direct.txt","r");
	FILE *out=fopen("config.txt","w");
	for(int i=0;i<4;++i){
		fscanf(in,"%lf",a+i);
		l+=a[i]*a[i];
	}
	fscanf(in,"%lf",&k);
	l=sqrt(l);
	for(int i=0;i<4;++i)a[i]/=l;
	if(a[3]>0){
		for(int i=0;i<4;++i)a[i]=-a[i];
		k=-k;
	}
	a[3]-=1;
	for(int i=0;i<3;++i){
		for(int j=0;j<4;++j){
			double p=a[i]*a[j]/a[3];
			if(i==j)p+=1;
			fprintf(out,"%.9lf ",p);
		}
		fprintf(out,"\n");
	}
	a[3]+=1;
	for(int i=0;i<4;++i)fprintf(out,"%.9lf ",a[i]*k);
	fclose(in);
	fclose(out);
	system("3d");
	return 0;
}
