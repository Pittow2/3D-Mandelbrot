#include"../headers/clform.h"
#include"../headers/gdiform.h"
#include<cstdio>
#include<ctime>
#include<cmath>

#define A __attribute__((aligned(4)))
typedef unsigned char byte;
typedef double real;

struct A complex{
	real x,y;
};
struct A iter{
	complex c,z,z0;
	int id,sta;
};
struct A pos{
	real a[3];
};
struct A octree{
	int s[2],r[3];
	int da,de,w;
};
struct A pixel{
	pos p,v,cu;
	real add;
	int id,ans;
};
struct A camera{
	pos p;
	real a[3][3];
	void move(pos v){
		for(int i=0;i<3;++i)
		for(int j=0;j<3;++j)p.a[i]+=v.a[j]*a[j][i];
	}
	void rotate(real rad,bool w){
		int i,j,k;
		real s=sin(rad),c=cos(rad);
		real ra[3][3],ca[3][3];
		for(i=0;i<3;++i)
		for(j=0;j<3;++j)ra[i][j]=0;
		int u=!w,v=w;
		ra[u][u]=1;
		ra[v][v]=ra[2][2]=c;
		ra[v][2]=-s;
		ra[2][v]=s;
		for(i=0;i<3;++i)
		for(j=0;j<3;++j)ca[i][j]=a[i][j],a[i][j]=0;
		for(i=0;i<3;++i)
		for(j=0;j<3;++j)
		for(k=0;k<3;++k)a[i][j]+=ra[i][k]*ca[k][j];
		for(i=0;i<3;++i){
			s=c=0;
			for(j=0;j<3;++j)s+=a[i][j]*a[i][j],c+=a[j][i]*a[j][i];
			s=sqrt(s);
			c=sqrt(c);
			for(j=0;j<3;++j)a[i][j]/=s,a[j][i]/=c;
		}
	}
};
struct A proj{
	real a[4][4];
	//x*a[0]+y*a[1]+z*a[2]+a[3]
};

int px,py,sx,sy,sz;
window w0;
source s0;
kernel init,iframe,trace,sum0,sum1,inode,anode,ipos,snode,
		iiter,aiter,resta,split,emu,render,copy;
memory node,nsum,ilist[2],isum,pix,cam,pj,bmp;
int nc,nb,ns,ic,ib[2],sc;//node count,node sum count,iter count,iter sum count
camera ca;
proj pr;

void initcl(){
	GetCLAPI();
	InitialCL();
	s0.Load("render.cl");
	init.Load(s0,"initial");
	iframe.Load(s0,"initframe");
	trace.Load(s0,"trace");
	sum0.Load(s0,"sumstage0");
	sum1.Load(s0,"sumstage1");
	inode.Load(s0,"initnodelist");
	anode.Load(s0,"allocnode");
	ipos.Load(s0,"initpos");
	snode.Load(s0,"setnodestatus");
	iiter.Load(s0,"inititerlist");
	aiter.Load(s0,"allociter");
	resta.Load(s0,"resetstatus");
	split.Load(s0,"splitcheck");
	emu.Load(s0,"emu");
	render.Load(s0,"render");
	copy.Load(s0,"copy");
}
void setnode(){
	trace.Setptr(1,node);
	inode.Setptr(0,node);
	anode.Setptr(2,node);
	ipos.Setptr(0,node);
	snode.Setptr(0,node);
	iiter.Setptr(0,node);
	resta.Setptr(0,node);
	split.Setptr(0,node);
}
void setnlist(){
	inode.Setptr(1,nsum);
	anode.Setptr(3,nsum);
}
void setiter(){
	anode.Setptr(4,ilist[1]);
	ipos.Setptr(1,ilist[0]);
	iiter.Setptr(1,ilist[0]);
	iiter.Setptr(2,isum);
	aiter.Setptr(0,ilist[0]);
	aiter.Setptr(1,isum);
	aiter.Setptr(2,ilist[1]);
	emu.Setptr(1,ilist[0]);
}
void initargs(){
	iframe.Setint(0,sx);
	iframe.Setint(1,sy);
	iframe.Setint(2,sz);
	iframe.Setptr(3,cam);
	iframe.Setptr(4,pix);
	trace.Setptr(0,pix);
	ipos.Setptr(3,pj);
	render.Setptr(0,cam);
	render.Setptr(1,pix);
	render.Setptr(2,bmp);
	setnode();
	setnlist();
	setiter();
}
void reset(){
	nc=1;
	ic=0;
	init.Setptr(0,node);
	init.call1(0,1);
}
void initdata(){
	for(int i=0;i<3;++i)ca.a[i][i]=1;
	FILE *fi=fopen("config.txt","r");
	for(int i=0;i<4;++i)
	for(int j=0;j<4;++j)fscanf(fi,"%lf",&pr.a[i][j]);
	fclose(fi);
	ca.p=(pos){1,-1,-1};
	nb=ns=1;
	ib[0]=ib[1]=sc=0;
	node.Buffout(sizeof(octree));
	nsum.Buffout(sizeof(int));
	pix.Buffout(sizeof(pixel)*sx*sy);
	cam.Buffin(sizeof(camera),&ca);
	pj.Buffin(sizeof(proj),&pr);
	bmp.Buffout(sx*sy*3);
	reset();
}
template<typename T>
void swap(T &a,T &b){
	T c=a;
	a=b;b=c;
}

void presum(memory s,int n,bool log=0){
	int ti=clock(),i;
	sum0.Setptr(0,s);
	sum1.Setptr(0,s);
	for(i=0;(1<<i)<=n;++i){
		sum0.Setint(1,i);
		sum0.call1(0,n>>(i+1),0,0);
	}
	for(--i;~i;--i){
		sum1.Setint(1,i);
		sum1.call1(1,((n>>i)-1)>>1,0,0);
	}
	kernel::wait();
	if(log)printf("presum time:%dms\n",clock()-ti);
}
void start(bool w=0){
	if(w)puts("Start:");
	cam.Buffwrite(&ca);
	iframe.call2(0,0,sx,sy,w);
}
void play(int step,bool wa=0,bool w=0){
	if(wa)puts("\nEmu:");
	int ti=clock(),cc,cd;
	trace.call1(0,sx*sy,w);
	split.call1(0,nc,w);
	if(nc>ns){
		nsum.Clear();
		while(ns<nc)ns<<=1;
		nsum.Buffout(ns*sizeof(int));
		setnlist();
	}
	inode.call1(0,nc,w);
	presum(nsum,nc,w);
	cc=nsum.Readint(nc-1)*2;//new node count
	if(nc+cc>nb){
		while(nb<nc+cc)nb<<=1;
		memory ni;
		ni.Buffout(nb*sizeof(octree));
		copy.Setptr(0,node);
		copy.Setptr(1,ni);
		copy.call1(0,nc);
		node.Clear();
		node=ni;
		setnode();
	}
	if(ic>sc){
		if(sc)isum.Clear();
		for(sc=1;sc<ic;sc<<=1);
		isum.Buffout(sc*sizeof(int));
		setiter();
	}
	iiter.call1(0,ic,w);
	presum(isum,ic,w);
	cd=ic?isum.Readint(ic-1):0;//running iter count
	if(cc+cd>ib[1]){
		if(ib[1])ilist[1].Clear();
		for(ib[1]=1;ib[1]<cc+cd;ib[1]<<=1);
		ilist[1].Buffout(ib[1]*sizeof(iter));
		setiter();
	}
	aiter.call1(0,ic,w);
	anode.Setint(0,nc);
	anode.Setint(1,cd);
	anode.call1(0,nc,w);
	ic=cc+cd;
	swap(ib[0],ib[1]);
	swap(ilist[0],ilist[1]);
	setiter();
	ipos.Setint(2,cd-nc);
	ipos.call1(nc,cc,w);
	snode.call1(0,nc,w);
	emu.Setint(0,step);
	emu.call1(0,ic,w);
	nc+=cc;
	if(wa)printf("Emu Total:%dms\n\n",clock()-ti);
}
void play0(int step,bool wa=0,bool w=0){
	if(wa)puts("\nRender:");
	int ti=clock(),cc=0,cd;
	trace.call1(0,sx*sy,w);
	iiter.call1(0,ic,w);
	presum(isum,ic,w);
	cd=ic?isum.Readint(ic-1):0;//running iter count
	if(cc+cd>ib[1]){
		if(ib[1])ilist[1].Clear();
		for(ib[1]=1;ib[1]<cc+cd;ib[1]<<=1);
		ilist[1].Buffout(ib[1]*sizeof(iter));
		setiter();
	}
	aiter.call1(0,ic,w);
	ic=cc+cd;
	swap(ib[0],ib[1]);
	swap(ilist[0],ilist[1]);
	setiter();
	emu.Setint(0,step);
	emu.call1(0,ic,w);
	if(wa)printf("Emu Total:%dms\n\n",clock()-ti);
}
void paste(bool w=0){
	render.call1(0,sx*sy,w);
	bmp.Buffread(w0.data);
	w0.paste();
}

bool crt=1;
real speed=1;
void click(int id,int x,int y,bool right){

}
void keys(int id,u64 wp){
	if(wp==VK_SPACE){
		crt=!crt;
		if(crt)resta.call1(0,nc,1);
	}
	if(wp=='c'||wp=='C')speed*=2;
	if(wp=='v'||wp=='V')speed/=2;
}

int main(){
	initcl();
	px=100;py=64;
	sx=1024;sy=768;sz=1000;
	w0.InitialGDI(px,py,sx,sy,"3D Mandelbrot");
	initdata();
	initargs();
	start();
	w0.sethide(1);
	w0.ready(1);
	while(w0.end){
		message();
		if(w0.fc){
			POINT p;
			GetCursorPos(&p);
			int dx=w0.px+(w0.sx>>1),dy=w0.py+(w0.sy>>1);
			int cx=p.x-dx,cy=p.y-dy;
			real add=0.03,mul=0.001;
			pos a={0,0,0};
			bool flag=0;
			if(Press(VK_SHIFT))add*=8;
			if(Press(VK_CONTROL))add/=8;
			add*=speed;
			if(Press('W'))a.a[2]+=add;
			if(Press('S'))a.a[2]-=add;
			if(Press('D'))a.a[0]+=add;
			if(Press('A'))a.a[0]-=add;
			if(Press('Q'))a.a[1]+=add;
			if(Press('E'))a.a[1]-=add;
			if(a.a[0]!=0||a.a[1]!=0||a.a[2]!=0)flag=1;
			if(cx||cy)flag=1;
			ca.move(a);
			ca.rotate(cx*mul,0);
			ca.rotate(cy*mul,1);
			if(flag){
				SetCursorPos(dx,dy);
				start(1);
			}
			if(crt)play(256,1);
			else play0(256,1);
			printf("node count=%d, buf=%d\n",nc,nb);
			printf("iter queue=%d, buf=%d+%d\n\n",ic,ib[0],ib[1]);
			paste(1);
		}
		else Sleep(30);
	}
	return 0;
}
