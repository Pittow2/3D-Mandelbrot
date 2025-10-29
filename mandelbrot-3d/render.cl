#define __OPENCL_C_VERSION__ 300

#define G __global
#define GE __generic
#define K __kernel void
#define A __attribute__((aligned(4)))
typedef unsigned char byte;
typedef double real;

typedef struct _complex{
	real x,y;
}A complex;
typedef struct _iter{
	complex c,z,z0;
	int id;//oct node id
	int sta;
	/**	+x:playing x-1 steps
		 0:caged
		-x:escaped after -1-x steps
	 */
}A iter;
typedef struct _pos{
	union{
		real a[3];
		struct{
			real x,y,z;
		};
	};
}A pos;
typedef struct _octree{
	int s[2];//subnode id,split x,y,z in order
	int r[3];//neighbor hoods in x,y,z direct, further
	int da,de;//dad id,depth
	union{
		int w;
		struct{
			byte f[2],gap,lz;
		};
	};
	/**	-1:playing
		-2:splitted
		-3:tag to split
		 x:escaped after x steps
		-
		caged lazy:lz=128
		f:pixel touch flag
	 */
}A octree;
typedef struct _pixel{
	pos p,v,cu;//trace pos,velo,oct node center
	real add;
	int id;//oct node id(-1 for end)
	int ans;
	/**	-1:caged
		 x:max escape after x steps
	 */
}A pixel;
typedef struct _camera{
	pos p;//center
	real a[3][3];//rotator
}A camera;
typedef struct _proj{
	real a[4][4];
	//x*a[0]+y*a[1]+z*a[2]+a[3]
}A proj;

inline real abs(real x){
	return x<0?-x:x;
}
inline void cadd(GE complex *a,GE complex *b){//a+=b;
	a->x+=b->x;
	a->y+=b->y;
}
inline bool csquare(GE complex *a){//a=a*a;
	real x=a->x*a->x;
	real y=a->y*a->y;
	if(x+y>=4)return 1;
	x-=y;
	y=a->x*a->y;
	a->x=x;
	a->y=y+y;
	return 0;
}
inline bool cequal(GE complex *a,GE complex *b){//a==b;
	return a->x==b->x&&a->y==b->y;
}

K initial(G octree *s){
	s[0].w=-3;
	s[0].r[0]=s[0].r[1]=s[0].r[2]=-1;
	s[0].da=-1;
	s[0].de=0;
}
K initframe(int sx,int sy,int sz,G camera *src,G pixel *tar){
	const int x=get_global_id(0),y=get_global_id(1);
	const int p=x+y*sx;
	tar[p].p=src->p;
	tar[p].v=(pos){0,0,0};
	pos w=(pos){x-(sx>>1),y-(sy>>1),sz};
	for(int i=0;i<3;++i)
	for(int j=0;j<3;++j)tar[p].v.a[j]+=w.a[i]*src->a[i][j];
	tar[p].cu=(pos){0,0,0};
	tar[p].add=1;
	tar[p].id=0;
	tar[p].ans=0;
}
#define jumpup {\
	int pid=map[p->id].da,we=map[pid].de%3;\
	if(we==2)p->add*=2;\
	if(map[pid].s[0]==p->id)p->cu.a[we]+=p->add;\
	else p->cu.a[we]-=p->add;\
	p->id=pid;\
	if(!p->id)p->id=-1;\
	if(p->id==-1)return;\
}
K trace(G pixel *w,G octree *map){
	pixel *p=&w[get_global_id(0)];
	while(p->id!=-1){
		octree *ma;
		while(ma=&map[p->id],ma->w==-2){
			int wi=ma->de%3;
			if(p->p.a[wi]<p->cu.a[wi]){
				p->id=ma->s[0];
				p->cu.a[wi]-=p->add;
			}
			else{
				p->id=ma->s[1];
				p->cu.a[wi]+=p->add;
			}
			if(wi==2)p->add/=2;
		}
		if(ma->w==-1)return;
		if(ma->lz==128){
			ma->f[get_global_id(0)&1]=1;
			p->ans=-1;
			return;
		}
		if(p->ans<ma->w)p->ans=ma->w;
		pos g;
		int wi=ma->de%3;
		for(int i=0;i<3;++i){
			g.a[i]=p->cu.a[i]-p->p.a[i];
			if(p->v.a[i]<0)g.a[i]=-g.a[i];
			g.a[i]+=p->add;
			if(i>=wi)g.a[i]+=p->add;
		}
		for(int i=0;i<3;++i)
			if(abs(g.a[i]*p->v.a[wi])<abs(g.a[wi]*p->v.a[i]))wi=i;
		bool pr=p->v.a[wi]>0;
		real dt=abs(g.a[wi]/p->v.a[wi]);
		for(int i=0;i<3;++i)p->p.a[i]+=dt*p->v.a[i];
		//optimize with neighbor - to do
		while((map[p->id].de+2)%3!=wi)jumpup;
		while((map[map[p->id].da].s[1]==p->id)==pr)
			for(int i=0;i<3;++i)jumpup;
		p->id=map[map[p->id].da].s[pr];
		real a=p->add*2;
		if(map[p->id].de%3==0)a+=a;
		if(!pr)a=-a;
		p->cu.a[wi]+=a;
	}
}
K sumstage0(G int *a,int r){
	int t=(((get_global_id(0)<<1)+1)<<r)-1;
	a[t+(1<<r)]+=a[t];
}
K sumstage1(G int *a,int r){
	int t=(((get_global_id(0)<<1)+1)<<r)-1;
	a[t]+=a[t-(1<<r)];
}
K initnodelist(G octree *s,G int *tar){
	const int id=get_global_id(0);
	if(s[id].w==-3)tar[id]=1;
	else tar[id]=0;
}
K allocnode(int tc,int wc,G octree *s,G int *tar,G iter *w){
	/**	tc:tree node count
		wc:new iter list count
		tar:leaf list presum
		w:new iter list
	 */
	const int id=get_global_id(0);
	if(s[id].w!=-3)return;
	//init node
	int add=tar[id]-1,ls=tc+(add<<1),ws=wc+(add<<1);
	s[id].s[0]=ls;
	s[id].s[1]=ls+1;
	s[ls].da=s[ls+1].da=id;
	s[ls].de=s[ls+1].de=s[id].de+1;
	s[ls].w=s[ls+1].w=-1;
	//set neighbor
	int wi=s[id].de%3;
	if(s[id].de>=3){
		int p0,p1,p2,pe;
		p0=s[id].da;
		p1=s[p0].da;
		p2=s[p1].da;
		bool w=s[p2].s[0]==p1;
		pe=s[p2].s[w];
		if(s[pe].w==-2)pe=s[pe].s[s[p1].s[1]==p0];
		if(s[pe].w==-2)pe=s[pe].s[s[p0].s[1]==id];
		if(s[pe].w==-2)pe=s[pe].s[!w];
		s[ls+w].r[wi]=pe;
		pe=s[id].r[wi];
		if(pe!=-1&&s[pe].de==s[id].de&&s[pe].w==-2)pe=s[pe].s[w];
		s[ls+!w].r[wi]=pe;
	}
	else s[ls].r[wi]=s[ls+1].r[wi]=-1;
	for(int i=0;i<3;++i)if(i!=wi){
		int pe=s[id].r[i];
		if(pe!=-1&&s[pe].de==s[id].de&&s[pe].w==-2){
			s[ls].r[i]=s[pe].s[0];
			s[ls+1].r[i]=s[pe].s[1];
		}
		else s[ls].r[i]=s[ls+1].r[i]=pe;
	}
}
K initpos(G octree *s,G iter *w,int dis,G proj *p){
	//init iter data
	int id=get_global_id(0);
	const int ws=id+dis;
	w[ws].id=id;
	w[ws].sta=1;
	real a[4]={0,0,0,1},b[4]={0,0,0,0};
	while(id){
		int pid=s[id].da;
		real *c=&a[s[pid].de%3];
		*c/=2;
		if(s[pid].s[0]==id)*c-=1;
		else *c+=1;
		id=pid;
	}
	for(int i=0;i<4;++i)
	for(int j=0;j<4;++j)b[j]+=a[i]*p->a[i][j];
	w[ws].c.x=b[0];
	w[ws].c.y=b[1];
	w[ws].z.x=b[2];
	w[ws].z.y=b[3];
	w[ws].z0=w[ws].z;
}
K setnodestatus(G octree *s){
	const int id=get_global_id(0);
	if(s[id].w==-3)s[id].w=-2;
}
K inititerlist(G octree *s,G iter *w,G int *tar){
	const int id=get_global_id(0);
	tar[id]=0;
	int *t=&s[w[id].id].w;
	if(*t==-2||*t==-3){
		w[id].sta=-1;
		return;
	}
	if(w[id].sta==0)*t=INT_MIN;//lz=128,f0=f1=0
	else if(w[id].sta<0)*t=-1-w[id].sta;
	else tar[id]=1;
}
K allociter(G iter *a,G int *tar,G iter *b){
	const int id=get_global_id(0);
	if(a[id].sta>0)b[tar[id]-1]=a[id];
}
K resetstatus(G octree *s){
	octree *p=&s[get_global_id(0)];
	if(p->lz==128)p->w=INT_MIN;
}
K splitcheck(G octree *s){
	int id=get_global_id(0);
	octree *pi=&s[id];
	if(pi->de<6&&pi->w!=-2){
		pi->w=-3;
		return;
	}
	if(pi->lz!=128)return;
	if(pi->f[0]&&pi->f[1])pi->w=-3;
	else{
		pi->w=INT_MIN;//lz=128,f0=f1=0
		return;
	}
	bool w[3];
	int p[4];
	p[0]=id;
	for(int i=1;i<4&&i<=pi->de;++i){
		p[i]=s[p[i-1]].da;
		w[i-1]=s[p[i]].s[1]==p[i-1];
		int pe=s[p[i]].s[!w[i-1]];
		for(int j=i-2;~j&&s[pe].w==-2;--j)pe=s[pe].s[w[j]];
		if(s[pe].de<pi->de&&s[pe].w!=-2)s[pe].w=-3;
	}
	for(int i=0;i<3;++i){
		int pe=pi->r[i];
		if(pe==-1)continue;
		while(s[pe].de<pi->de&&s[pe].w==-2){
			int p=id;
			for(int i=pi->de-s[pe].de-1;i;--i)p=s[p].da;
			int pa=s[p].da;
			bool w=s[pa].s[1]==p;
			w^=s[pa].de%3==i;
			pe=s[pe].s[w];
		}
		pi->r[i]=pe;
		if(s[pe].de<pi->de&&s[pe].w!=-2)s[pe].w=-3;
	}
}
K emu(int step,G iter *d){
	const int i=get_global_id(0);
	while(d[i].sta>0&&step--){
		++d[i].sta;
		if(csquare(&d[i].z))d[i].sta=-d[i].sta;
		else{
			cadd(&d[i].z,&d[i].c);
			if(d[i].sta&1){
				csquare(&d[i].z0);
				cadd(&d[i].z0,&d[i].c);
				if(cequal(&d[i].z,&d[i].z0))d[i].sta=0;
			}
		}
		// if(d[i].sta>256)d[i].sta=0;
	}
}
K render(G camera *c,G pixel *w,G byte *t){
	const int id=get_global_id(0);
	if(w[id].ans==-1){
		real d=0,rate=75;
		for(int i=0;i<3;++i){
			real ca=w[id].p.a[i]-c->p.a[i];
			d+=ca*ca;
		}
		d=log(d)*7;
		t[id*3]=rate*(sin(d+6.28*2/3)+1);
		t[id*3+1]=rate*(sin(d+6.28/3)+1);
		t[id*3+2]=rate*(sin(d)+1);
	}
	else{
		real d=w[id].ans/16.0,rate=128;
		t[id*3]=rate*(sin(d)+1);
		t[id*3+1]=rate*(sin(d+6.28/3)+1);
		t[id*3+2]=rate*(sin(d+6.28*2/3)+1);
	}
}
K copy(G octree *a,G octree *b){
	b[get_global_id(0)]=a[get_global_id(0)];
}