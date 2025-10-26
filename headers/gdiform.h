#ifndef _GDIFORM_H
#define _GDIFORM_H

#undef UNICODE

#include<windows.h>
#define Press GetAsyncKeyState

typedef unsigned char byte;
typedef unsigned int i32;
typedef unsigned long long u64;

struct rgb{
	byte r,g,b;
};
struct x2d{
	double x,y;
};
void swap(double &a,double &b){
	double c=a;a=b;b=c;
}

struct window{
	bool end,focus,hide;
	bool fc;
	int px,py,ln,sx,sy;
	byte *data;
	char *Title;
	WNDCLASSEXA wc;
	HWND hw;
	HDC hdc,mdc;
	HBITMAP map;
	BITMAPINFO bi;
	int id;
	void InitialGDI(int,int,int,int,const char*);
	void ready(bool);
	void sethide(bool);
	void paste();
	void drawpoint(int,int,rgb);
	void drawline(x2d,x2d,rgb);
}*wi[256];
int wcnt;

rgb getrgb(double h){
	while(h>=1)h-=1;
	double x=h*6,r=0,g=0,b=0;
	while(x>=2)x-=2;
	if(x>1)x=2-x;
	int a=h*6;
	if(a==0||a==5)r=1;
	else if(a==1||a==2)g=1;
	else b=1;
	if(a==1||a==4)r=x;
	else if(a==0||a==3)g=x;
	else b=x;
	byte r0=r*255,g0=g*255,b0=b*255;
	return (rgb){r0,g0,b0};
}

void keys(int,u64);
void click(int,int,int,bool);
long long callback(HWND hw,i32 msg,u64 wp,long long lp){
	int id;
	for(id=0;id<wcnt&&wi[id]->hw!=hw;++id);
	if(id==wcnt)return DefWindowProcA(hw,msg,wp,lp);
	window *g=wi[id];
	if(msg==WM_PAINT)BitBlt(g->hdc,0,0,g->sx,g->sy,g->mdc,0,0,SRCCOPY);
	else if(msg==WM_SETFOCUS){
		g->focus=1;
		if(g->hide)ShowCursor(SW_HIDE);
	}
	else if(msg==WM_KILLFOCUS){
		g->focus=0;
		if(g->hide)ShowCursor(SW_SHOW);
	}
	else if(msg==WM_MOVING){
		RECT *prc=(RECT*)lp;
		g->px=prc->left;
		g->py=prc->top;
	}
	else if(msg==WM_SIZING){
		RECT *prc=(RECT*)lp;
		prc->left=g->px;
		prc->top=g->py;
		prc->right=g->px+g->sx;
		prc->bottom=g->py+g->sy;
		return 1;
	}
	else if(msg==WM_LBUTTONDOWN||msg==WM_RBUTTONDOWN){
		int y=lp>>16,x=lp-(y<<16);
		bool r=0;
		if(msg==WM_RBUTTONDOWN)r=1;
		if(g->fc)click(id,x,y,r);
	}
	else if(msg==WM_CLOSE||(msg==WM_CHAR&&wp==VK_ESCAPE)){
		if(MessageBoxA(hw,"Exit?",g->Title,MB_ICONQUESTION|MB_YESNO|MB_APPLMODAL)==6)
			DestroyWindow(hw);
		return 0;
	}
	else if(msg==WM_DESTROY)g->end=0;
	else if(msg==WM_CHAR)keys(id,wp);
	return DefWindowProcA(hw,msg,wp,lp);
}
void window::paste(){
	SetDIBits(mdc,map,0,sy,data,&bi,DIB_RGB_COLORS);
	BitBlt(hdc,0,0,sx,sy,mdc,0,0,SRCCOPY);
}
void window::ready(bool w){
	if(w){
		ShowWindow(hw,SW_SHOW);
		SetForegroundWindow(hw);
	}else ShowWindow(hw,SW_HIDE);
}
void window::sethide(bool w){hide=w;}
void window::InitialGDI(int _px,int _py,int _sx,int _sy,const char *title){
	Title=(char*)title;
	px=_px;py=_py;sx=_sx;sy=_sy;
	hide=0;
	wc.cbSize=sizeof(WNDCLASSEX);
	wc.lpfnWndProc=callback;
	wc.hCursor=LoadCursorA(0,IDC_ARROW);
	wc.hbrBackground=(HBRUSH)(COLOR_WINDOW+1);
	wc.lpszClassName="WindowClass";
	wc.hIcon=LoadIconA(0,IDI_APPLICATION);
	wc.hIconSm=LoadIconA(0,IDI_APPLICATION);
	RegisterClassExA(&wc);
	hw=CreateWindowExA(0,wc.lpszClassName,Title,
		WS_SYSMENU|WS_MINIMIZEBOX|WS_THICKFRAME,
		px,py,sx,sy,0,0,0,0);
	hdc=GetDC(hw);
	mdc=CreateCompatibleDC(hdc);
	map=CreateCompatibleBitmap(hdc,sx,sy);
	SelectObject(mdc,map);
	bi.bmiHeader.biSize=sizeof(BITMAPINFOHEADER);
	bi.bmiHeader.biWidth=sx;
	bi.bmiHeader.biHeight=-sy;
	bi.bmiHeader.biPlanes=1;
	bi.bmiHeader.biBitCount=24;
	bi.bmiHeader.biCompression=BI_RGB;
	bi.bmiHeader.biSizeImage=0;
	bi.bmiHeader.biClrUsed=0;
	bi.bmiHeader.biClrImportant=0;
	ln=sx*3;
	while(ln&3)++ln;
	data=new byte[ln*sy];
	end=1;
	id=wcnt;
	wi[wcnt++]=this;
}
void message(){
	MSG msg;
	while(PeekMessageA(&msg,0,0,0,PM_REMOVE)){
		TranslateMessage(&msg);
		DispatchMessageA(&msg);
	}
	for(int i=0;i<wcnt;++i)wi[i]->fc=wi[i]->focus;
}

inline void window::drawpoint(int x,int y,rgb c){
	if(x<0||x>=sx||y<0||y>=sy)return;
	data[x*3+y*ln+2]=c.r;
	data[x*3+y*ln+1]=c.g;
	data[x*3+y*ln+0]=c.b;
}
inline void window::drawline(x2d f,x2d t,rgb c){
	double dx=t.x-f.x,dy=t.y-f.y,a,b;
	bool flag=0;
	int lx=sx,ly=sy;
	if(abs(dx)<abs(dy)){
		flag=1;
		swap(f.x,f.y);
		swap(t.x,t.y);
		swap(dx,dy);
		lx=sy;ly=sx;
	}
	a=f.y-f.x*dy/dx;
	if(f.x>t.x)swap(f.x,t.x),swap(f.y,t.y);
	if(f.x<0)f.x=0;
	if(t.x>=lx)t.x=lx-1;
	for(int i=f.x;i<=t.x;++i){
		b=i*dy/dx+a;
		if(b>=0&&b<ly){
			if(flag)drawpoint(b,i,c);
			else drawpoint(i,b,c);
		}
	}
}

#endif
