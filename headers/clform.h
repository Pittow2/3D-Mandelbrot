#ifndef _CLFORM_H
#define _CLFORM_H

#include<cstdio>
#include<windows.h>
#define CL_TARGET_OPENCL_VERSION 300
#include"CL/cl.h"

typedef unsigned long u32;
typedef unsigned long long u64;
typedef cl_int(*a0)(cl_uint,cl_platform_id*,cl_uint*);
a0 DclGetPlatformIDs;
typedef cl_int(*a1)(cl_platform_id,cl_platform_info,size_t,void*,size_t*);
a1 DclGetPlatformInfo;
typedef cl_int(*a2)(cl_platform_id,cl_device_type,cl_uint,cl_device_id*,cl_uint*);
a2 DclGetDeviceIDs;
typedef cl_int(*a3)(cl_device_id,cl_device_info,size_t,void*,size_t*);
a3 DclGetDeviceInfo;
typedef cl_context(*a4)(const cl_context_properties*,cl_uint,const cl_device_id*,
	void*(const char*,const void*,size_t,void*),
	void*,cl_int*);
a4 DclCreateContext;
typedef cl_command_queue(*a5)(cl_context,cl_device_id,cl_command_queue_properties,cl_int*);
a5 DclCreateCommandQueue;
typedef cl_program(*a6)(cl_context,cl_uint,const char**,const size_t*,cl_int*);
a6 DclCreateProgramWithSource;
typedef cl_int(*a7)(cl_program,cl_uint,const cl_device_id*,const char*,
	void*(cl_program,void*),
	void*);
a7 DclBuildProgram;
typedef cl_int(*a8)(cl_program,cl_device_id,cl_program_build_info,size_t,void*,size_t*);
a8 DclGetProgramBuildInfo;
typedef cl_kernel(*a9)(cl_program,const char*,cl_int*);
a9 DclCreateKernel;
typedef cl_mem(*b0)(cl_context,cl_mem_flags,size_t,void*,cl_int*);
b0 DclCreateBuffer;
typedef cl_int(*b1)(cl_kernel,cl_uint,size_t,const void*);
b1 DclSetKernelArg;
typedef cl_int(*b2)(cl_command_queue,cl_kernel,cl_uint,const size_t*,const size_t*,const size_t*,cl_uint,const cl_event*,cl_event*);
b2 DclEnqueueNDRangeKernel;
typedef cl_int(*b3)(cl_uint,const cl_event*);
b3 DclWaitForEvents;
typedef cl_int(*b4)(cl_event);
b4 DclReleaseEvent;
typedef cl_int(*b5)(cl_command_queue,cl_mem,cl_bool,size_t,size_t,void*,cl_uint,const cl_event*,cl_event*);
b5 DclEnqueueReadBuffer;
typedef cl_int(*b6)(cl_command_queue,cl_mem,cl_bool,size_t,size_t,const void*,cl_uint,const cl_event*,cl_event*);
b6 DclEnqueueWriteBuffer;
typedef cl_int(*b7)(cl_mem);
b7 DclReleaseMemObject;
typedef cl_int(*b8)(cl_event,cl_profiling_info,size_t,void*,size_t*);
b8 DclGetEventProfilingInfo;
typedef void*(*b9)(cl_command_queue,cl_mem,cl_bool,cl_map_flags,size_t,size_t,cl_uint,const cl_event*,cl_event*,cl_int*);
b9 DclEnqueueMapBuffer;
typedef cl_int(*c0)(cl_command_queue);
c0 DclFinish;

void GetCLAPI(){
	HMODULE mod=LoadLibraryA("OpenCL.dll");
	DclGetPlatformIDs=(a0)GetProcAddress(mod,"clGetPlatformIDs");
	DclGetPlatformInfo=(a1)GetProcAddress(mod,"clGetPlatformInfo");
	DclGetDeviceIDs=(a2)GetProcAddress(mod,"clGetDeviceIDs");
	DclGetDeviceInfo=(a3)GetProcAddress(mod,"clGetDeviceInfo");
	DclCreateContext=(a4)GetProcAddress(mod,"clCreateContext");
	DclCreateCommandQueue=(a5)GetProcAddress(mod,"clCreateCommandQueue");
	DclCreateProgramWithSource=(a6)GetProcAddress(mod,"clCreateProgramWithSource");
	DclBuildProgram=(a7)GetProcAddress(mod,"clBuildProgram");
	DclGetProgramBuildInfo=(a8)GetProcAddress(mod,"clGetProgramBuildInfo");
	DclCreateKernel=(a9)GetProcAddress(mod,"clCreateKernel");
	DclCreateBuffer=(b0)GetProcAddress(mod,"clCreateBuffer");
	DclSetKernelArg=(b1)GetProcAddress(mod,"clSetKernelArg");
	DclEnqueueNDRangeKernel=(b2)GetProcAddress(mod,"clEnqueueNDRangeKernel");
	DclWaitForEvents=(b3)GetProcAddress(mod,"clWaitForEvents");
	DclReleaseEvent=(b4)GetProcAddress(mod,"clReleaseEvent");
	DclEnqueueReadBuffer=(b5)GetProcAddress(mod,"clEnqueueReadBuffer");
	DclEnqueueWriteBuffer=(b6)GetProcAddress(mod,"clEnqueueWriteBuffer");
	DclReleaseMemObject=(b7)GetProcAddress(mod,"clReleaseMemObject");
	DclGetEventProfilingInfo=(b8)GetProcAddress(mod,"clGetEventProfilingInfo");
	DclEnqueueMapBuffer=(b9)GetProcAddress(mod,"clEnqueueMapBuffer");
	DclFinish=(c0)GetProcAddress(mod,"clFinish");
}
void transprint(u64 w){
	printf("%4lluGiB",w>>30);
	printf("%4lluMiB",(w>>20)&1023);
	printf("%4lluKiB",(w>>10)&1023);
	printf("%4lluB\n",w&1023);
}
cl_device_id Queryinfo(){
	cl_uint n;
	DclGetPlatformIDs(0,0,&n);
	cl_platform_id *plat=new cl_platform_id[n];
	DclGetPlatformIDs(n,plat,0);
	cl_device_id *dev=new cl_device_id[n],ca;
	u64 ret;
	char name[260],profile[260],version[260],vendor[260],*extensions;
	char devname[260];
	int Clock=0;
	u64 ComputeUnits=0,WorkItemDimensions=0,WorkGroupSize=0,*WorkItemSizes,LocalMemSize=0,ConstantBufferSize=0,GlobalMemSize=0;
	printf("%d platform found:\n\n",n);
	for(int i=0;i<n;++i){
		printf("Number %d:\n",i);
		DclGetPlatformInfo(plat[i],CL_PLATFORM_NAME   ,260,name   ,&ret);
		DclGetPlatformInfo(plat[i],CL_PLATFORM_PROFILE,260,profile,&ret);
		DclGetPlatformInfo(plat[i],CL_PLATFORM_VERSION,260,version,&ret);
		DclGetPlatformInfo(plat[i],CL_PLATFORM_VENDOR ,260,vendor ,&ret);
		DclGetPlatformInfo(plat[i],CL_PLATFORM_EXTENSIONS,0,0,&ret);
		extensions=new char[ret];
		DclGetPlatformInfo(plat[i],CL_PLATFORM_EXTENSIONS,ret,extensions,&ret);
		
		DclGetDeviceIDs(plat[i],CL_DEVICE_TYPE_GPU,1,dev+i,0);
		DclGetDeviceInfo(dev[i],CL_DEVICE_NAME                    ,260, devname           ,&ret);
		DclGetDeviceInfo(dev[i],CL_DEVICE_MAX_CLOCK_FREQUENCY     ,  4,&Clock             ,&ret);
		DclGetDeviceInfo(dev[i],CL_DEVICE_MAX_COMPUTE_UNITS       ,  8,&ComputeUnits      ,&ret);
		DclGetDeviceInfo(dev[i],CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,  8,&WorkItemDimensions,&ret);
		DclGetDeviceInfo(dev[i],CL_DEVICE_MAX_WORK_GROUP_SIZE     ,  8,&WorkGroupSize     ,&ret);
		WorkItemSizes=new u64[WorkItemDimensions];
		u64 k=sizeof(u64)*WorkItemDimensions;
		DclGetDeviceInfo(dev[i],CL_DEVICE_MAX_WORK_ITEM_SIZES     ,  k, WorkItemSizes     ,&ret);
		DclGetDeviceInfo(dev[i],CL_DEVICE_LOCAL_MEM_SIZE          ,  8,&LocalMemSize      ,&ret);
		DclGetDeviceInfo(dev[i],CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,  8,&ConstantBufferSize,&ret);
		DclGetDeviceInfo(dev[i],CL_DEVICE_GLOBAL_MEM_SIZE         ,  8,&GlobalMemSize     ,&ret);
		
		puts(name);puts(profile);puts(version);puts(vendor);puts(extensions);
		puts("");
		puts(devname);
		printf("ClockFrequency    :%dMHz\n",Clock);
		printf("ComputeUnits      :%llu\n",ComputeUnits);
		printf("WorkItemDimensions:%llu\n",WorkItemDimensions);
		printf("WorkGroupSize     :%llu\n",WorkGroupSize);
		printf("WorkItemSizes     :");
		for(int j=0;j<WorkItemDimensions;++j)printf("%llu ",WorkItemSizes[j]);
		puts("");
		printf("LocalMemSize      :");transprint(LocalMemSize);
		printf("ConstantBufferSize:");transprint(ConstantBufferSize);
		printf("GlobalMemSize     :");transprint(GlobalMemSize);
		puts("");
		delete extensions;
		delete WorkItemSizes;
	}
	int id=0;
	if(n>1){
		puts("Choose a device:");
		for(int i=0;i<n;++i){
			DclGetDeviceInfo(dev[i],CL_DEVICE_NAME,260,devname,&ret);
			printf("%2d: %s\n",i,devname);
		}
		scanf("%d",&id);
		if(id<0||id>=n)id=0;
	}
	ca=dev[id];
	delete[] plat;
	delete[] dev;
	return ca;
}
inline void printlog(cl_event event){
	cl_ulong t0,t1,t2,t3,t4;
	DclGetEventProfilingInfo(event,CL_PROFILING_COMMAND_QUEUED  ,sizeof(cl_ulong),&t0,0);
	DclGetEventProfilingInfo(event,CL_PROFILING_COMMAND_SUBMIT  ,sizeof(cl_ulong),&t1,0);
	DclGetEventProfilingInfo(event,CL_PROFILING_COMMAND_START   ,sizeof(cl_ulong),&t2,0);
	DclGetEventProfilingInfo(event,CL_PROFILING_COMMAND_END     ,sizeof(cl_ulong),&t3,0);
	DclGetEventProfilingInfo(event,CL_PROFILING_COMMAND_COMPLETE,sizeof(cl_ulong),&t4,0);
	// printf("Submit    time:%lf ms\n",(t1-t0)/1000000.0);
	printf("Execution time:%lf ms\n",(t3-t2)/1000000.0);
	// printf("Total     time:%lf ms\n",(t4-t0)/1000000.0);
}
char *LoadCLscript(const char *path,u64 *sz){
	u32 s,ret;
	HANDLE h=CreateFileA(path,GENERIC_READ,0,0,OPEN_EXISTING,0,0);
	if(h==INVALID_HANDLE_VALUE)return 0;
	s=GetFileSize(h,&ret);
	if(s>1048576)return 0;
	*sz=s;
	char *buf=new char[s];
	ReadFile(h,buf,s,&ret,0);
	CloseHandle(h);
	return buf;
}
cl_device_id dev;
cl_context context;
cl_command_queue queue;

struct source{
	cl_program p;
	bool Load(const char *srcpath){
		int ret;
		u64 sz;
		char *src=LoadCLscript(srcpath,&sz);
		p=DclCreateProgramWithSource(context,1,(const char**)&src,(const u64*)&sz,&ret);
		ret=DclBuildProgram(p,1,&dev,"-cl-std=CL3.0",0,0);
		if(ret!=CL_SUCCESS){
			u64 sz;
			DclGetProgramBuildInfo(p,dev,CL_PROGRAM_BUILD_LOG,0,0,&sz);
			char *log=new char[sz];
			DclGetProgramBuildInfo(p,dev,CL_PROGRAM_BUILD_LOG,sz,log,0);
			printf("Compile Failed!\n\n%s",log);
			return 1;
		}
		return 0;
	}
};
struct memory{
	cl_mem m;
	u64 len;
	void Buffin(u64 l,void *ptr){
		len=l;
		int ret;
		m=DclCreateBuffer(context,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,len,ptr,&ret);
	}
	void Buffout(u64 l){
		len=l;
		int ret;
		m=DclCreateBuffer(context,CL_MEM_READ_WRITE,len,0,&ret);
	}
	void Clear(){DclReleaseMemObject(m);}
	int Readint(int p){
		int ret;
		DclEnqueueReadBuffer(queue,m,CL_TRUE,p*sizeof(int),sizeof(int),&ret,0,0,0);
		return ret;
	}
	void Buffread(void *ptr,bool log=0){
		cl_event event;
		DclEnqueueReadBuffer(queue,m,CL_TRUE,0,len,ptr,0,0,&event);
		DclWaitForEvents(1,&event);
		if(log)printlog(event);
		DclReleaseEvent(event);
	}
	void Buffwrite(void *ptr,bool log=0){
		cl_event event;
		DclEnqueueWriteBuffer(queue,m,CL_TRUE,0,len,ptr,0,0,&event);
		DclWaitForEvents(1,&event);
		if(log)printlog(event);
		DclReleaseEvent(event);
	}
};
struct kernel{
	cl_kernel k;
	void Load(source pro,const char *entry){
		int ret;
		k=DclCreateKernel(pro.p,entry,&ret);
	}
	void Setptr(int cid,memory v){
		DclSetKernelArg(k,cid,sizeof(cl_mem),(void*)&(v.m));
	}
	void Setint(int cid,int v){
		DclSetKernelArg(k,cid,sizeof(int),(void*)&v);
	}
	void Setfloat(int cid,float v){
		DclSetKernelArg(k,cid,sizeof(float),(void*)&v);
	}
	void call2(u64 o0,u64 o1,u64 g0,u64 g1,bool log=0,bool wait=1){//dim=2,offset,global,item
		if(g0==0||g1==0)return;
		u64 offset[2]={o0,o1},group[2]={g0,g1};
		if(!wait){
			DclEnqueueNDRangeKernel(queue,k,2,offset,group,0,0,0,0);
			return;	
		}
		cl_event event;
		DclEnqueueNDRangeKernel(queue,k,2,offset,group,0,0,0,&event);
		DclWaitForEvents(1,&event);
		if(log)printlog(event);
		DclReleaseEvent(event);
	}
	void call1(u64 o,u64 g,bool log=0,bool wait=1){//dim=1,offset,global,item
		if(!g)return;
		if(!wait){
			DclEnqueueNDRangeKernel(queue,k,1,&o,&g,0,0,0,0);
			return;
		}
		cl_event event;
		DclEnqueueNDRangeKernel(queue,k,1,&o,&g,0,0,0,&event);
		DclWaitForEvents(1,&event);
		if(log)printlog(event);
		DclReleaseEvent(event);
	}
	static void wait(){
		DclFinish(queue);
	}
};
void InitialCL(){
	int ret;
	dev=Queryinfo();
	context=DclCreateContext(0,1,&dev,0,0,&ret);
	queue=DclCreateCommandQueue(context,dev,CL_QUEUE_PROFILING_ENABLE,&ret);
}

#endif
