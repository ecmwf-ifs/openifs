/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifdef PAPI

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>  /* sleep(1) */
#include <pthread.h>
#include <signal.h>
#include <string.h>
#include <sched.h>
#include <papi.h>

struct eval {
    double val;
    int rank;
};

int rank; //rank
int numranks; //total number of ranks
char affinity[128];//core affinity
pid_t pid; //pid of process
char nodename[128]; //nodename
FILE *fp; //outputfil

#define NUM_EV 4

long_long papi_values[NUM_EV];
long_long start_usec_p;
long_long start_usec_r;

int init_error;

struct  thread_data{
    pthread_t parent_thread;
    float report_interval;
};

//global as it needs to live on also if master thread leaves MPI_Init...
struct thread_data td;


/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
static
char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
    char *ptr = str;
    int i, j, entry_made = 0;
    for (i = 0; i < CPU_SETSIZE; i++) {
        if (CPU_ISSET(i, mask)) {
            int run = 0;
            entry_made = 1;
            for (j = i + 1; j < CPU_SETSIZE; j++) {
                if (CPU_ISSET(j, mask)) run++;
                else break;
            }
            if (!run)
                sprintf(ptr, "%d,", i);
            else if (run == 1) {
                sprintf(ptr, "%d,%d,", i, i + 1);
                i++;
            } else {
                sprintf(ptr, "%d-%d,", i, i + run);
                i += run;
            }
            while (*ptr != 0) ptr++;
        }
    }
    ptr -= entry_made;
    *ptr = 0;
    return(str);
}

static
void getProcessIdentity(char *hostname,int hostnamesize,char *affinity,int affinitysize,int *rank,int *numranks,pid_t *pid){
    cpu_set_t coremask;

    (void)gethostname(hostname,hostnamesize);
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);
    cpuset_to_cstr(&coremask, affinity);
 
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
    MPI_Comm_size(MPI_COMM_WORLD,numranks);
        
    //pid
    *pid=getpid();
}

static
int parseMeminfo(const char *label){
    int i;
    int done=0;
    FILE *fp=fopen("/proc/meminfo","r");
    char line[120];
    char *subline;
    int val=-1;
    int label_length=strlen(label);
    
    while(!done && fp!=NULL){
      if (fgets(line,sizeof(line),fp)) {
        subline=strstr(line,label);
        if(subline!=NULL){
            val=atoi(subline+label_length);
            fclose(fp);
            return val;
        }
      }
    }
    fclose(fp);
    return val;
}

static
int report_init(int periodicreport){
    double rtime,ptime;
    int events[NUM_EV];
    int rc, num;
    static int first_time = 1;
    
    /* open file if we use files and write headers*/
    if(periodicreport){
        char fname[40];
        sprintf(fname,"perfmon_report_%d.dat",rank);
        fp=fopen(fname,"w");
        if(fp==NULL){
            return 2;
        }
        //header
        fprintf(fp,"# rank: %d nodename: %s core-affinity: %s \n",rank,nodename,affinity);
        fprintf(fp,"#rtime ptime gflops l1-hit memusage(MB) freeMem(MB)\n");
    }

    if (!first_time) return 0;
    first_time = 0;
    
    //get papi info, first time it intializes PAPI counters/library
    events[0]=PAPI_L1_DCM;
    events[1]=PAPI_L1_DCH;
    events[2]=PAPI_FP_OPS;
    events[3]=PAPI_TOT_INS;

    rc = (num = PAPI_num_counters());
    if (rc != PAPI_OK) {
      PAPI_perror(rc, "PAPI_num_counters", strlen("PAPI_num_counters"));
    }

    //fprintf(stderr,"PAPI_num_counters = %d\n",num);
    
    rc = PAPI_start_counters(events, NUM_EV);
    if (rc != PAPI_OK) {
      return rc;
    }
    start_usec_r=PAPI_get_real_usec();
    start_usec_p=PAPI_get_virt_usec();
    return 0;
}

static
void report_periodic(){
    double rtime,ptime;
    static double prevrtime=0.0;
    double freemem;
    double memuse;
    double gflops,l1hitratio;
    PAPI_dmem_info_t dmem;
    int i;
    //get papi info,
    long_long end_usec_r,end_usec_p;
    long_long prev_values[NUM_EV];

    end_usec_r = PAPI_get_real_usec();
    end_usec_p = PAPI_get_virt_usec();
        
    rtime=(double)(end_usec_r-start_usec_r)/1e6;
    ptime=(double)(end_usec_p-start_usec_p)/1e6;
    
    //get memdata
    PAPI_get_dmem_info(&dmem);
    for(i=0;i<NUM_EV;i++)
        prev_values[i]=papi_values[i];
    PAPI_accum_counters(papi_values, NUM_EV);
    
    l1hitratio=100.0*(double)(papi_values[1]-prev_values[1])
        /(papi_values[0]+papi_values[1]-prev_values[0]-prev_values[1]);
    gflops=(double)(papi_values[2]-prev_values[2])/(rtime-prevrtime)/1e9;
    

    memuse=(double)dmem.size/1.0e3;
    //get free memory
    freemem=(double)parseMeminfo("MemFree:")/1.0e3;
    
    fprintf(fp," %f %f %f  %f %f %f\n",
            rtime,ptime,gflops,l1hitratio,memuse,freemem);
    fflush(fp); //we should flush more seldomly?
    prevrtime=rtime;
}

static
void coll_print(FILE *fp, const char *label,double val,int print_aggregate,MPI_Comm Comm){
    struct eval in;
    struct eval out;
    double sum;
    in.val=val;
    in.rank=rank;
    MPI_Reduce(&val,&sum,1,MPI_DOUBLE,MPI_SUM,0,Comm);
    if(rank==0){
        if(print_aggregate)
            fprintf(fp,"#%19s %14.3f %10.3f ",label,sum,sum/numranks);
        else
            fprintf(fp,"#%19s                %10.3f ",label,sum/numranks);
    }

    MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MINLOC,0,Comm);
    if(rank==0){
        fprintf(fp,"%4d %10.3f ",out.rank,out.val);
    }
    MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,0,Comm);
    if(rank==0){
        fprintf(fp,"%4d %10.3f\n",out.rank,out.val);
    }
}

static
void report_final(FILE *fp, const MPI_Comm *comm){
    double rtime,ptime;
    double memuse,peakmem;
    double gflops,avegflops;
    double gflop_opers;
    PAPI_dmem_info_t dmem;
    double val;
    struct eval in;
    struct eval out;
    int error=0;
    double l1hitratio;
    long_long end_usec_p;
    long_long end_usec_r;
    MPI_Comm Comm = comm ? *comm : MPI_COMM_WORLD;
    
    //get papi info, first time it intializes PAPI counters
    end_usec_r = PAPI_get_real_usec();
    end_usec_p = PAPI_get_virt_usec();
    
    if(PAPI_accum_counters(papi_values, NUM_EV) != PAPI_OK)
        error++;

    if(PAPI_get_dmem_info(&dmem)!=PAPI_OK)
        error++;
    
    rtime=(double)(end_usec_r-start_usec_r)/1e6;
    ptime=(double)(end_usec_p-start_usec_p)/1e6;
    l1hitratio=100.0*(double)papi_values[1]/(papi_values[0]+papi_values[1]);
    avegflops=(double)papi_values[2]/rtime/1e9;
    gflop_opers = (double)papi_values[2]/1e9;
    
    if(rank==0 ){
        fprintf(fp,"####### CSC PERFMON REPORT  \n");
        fprintf(fp,"# MPI tasks   %d\n",numranks);
        fprintf(fp,"#                        aggregated    average    min(rank/val)   max(rank/val) \n");
    }
    coll_print(fp,"Real time (s)",rtime,1,Comm);
    coll_print(fp,"Process time (s)",ptime,1,Comm);
    coll_print(fp,"Flops (GFlop/s)",avegflops,1,Comm);
    coll_print(fp,"Flp-opers (10^9)",gflop_opers,1,Comm);
    coll_print(fp,"L1 hit ratio (%)",l1hitratio,0,Comm);
    coll_print(fp,"Peak mem size (MB)", (double)dmem.peak/1.0e3,0,Comm );
    coll_print(fp,"Peak resident (MB)", (double)dmem.high_water_mark/1.0e3 ,0,Comm);
    if(rank==0)  fprintf(fp,"#######           \n");
    fflush(fp);
}

static
void *thread_worker(void *threadarg){
    //we could of course also directly read global td
    struct thread_data *tdloc= (struct thread_data *) threadarg;
//if zero then we do not report data
    if(tdloc->report_interval>0){ 
        //find out how many total seconds and mikroseconds to sleep
        int sec=(int)tdloc->report_interval;  
        int usec=(int)((tdloc->report_interval-sec)*1.0e6);
        //do not allow busy loop  when interval less than usec
        while(sec+usec>0){
            if(sec>0)  sleep(sec);
            if(usec>0) usleep(usec);
            //send signal to parent thread that it should report flops etc.
            pthread_kill(tdloc->parent_thread,SIGUSR1);
        }

    }
    
    return NULL;
}


static
void common_inits()
{
    pthread_t t;
    pthread_attr_t thread_attr;
    int thread_id,thread_create_return;
    int temp;
    char *envvar;
        
    //initialize global values for process identity
    getProcessIdentity(nodename,sizeof(nodename),affinity,sizeof(affinity),&rank,&numranks,&pid);
    
    //init parameters
  

    envvar=getenv("PERFMON_INTERVAL");
    if(envvar==NULL){
        td.report_interval=-1; //default never
    }
    else {
        td.report_interval=atof(envvar);
    }

  

    /* if report interval is larger or equal to 10ms then do "unsafe"
       periodic reporting stuff; start up signal handling and launch
       thread
       
    */
    
    if(td.report_interval>=0.009999) {
        //initialize PAPI counters with periodic reporting
        init_error=report_init(1);
        if(!init_error){
            //print PAPI counters when receiving signal USR1
            signal(SIGUSR1,report_periodic); 
            
            //launch sampling thread 
            td.parent_thread=pthread_self();
            temp=pthread_create(&t,NULL,thread_worker,(void *)&td);
        }
    }
    else {
       //initialize PAPI counters without periodic reporting
      fprintf(stderr,"Calling report_init(0)\n");
      init_error=report_init(0);
      if (init_error) fprintf(stderr,
			      "Unable to init PAPI counters (init_error=%d) : %s\n",
			      init_error,
			      PAPI_strerror(init_error));
    }
}


int __wrap_MPI_Init(int *argc, char ***argv)
{
  //call true MPI_Init
  int ret= __real_MPI_Init(argc, argv);
  common_inits();
  return ret;
}


int __wrap_MPI_Init_thread(int *argc, char ***argv, int required, int *provided)
{
  //call true MPI_Init_thread
  int ret= __real_MPI_Init(argc, argv, required, provided);
  common_inits();
  return ret;
}



int __wrap_MPI_Finalize()
{

    //only print if startup was without errors
    if(!init_error)
        report_final(stdout,NULL);
    return __real_MPI_Finalize();
}


void csc_perfmon_begin_  () { common_inits();       }
void csc_perfmon_begin__ () { csc_perfmon_begin_(); }
void csc_perfmon_begin   () { csc_perfmon_begin_(); }
void CSC_PERFMON_BEGIN   () { csc_perfmon_begin_(); }

void csc_perfmon_end_  ()   { if(!init_error) report_final(stdout,NULL); }
void csc_perfmon_end__ ()   { csc_perfmon_end_(); }
void csc_perfmon_end   ()   { csc_perfmon_end_(); }
void CSC_PERFMON_END   ()   { csc_perfmon_end_(); }

#else

void dummy_Opfla_PERFMOM() { }

#endif
