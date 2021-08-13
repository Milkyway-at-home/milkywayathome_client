//
// Created by meyerc5 on 5/14/20.
//

#include <popt.h>
#include <stdio.h>
#include "nbody_lua.h"
#include <omp.h>
#include <nbody_likelihood.h>
#include <nbody_histogram.h>
#include "nbody_cl.h"
#include "milkyway_util.h"
#include "milkyway_cl.h"
#include "nbody.h"
#include "nbody_io.h"
#include "nbody_plain.h"
#include "milkyway_extra.h"
#include "nbody_types.h"
#include "nbody_tree.h"
#include "nbody_defaults.h"
#include "milkyway_git_version.h"



//static NBodyCtx _ctx = EMPTY_NBODYCTX;
//static NBodyState _st = EMPTY_NBODYSTATE;
static NBodyCtx Cctx = EMPTY_NBODYCTX;
static NBodyState Cst = EMPTY_NBODYSTATE;
static NBodyCtx Gctx = EMPTY_NBODYCTX;
static NBodyState Gst = EMPTY_NBODYSTATE;

int steps = 3000;

static void CLR(CLRequest* clr, const NBodyFlags* nbf){

    memset(clr, 0, sizeof(*clr));

    clr->platform = nbf->platform;
    clr->devNum = nbf->devNum;    NBodyStatus rc = NBODY_SUCCESS;
    clr->enableProfiling = TRUE;
}


NBodyState * runGPU(const NBodyFlags* nbf, const HistogramParams* hp){
    NBodyCtx* ctx = &Gctx;
    NBodyState* st = &Gst;
    CLRequest clr;
    lua_State* luaSt;
    st->checkpointResolved = "checkpoint.dat";
    CLR(&clr,nbf);
    int rc = nbInitCL(st, ctx, &clr);
    if(rc){
        mw_printf("OpenCL failed to initialize. Make sure OpenCL is supported on this device.\n");
        exit(-2);
    }
    nbSetup(ctx,st,nbf);
    ctx->nStep = steps;
    rc = nbInitNBodyStateCL(st, ctx);
    if(rc){
        exit(-1);
    }
    nbRunSystemCL(ctx,st);
    return st;
}

NBodyState* runCPU(const NBodyFlags* nbf, const HistogramParams* hp){
    NBodyCtx* ctx = &Cctx;
    NBodyState* st = &Cst;
    st->checkpointResolved = "checkpoint.dat";
    nbSetup(ctx,st,nbf);
    ctx->nStep = steps;
    nbRunSystemPlain(ctx,st,nbf);
    return st;
}

int main(int argc, char *argv[]){
    NBodyFlags nbf = EMPTY_NBODY_FLAGS;
    int rc = 0;
    omp_set_num_threads(omp_get_max_threads());
    steps = 3000;
    nbf.debugLuaLibs = 0;
    nbf.outputlbrCartesian = 1;
	
    nbf.inputFile = argv[1];
    nbf.checkpointPeriod = 0;
    //nbf.checkpointFileName = "checkpoint.dat";
    nbf.printHistogram = 1;
    nbf.seed = 1459;
    const char** a = mwCalloc(6,sizeof(char*));
    int i = 0;
    while(a[i] == NULL){
        a[i++] = mwCalloc(20,1);
    }
    nbf.numForwardedArgs = 6;
    strcpy(a[0],"0.1");
    strcpy(a[1],"1.0");
    strcpy(a[2],"0.2");
    strcpy(a[3],"0.2");
    strcpy(a[4],"12.0");
    strcpy(a[5],"0.2");
    a[6] = NULL;
    HistogramParams hp;
    NBodyLikelihoodMethod method;
    nbf.forwardedArgs = a;
    //Test 1 step
    
    nbGetLikelihoodInfo(&nbf,&hp,&method);
    
    nbf.outFileName = "test_outputGPU.txt";
    nbf.histoutFileName = "GPU_hist.dat";
    NBodyState G_ST = *runGPU(&nbf,&hp);
    
    nbf.outFileName = "test_outputCPU.txt";
    nbf.histoutFileName = "CPU_hist.dat";
    NBodyState  C_ST = *runCPU(&nbf,&hp);

    Body * p = C_ST.bodytab;
    Body * q = G_ST.bodytab;
    assert(C_ST.nbody == G_ST.nbody);
    double* data = mwMalloc(sizeof(double)*C_ST.nbody);
    double* dataAvg = mwMalloc(sizeof(double)*C_ST.nbody);
    while(p < C_ST.bodytab + C_ST.nbody && q < G_ST.bodytab + G_ST.nbody){
        double x_dif = X(Pos(p))-X(Pos(q));
        double y_dif = Y(Pos(p))-Y(Pos(q));
        double z_dif = Z(Pos(p))-Z(Pos(q));
        data[p-C_ST.bodytab] = mw_sqrt(sqr(x_dif)+sqr(y_dif)+sqr(z_dif));
        p++;q++;
    }
    for(int i = 0; i < C_ST.nbody; i++) {
        double Avgdistance = 0.0;
    	for(int j = 0; j < C_ST.nbody; j++) {
    	   if(i != j) {
    	       Body* first = C_ST.bodytab;
    	       first += i;
    	       Body* second = C_ST.bodytab;
    	       second += j;
    	       double x_dif = X(Pos(first))-X(Pos(second));
    	       double y_dif = Y(Pos(first))-Y(Pos(second));
               double z_dif = Z(Pos(first))-Z(Pos(second));
               Avgdistance += mw_sqrt(sqr(x_dif)+sqr(y_dif)+sqr(z_dif));
    	   }
    	}
    	Avgdistance /= (C_ST.nbody - 1);
    	dataAvg[i] = Avgdistance;
    }
    
    double sum = 0.0;
    double sumAvg = 0.0;
    for (int j = 0; j < C_ST.nbody; ++j) {
        sum += data[j];
        if(j + 1 < C_ST.nbody) {
           sumAvg += dataAvg[j];
        }
    }
    sum /= C_ST.nbody;
    sumAvg /= C_ST.nbody;
    mw_printf("The average distance between the GPU and the CPU results is: %.8lf\n",sum);
    mw_printf("The average distance between datapoints is: %.8lf\n",sumAvg);
    remove("checkpoint.dat");
    remove("CPU_hist.dat");
    remove("GPU_hist.dat");
    remove("test_outputCPU.txt");
    
    if(sum >= .15){
        mw_printf("Test failed\n");
        return -1;
    } else {
        mw_printf("Test passed\n");
        return 0;
    }
}
