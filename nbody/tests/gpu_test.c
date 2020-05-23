//
// Created by meyerc5 on 5/14/20.
//

#include <popt.h>
#include "nbody_lua.h"
#include <omp.h>
#include <nbody_likelihood.h>
#include <nbody_histogram.h>
#include "nbody_cl.h"
#include "milkyway_util.h"
#include "nbody.h"
#include "nbody_io.h"
#include "nbody_plain.h"
#include "milkyway_extra.h"
#include "nbody_types.h"
#include "nbody_tree.h"
#include "nbody_defaults.h"
#include "milkyway_git_version.h"



static NBodyCtx _ctx = EMPTY_NBODYCTX;
static NBodyState _st = EMPTY_NBODYSTATE;
static NBodyCtx Cctx = EMPTY_NBODYCTX;
static NBodyState Cst = EMPTY_NBODYSTATE;
static NBodyCtx Gctx = EMPTY_NBODYCTX;
static NBodyState Gst = EMPTY_NBODYSTATE;

int steps = 2;

static void CLR(CLRequest* clr, const NBodyFlags* nbf){

    memset(clr, 0, sizeof(*clr));

    clr->platform = nbf->platform;
    clr->devNum = nbf->devNum;    NBodyStatus rc = NBODY_SUCCESS;
    clr->enableProfiling = TRUE;
}
int do_nothing(const NBodyFlags* nbf){
    NBodyCtx* ctx = &_ctx;
    NBodyState* st = &_st;
    st->checkpointResolved = "checkpoint.dat";
    nbSetup(ctx,st,nbf);
    nbMakeTree(ctx,st);
    nbWriteBodies(ctx, st, nbf);
}

int runGPU(const NBodyFlags* nbf, const HistogramParams* hp){
    NBodyCtx* ctx = &Gctx;
    NBodyState* st = &Gst;
    CLRequest clr;
    lua_State* luaSt;
    st->checkpointResolved = "checkpoint.dat";
    for(int i = 0;i < 6;i++){
        printf("%s\n",nbf->forwardedArgs[i]);
    }
    CLR(&clr,nbf);
    int rc = nbInitCL(st, ctx, &clr);
    if(rc){
        exit(-1);
    }
    nbSetup(ctx,st,nbf);
//    ctx->useQuad = TRUE;
//    ctx->criterion = TreeCode;
//    ctx->nStep = steps;
    printf("%d\n",ctx->nStep);
    rc = nbInitNBodyStateCL(st, ctx);

    if(rc){
        exit(-1);
    }
    nbRunSystemCL(ctx,st);
    nbWriteBodies(ctx, st, nbf);
    MainStruct* hist = nbCreateHistogram(ctx,st,hp);
    nbWriteHistogram(nbf->histoutFileName,ctx,st,hist);
    return 0;
}

int runCPU(const NBodyFlags* nbf, const HistogramParams* hp){
    NBodyCtx* ctx = &Cctx;
    NBodyState* st = &Cst;
    st->checkpointResolved = "checkpoint.dat";
    nbSetup(ctx,st,nbf);
//    ctx->useQuad = TRUE;
//    ctx->criterion = TreeCode;
//    ctx->nStep = steps;
    printf("%d\n",ctx->nStep);
    nbRunSystemPlain(ctx,st,nbf);
    nbWriteBodies(ctx, st, nbf);
    MainStruct* hist = nbCreateHistogram(ctx,st,hp);
    nbWriteHistogram(nbf->histoutFileName,ctx,st,hist);
    return 0;
}

int main(int argc, char* argv[]){
    NBodyFlags nbf = EMPTY_NBODY_FLAGS;
    int rc = 0;
    omp_set_num_threads(16);
    steps = 1000;
    nbf.debugLuaLibs = 0;
    nbf.outputlbrCartesian = 1;
    nbf.inputFile = "four_developers.lua";
    nbf.checkpointFileName = "checkpoint.dat";
    nbf.printHistogram = 1;
    nbf.seed = 1459;
    const char** a = mwCalloc(6,sizeof(char*));
    int i = 0;
    while(a[i] == NULL){
        a[i++] = mwCalloc(20,1);
    }
    nbf.numForwardedArgs = 6;
    strcpy(a[0],"4.0");
    strcpy(a[1],"1.0");
    strcpy(a[2],"0.2");
    strcpy(a[3],"0.2");
    strcpy(a[4],"12.0");
    strcpy(a[5],"0.2");
    a[6] = NULL;
    HistogramParams hp;
    NBodyLikelihoodMethod method;
    nbf.forwardedArgs = a;

    nbGetLikelihoodInfo(&nbf,&hp,&method);
    nbf.outFileName = "test_outputCPU.txt";
    nbf.histoutFileName = "CPU_hist.dat";
    runCPU(&nbf,&hp);
    nbf.outFileName = "test_outputGPU.txt";
    nbf.histoutFileName = "GPU_hist.dat";
    runGPU(&nbf,&hp);

//    FILE * cpu_results= fopen("test_outputCPU.txt","r");
//    char buff[256];
//    while(fscanf(cpu_results,"")
}