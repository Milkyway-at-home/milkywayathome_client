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

NBodyState * runGPU(const NBodyFlags* nbf, const HistogramParams* hp){
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
    MainStruct* hist = nbCreateHistogram(ctx,st,hp);
    nbWriteHistogram(nbf->histoutFileName,ctx,st,hist);
    return st;
}

NBodyState* runCPU(const NBodyFlags* nbf, const HistogramParams* hp){
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
    return st;
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

    nbGetLikelihoodInfo(&nbf,&hp,&method);
    nbf.outFileName = "test_outputCPU.txt";
    nbf.histoutFileName = "CPU_hist.dat";
    NBodyState  C_ST = *runCPU(&nbf,&hp);

    nbf.outFileName = "test_outputGPU.txt";
    nbf.histoutFileName = "GPU_hist.dat";
    NBodyState G_ST = *runGPU(&nbf,&hp);
    Body * p = C_ST.bodytab;
    Body * q = G_ST.bodytab;
    assert(C_ST.nbody == G_ST.nbody);
    double* data = mwMalloc(sizeof(double)*C_ST.nbody);
    while(p < C_ST.bodytab + C_ST.nbody && q < G_ST.bodytab + G_ST.nbody){
        double x_dif = X(Pos(p))-X(Pos(q));
        double y_dif = Y(Pos(p))-Y(Pos(q));
        double z_dif = Z(Pos(p))-Z(Pos(q));
        data[p-C_ST.bodytab] = mw_sqrt(sqr(x_dif)+sqr(y_dif)+sqr(z_dif));
        p++;q++;
    }
    double sum = 0.0;
    for (int j = 0; j < C_ST.nbody; ++j) {
        sum += data[j];
    }
    sum /= C_ST.nbody;
    printf("The average distance between the GPU and the CPU results is: %.8lf\n",sum);


}