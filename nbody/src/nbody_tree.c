/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
   Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "nbody_priv.h"
#include "nbody_tree.h"
#include <time.h>

#include <lua.h>
#include <lauxlib.h>
#include "nbody_lua_types.h"
#include "milkyway_util.h"

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

/* subIndex: compute subcell index for body p in cell q. */

/* Find k nearest neighbors (k=32) for body p using the NBodyTree.
 * Results are written into out[], sorted by distance (closest first).
 * Returns the number of neighbors actually found (≤32).
 */
#define K 32

/* Priority queue struct for nearest neighbor calculations */
typedef struct {
    Body* neighbors[K];
    real dists[K];
    int size;           /* Current number of neighbors */
    real worstDist;     /* Current worst distance in the heap */
} NeighborPQ;

/* create  priority queue */
static void initPQ(NeighborPQ* pq) {
    for (int i = 0; i < K; i++) {
        pq->neighbors[i] = NULL;
        pq->dists[i] = HUGE_VAL;
    }
    pq->size = 0;
    pq->worstDist = HUGE_VAL;
}

/* Insert neighbor into priority queue if it's closer than current worst */
static inline void insertNeighbor(NeighborPQ* pq, Body* b, real dist2) {
    /* Skip dark matter or invalid bodies */
    if (b == NULL || b->bodynode.type > -1) {
        return;
    }
    
    /* If queue isn't full, just add it */
    if (pq->size < K) {
        int i = pq->size++;
        pq->neighbors[i] = b;
        pq->dists[i] = dist2;
        
        /* Bubble up */
        while (i > 0) {
            int parent = (i - 1) / 2;
            if (pq->dists[parent] >= pq->dists[i]) break;
            
            /* Swap */
            Body* tempB = pq->neighbors[i];
            pq->neighbors[i] = pq->neighbors[parent];
            pq->neighbors[parent] = tempB;
            
            real tempD = pq->dists[i];
            pq->dists[i] = pq->dists[parent];
            pq->dists[parent] = tempD;
            
            i = parent;
        }
        
        if (pq->size == K) {
            pq->worstDist = pq->dists[0]; /* Max-heap, root is largest */
        }
        return;
    }
    
    /* Queue is full, check if new neighbor is better than worst */
    if (dist2 >= pq->worstDist) {
        return;
    }
    
    /* Replace worst neighbor (root of max-heap) */
    pq->neighbors[0] = b;
    pq->dists[0] = dist2;
    
    /* Heapify down */
    int i = 0;
    while (1) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        int largest = i;
        
        if (left < K && pq->dists[left] > pq->dists[largest]) {
            largest = left;
        }
        if (right < K && pq->dists[right] > pq->dists[largest]) {
            largest = right;
        }
        
        if (largest == i) break;
        
        /* Swap */
        Body* tempB = pq->neighbors[i];
        pq->neighbors[i] = pq->neighbors[largest];
        pq->neighbors[largest] = tempB;
        
        real tempD = pq->dists[i];
        pq->dists[i] = pq->dists[largest];
        pq->dists[largest] = tempD;
        
        i = largest;
    }
    
    pq->worstDist = pq->dists[0];
}

/* Sort neighbors by distance (ascending) using selection sort on small array */
static void sortNeighbors(NeighborPQ* pq) {
    for (int i = 0; i < pq->size - 1; i++) {
        int minIdx = i;
        for (int j = i + 1; j < pq->size; j++) {
            if (pq->dists[j] < pq->dists[minIdx]) {
                minIdx = j;
            }
        }
        if (minIdx != i) {
            /* Swap */
            Body* tempB = pq->neighbors[i];
            pq->neighbors[i] = pq->neighbors[minIdx];
            pq->neighbors[minIdx] = tempB;
            
            real tempD = pq->dists[i];
            pq->dists[i] = pq->dists[minIdx];
            pq->dists[minIdx] = tempD;
        }
    }
}

/* Optimized bounding box distance calculation */
static inline real boxDistanceSq(const mwvector* pPos, const mwvector* cPos, real halfSize) {
    real dx = X(*pPos) - X(*cPos);
    real dy = Y(*pPos) - Y(*cPos);
    real dz = Z(*pPos) - Z(*cPos);
    
    dx = (dx > 0) ? mw_fmax(0.0, dx - halfSize) : mw_fmax(0.0, -dx - halfSize);
    dy = (dy > 0) ? mw_fmax(0.0, dy - halfSize) : mw_fmax(0.0, -dy - halfSize);
    dz = (dz > 0) ? mw_fmax(0.0, dz - halfSize) : mw_fmax(0.0, -dz - halfSize);
    
    return dx*dx + dy*dy + dz*dz;
}

/* Recursive search with better pruning */
static void searchNodeOptimized(const NBodyTree* t, const NBodyNode* node, 
                                const mwvector* pPos, NeighborPQ* pq, real nodeSize) {
    if (node == NULL) return;
    
    /* Early prune using bounding sphere of node */
    if (isCell(node)) {
        const NBodyCell* c = (const NBodyCell*) node;
        real boxDist2 = boxDistanceSq(pPos, &Pos(c), nodeSize * 0.5f);
        
        if (boxDist2 >= pq->worstDist) {
            return; /* Entire node is too far */
        }
        
        /* Search children in order of increasing minimum distance */
        const NBodyNode* children[NSUB];
        real childDists[NSUB];
        int childCount = 0;
        
        /* Collect valid children and compute their distances */
        for (int i = 0; i < NSUB; i++) {
            NBodyNode* child = Subp(c)[i];
            if (child) {
                children[childCount] = child;
                if (isCell(child)) {
                    /* For cells, use bounding box distance */
                    childDists[childCount] = boxDistanceSq(pPos, &Pos(child), nodeSize * 0.25f);
                } else {
                    /* For bodies, compute actual distance */
                    real dx = X(Pos(child)) - X(*pPos);
                    real dy = Y(Pos(child)) - Y(*pPos);
                    real dz = Z(Pos(child)) - Z(*pPos);
                    childDists[childCount] = dx*dx + dy*dy + dz*dz;
                }
                childCount++;
            }
        }
        
        /* Sort children by distance (simple insertion sort for small NSUB=8) */
        for (int i = 1; i < childCount; i++) {
            int j = i;
            while (j > 0 && childDists[j-1] > childDists[j]) {
                /* Swap distances */
                real tempDist = childDists[j-1];
                childDists[j-1] = childDists[j];
                childDists[j] = tempDist;
                
                /* Swap children */
                const NBodyNode* tempChild = children[j-1];
                children[j-1] = children[j];
                children[j] = tempChild;
                j--;
            }
        }
        
        /* Search closest children first (better pruning) */
        for (int i = 0; i < childCount; i++) {
            if (childDists[i] >= pq->worstDist) {
                break; /* Remaining children are even farther */
            }
            searchNodeOptimized(t, children[i], pPos, pq, nodeSize * 0.5f);
        }
    } 
    else if (isBody(node)) {
        const Body* b = (const Body*) node;
        
        /* Skip self and non-dark matter */
        if (b->bodynode.type > -1) {
            return;
        }
        
        /* Compute squared distance */
        real dx = X(Pos(b)) - X(*pPos);
        real dy = Y(Pos(b)) - Y(*pPos);
        real dz = Z(Pos(b)) - Z(*pPos);
        real dist2 = dx*dx + dy*dy + dz*dz;
        if (dist2 > 0.1) {
               return;   // outside interaction range
    }

        insertNeighbor(pq, (Body*)b, dist2);
    }
}

/* Main function to find nearest neighbors */
static int nbNearestNeighborsOptimized(const NBodyTree* t, const Body* p, Body* out[K]) {
    NeighborPQ pq;
    mwvector pPos = Pos(p);
    
    initPQ(&pq);
    
    /* Start recursive search from root */
    searchNodeOptimized(t, (const NBodyNode*)t->root, &pPos, &pq, t->rsize);
    
    /* Sort neighbors by distance */
    sortNeighbors(&pq);
    
    /* Copy to output array */
    for (int i = 0; i < pq.size; i++) {
        out[i] = pq.neighbors[i];
    }
    
    return pq.size;
}

/* improved print functions */
int nbPrintNearestNeighborsOptimized(const NBodyTree* t, const Body* p, int particle_id, 
                                     FILE* fdata, FILE* fvel, FILE* fdens, FILE* fhalo, const NBodyCtx* ctx) {
    Body* neighbors[K] = {NULL};
    int found = nbNearestNeighborsOptimized(t, p, neighbors);
    
    if (found == 0) {
        return 0;
    }
    
    const mwvector pos0 = Pos(p);
    
    /* Compute max distance in one pass */
    real maxDistanceSq = 0.0;
    if(found<32){
	maxDistanceSq = 0.1;
    }
    for (int i = 0; i < found; i++) {

        if (neighbors[i] == NULL) continue;
        
        mwvector npos = Pos(neighbors[i]);
        real dx = X(npos) - X(pos0);
        real dy = Y(npos) - Y(pos0);
        real dz = Z(npos) - Z(pos0);
        real dist2 = dx*dx + dy*dy + dz*dz;
        
        if (dist2 > maxDistanceSq) {
            maxDistanceSq = dist2;
        }
    }
    
    real maxDistance = mw_sqrt(maxDistanceSq);
    real volume = (4.0f / 3.0f) * 3.14159265358979323846f * maxDistance * maxDistanceSq;
    //real volume = (4.0f/ 3.0f) * 3.14159265358979323846f * maxDistance * maxDistance;
    /* Precompute mass of one neighbor (assuming all have same mass) */
    real neighborMass = Mass(neighbors[0]);
    real total_mass = found * neighborMass;
    real local_density = total_mass / volume;
    real radius = mw_sqrt(X(Pos(p))*X(Pos(p)) + Y(Pos(p))*Y(Pos(p)) + Z(Pos(p))*Z(Pos(p)));
    //real expected_density = get_density(luaSt->comp,mw_sqrt(X(Pos(p))*X(Pos(p)) + Y(Pos(p))*Y(Pos(p)) + Z(Pos(p))*Z(Pos(p))));
    //real difference = expected_density/local_density;
    //mw_printf("density difference: %d\n",difference);
    /* Compute probability */
    real total_probability = 1.0f;
    //mw_printf("%s\n",ctx->pot.sphere[0].type);
    mw_printf("Potential: %s\n", showPotential(&ctx->pot));
mw_printf("pos0: %f %f %f\n", pos0.x, pos0.y, pos0.z);

    const real vhalo  = 74.61;
    const real q  = 1;
    const real a  = 12;
    const real R2 = mw_pow(X(pos0),2.0) + mw_pow(Y(pos0),2.0);

    real numer = (2.0*q*q + 1.0)*a*a + R2 + mw_pow(Z(pos0),2.0)*(2.0-1.0/q/q);
    real denom = q*q*mw_pow(R2 + a*a + mw_pow(Z(pos0)/q,2.0),2.0);

    real rho =  vhalo*vhalo*numer/2.0/3.1415926535/denom;


    //real rho = nbExtDensity(&ctx->pot,pos0,100);
    mw_printf("rho: %f\n",rho);
    real v = mw_sqrt(X(p->vel)*X(p->vel) + Y(p->vel)*Y(p->vel) + Z(p->vel)*Z(p->vel));
    //real sigma = 0.00010220187;
    real sigma=102.27;
    real pi = 3.1415926535;
    real term1 = sqrt(2.0 / pi) * sigma * exp(-(v*v) / (2.0 * sigma * sigma));
    real term2 = v * (1.0 + (sigma * sigma) / (v * v)) * erf(v / (sqrt(2.0) * sigma));
    real v_rel = term1*term2;
    const real cross_section = 0.00004639155;
    const mwvector p_vel = p->vel;
    real halo_probability = v_rel * cross_section * rho * 0.004;
    mw_printf("halo prob: %.6f %.6f %.6f",v_rel, cross_section, rho);
    for (int i = 0; i < found; i++) {
        if (neighbors[i] == NULL) continue;
        
        const Body* nb = neighbors[i];
        const mwvector nb_vel = nb->vel;
        
        real dvx = X(nb_vel) - X(p_vel);
        real dvy = Y(nb_vel) - Y(p_vel);
        real dvz = Z(nb_vel) - Z(p_vel);
        
        /* Use faster approximation for sqrt if precision allows */
        real rel_speed = mw_sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        real probability = cross_section * rel_speed * 0.004;
        
        total_probability *= (1.0f - probability);
        
        /* Write velocity if file is open */
        if (fvel) {
            real speed = mw_sqrt(mw_pow(nb_vel.x, 2) + mw_pow(nb_vel.y, 2) + mw_pow(nb_vel.z, 2));
            fprintf(fvel, "%.6f\n", speed);
        }
    }
    //total_probability *= (1.0f - halo_probability); //include probability of scattering off of MW halo 
    total_probability = 1.0f - total_probability;
    
    /* Write data if files are open */
    if (fdata) {
        fprintf(fdata, "%.6f\n", total_probability);
    }
    
    if (fdens) {
        fprintf(fdens, "%.6f %.6f\n", local_density * 222288.47f, radius);
    } 
    if (fhalo) {
	fprintf(fhalo,"%.6f\n", halo_probability);
    }
    /* Random scattering decision */
    /* Commented out as per original, but kept structure */
    /*
    double r = (double)rand() / RAND_MAX;
    return (r < total_probability) ? 1 : 0;
    */
    return 0; /* Or return scattering result if needed */
}

static inline int nbSubIndex(Body* p, NBodyCell* q)
{
    int ind = 0;

    /* accumulate subcell index */
    /* loop over dimensions */
    if (X(Pos(q)) <= X(Pos(p)))     /* if beyond midpoint */
        ind += NSUB >> (0 + 1);     /* skip over subcells */

    if (Y(Pos(q)) <= Y(Pos(p)))
        ind += NSUB >> (1 + 1);

    if (Z(Pos(q)) <= Z(Pos(p)))
        ind += NSUB >> (2 + 1);

    return ind;
}

static inline void nbIncAddNBodyQuadMatrix(NBodyQuadMatrix* RESTRICT a, NBodyQuadMatrix* RESTRICT b)
{
    a->xx += b->xx;
    a->xy += b->xy;
    a->xz += b->xz;

    a->yy += b->yy;
    a->yz += b->yz;

    a->zz += b->zz;
}

/* hackQuad: descend tree, evaluating quadrupole moments.  Note that this
 * routine is coded so that the Subp() and Quad() components of a cell can
 * share the same memory locations.
 */
static void hackQuad(NBodyCell* p)
{
    unsigned int ndesc, i;
    NBodyNode* desc[NSUB];
    NBodyNode* q;
    mwvector dr;
    real drsq;
    NBodyQuadMatrix quad = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    ndesc = 0;                                  /* count occupied subnodes  */
    for (i = 0; i < NSUB; ++i)                  /* loop over all subnodes   */
    {
        if (Subp(p)[i] != NULL)                 /* if this one's occupied   */
        {
            desc[ndesc++] = Subp(p)[i];         /* copy it to safety        */
        }
    }

    for (i = 0; i < ndesc; ++i)                 /* loop over real subnodes  */
    {
        q = desc[i];                            /* access each one in turn  */
        if (isCell(q))                          /* if it's also a cell      */
        {
            hackQuad((NBodyCell*) q);           /* then process it first    */
        }

        dr = mw_subv(Pos(q), Pos(p));           /* find displacement vect.  */
        drsq = mw_sqrv(dr);                     /* and dot prod. (dr . dr)  */

        /* Outer product scaled by 3, then subtract drsq off the
         * diagonal to form quad moment*/
        {
            real m = Mass(q);   /* from CM of subnode */

            quad.xx = m * (3.0 * (X(dr) * X(dr)) - drsq);
            quad.xy = m * (3.0 * (X(dr) * Y(dr)));
            quad.xz = m * (3.0 * (X(dr) * Z(dr)));

            quad.yy = m * (3.0 * (Y(dr) * Y(dr)) - drsq);
            quad.yz = m * (3.0 * (Y(dr) * Z(dr)));

            quad.zz = m * (3.0 * (Z(dr) * Z(dr)) - drsq);
        }

        if (isCell(q)) /* if subnode is cell       */
        {
            nbIncAddNBodyQuadMatrix(&quad, &Quad(q));     /* then include its moment  */
        }

        nbIncAddNBodyQuadMatrix(&Quad(p), &quad); /* increment moment of cell */
    }
}


/* threadTree: do a recursive treewalk starting from node p,
 * with next stop n, installing Next and More links.
 */
static void threadTree(NBodyNode* p, NBodyNode* n)
{
    unsigned int ndesc, i;
    NBodyNode* desc[NSUB+1];

    Next(p) = n;                                /* link to next node */
    if (isCell(p))                              /* any children to thread? */
    {
        ndesc = 0;                              /* count extant children */
        for (i = 0; i < NSUB; ++i)              /* loop over subnodes */
        {
            if (Subp(p)[i] != NULL)             /* found a live one? */
            {
                desc[ndesc++] = Subp(p)[i];     /* store in table */
            }
        }
        More(p) = desc[0];                      /* link to first child */
        desc[ndesc] = n;                        /* end table with next */
        for (i = 0; i < ndesc; i++)             /* loop over children */
        {
            threadTree(desc[i], desc[i + 1]);     /* thread each w/ next */
        }
    }
}

/* expandBox: find range of coordinate values (with respect to root)
 * and expand root cell to fit. The size is doubled at each step to
 * take advantage of exact representation of powers of two.
 */
static void expandBox(NBodyTree* t, const Body* btab, int nbody)
{
    real xyzmax;
    const Body* p;
    const NBodyCell* root = t->root;

    assert(t->rsize > 0.0);

    xyzmax = 0.0;
    for (p = btab; p < btab + nbody; ++p)
    {
        xyzmax = mw_fmax(xyzmax, mw_abs(X(Pos(p)) - X(Pos(root))));
        xyzmax = mw_fmax(xyzmax, mw_abs(Y(Pos(p)) - Y(Pos(root))));
        xyzmax = mw_fmax(xyzmax, mw_abs(Z(Pos(p)) - Z(Pos(root))));
    }

    while (t->rsize < 2.0 * xyzmax)
    {
        t->rsize *= 2.0;
    }
}

/* makecell: return pointer to free cell. */
static NBodyCell* nbMakeCell(NBodyState* st, NBodyTree* t)
{
    NBodyCell* c;

    if (st->freeCell == NULL)                   /* no free cells left? */
    {
        c = (NBodyCell*) mwMallocA(sizeof(*c)); /* allocate a new one */
    }
    else                                        /* use existing free cell */
    {
        c = (NBodyCell*) st->freeCell;          /* take one on front */
        st->freeCell = Next(c);                 /* go on to next one */
    }
    Type(c) = CELL(0);                          /* initialize cell type */
    More(c) = NULL;
    memset(&c->stuff, 0, sizeof(c->stuff));     /* empty sub cells */
    t->cellUsed++;                              /* count one more cell */
    return c;
}

/* reclaim cells in tree, prepare to build new one. */
static void nbNewTree(NBodyState* st, NBodyTree* t)
{
    NBodyNode* p = (NBodyNode*) t->root;              /* start with the root */

    while (p != NULL)                       /* loop scanning tree */
    {
        if (isCell(p))                      /* found cell to free? */
        {
            Next(p) = st->freeCell;         /* link to front of */
            st->freeCell = p;               /* ...existing list */
            p = More(p);                    /* scan down tree */
        }
        else                                /* skip over bodies */
        {
            p = Next(p);                    /* go on to next */
        }
    }

    t->cellUsed = 0;   /* init count of cells, levels */
    t->maxDepth = 0;

    t->root = nbMakeCell(st, t);      /* allocate the root cell */
    mw_zerov(Pos(t->root));           /* initialize the midpoint */
}


ALWAYS_INLINE
static inline real calcOffset(real pPos, real qPos, real qsize)
{
    /* offset from parent */
    return qPos + 0.25 * (pPos < qPos ? -qsize : qsize);
}

ALWAYS_INLINE
static inline void nbInitMidpoint(NBodyCell* c, const Body* p, const NBodyCell* q, real qsize)
{
    X(Pos(c)) = calcOffset(X(Pos(p)), X(Pos(q)), qsize);
    Y(Pos(c)) = calcOffset(Y(Pos(p)), Y(Pos(q)), qsize);
    Z(Pos(c)) = calcOffset(Z(Pos(p)), Z(Pos(q)), qsize);
}

/* loadBody: descend tree and insert body p in appropriate place. */
static void nbLoadBody(NBodyState* st, NBodyTree* t, Body* p)
{
    NBodyCell* q;
    NBodyCell* c;
    size_t qind;
    unsigned int lev;
    real qsize;

    q = t->root;                                /* start with tree t.root */
    qind = nbSubIndex(p, q);                    /* get index of subcell */
    qsize = t->rsize;                           /* keep track of cell size */
    lev = 0;                                    /* count levels descended */
    while (Subp(q)[qind] != NULL)               /* loop descending tree */
    {
        if (qsize <= REAL_EPSILON)
        {
            if (!t->structureError)
            {
                t->structureError = TRUE;
            }
            return;
        }

        if (isBody(Subp(q)[qind]))              /* reached a "leaf"? */
        {
            c = nbMakeCell(st, t);             /* allocate new cell */
            nbInitMidpoint(c, p, q, qsize);    /* initialize midpoint */

            Subp(c)[nbSubIndex((Body*) Subp(q)[qind], c)] = Subp(q)[qind];
            /* put body in cell */
            Subp(q)[qind] = (NBodyNode*) c;    /* link cell in tree */
        }
        q = (NBodyCell*) Subp(q)[qind];   /* advance to next level */
        qind = nbSubIndex(p, q);          /* get index to examine */
        qsize *= 0.5;                     /* shrink current cell */
        ++lev;                            /* count another level */
    }
    Subp(q)[qind] = (NBodyNode*) p;            /* found place, store p */
    t->maxDepth = MAX(t->maxDepth, lev);  /* remember maximum level */
}

ALWAYS_INLINE
static inline real bmax2Inc(real cmPos, real pPos, real psize)
{
    real dmin;
    dmin = cmPos - (pPos - 0.5 * psize);         /* dist from 1st corner */
    return sqr(mw_fmax(dmin, psize - dmin));      /* sum max distance^2 */
}

ALWAYS_INLINE
static inline real calcSW93MaxDist2(const NBodyCell* p, const mwvector cmpos, real psize)
{
    real bmax2;

    /* compute max distance^2 */
    /* loop over dimensions */
    bmax2 = bmax2Inc(X(cmpos), X(Pos(p)), psize);
    bmax2 += bmax2Inc(Y(cmpos), Y(Pos(p)), psize);
    bmax2 += bmax2Inc(Z(cmpos), Z(Pos(p)), psize);

    return bmax2;
}

/* assign critical radius for cell p, using center-of-mass position
 * cmpos and cell size psize. */
static inline real findRCrit(const NBodyCtx* ctx, const NBodyCell* p, real treeRSize, mwvector cmpos, real psize)
{
    real rc, bmax2;

    if (mw_unlikely(ctx->theta == 0.0))
    {
        /* Do an exact force calculation by always opening cells */
        rc = 2.0 * treeRSize;
        return sqr(rc);
    }

    /* return square of radius */
    switch (ctx->criterion)
    {
        case TreeCode:
            /* use size plus offset */
            rc = psize / ctx->theta + mw_distv(cmpos, Pos(p));
            return sqr(rc);

        case SW93:                           /* use S&W's criterion? */
            /* compute max distance^2 */
            bmax2 = calcSW93MaxDist2(p, cmpos, psize);
            return bmax2 / sqr(ctx->theta);      /* using max dist from cm */

        case BH86:                          /* use old BH criterion? */
            rc = psize / ctx->theta;        /* using size of cell */
            return sqr(rc);

        case InvalidCriterion:
        case Exact: /* Uses separate path */
        default:
            rc = 0.0; /* Stop clang static analysis warning */
            mw_fail("Invalid criterion: %s (%d)\n", showCriterionT(ctx->criterion), ctx->criterion);
	    return 0.0;
    }
}

static inline void nbCheckTreeDim(NBodyTree* tree, real pPos, real cmPos, real halfPsize)
{
    /* CHECKME: Precision: This gets angry as N gets big, and the divisions get small */
    if (   cmPos < pPos - halfPsize       /* if out of bounds */
        || cmPos > pPos + halfPsize)      /* in either direction */
    {
        if (!tree->structureError)
        {
            /* Only print if we don't know about the error
             * already. The error will be caught later and this
             * otherwise could print a lot of noise */
            tree->structureError = TRUE;
        }
    }
}

static inline void nbCheckTreeStructure(NBodyTree* tree, const mwvector pPos, const mwvector cmPos, const real psize)
{
    real halfPsize = 0.5 * psize;

    nbCheckTreeDim(tree, X(pPos), X(cmPos), halfPsize);
    nbCheckTreeDim(tree, Y(pPos), Y(cmPos), halfPsize);
    nbCheckTreeDim(tree, Z(pPos), Z(cmPos), halfPsize);
}


/* hackCofM: descend tree finding center-of-mass coordinates and
 * setting critical cell radii.
 */
static void hackCofM(const NBodyCtx* ctx, NBodyTree* tree, NBodyCell* p, real psize)
{
    int i;
    NBodyNode* q;
    mwvector cmpos = ZERO_VECTOR;                /* init center of mass */

    assert(psize >= REAL_EPSILON);

    Mass(p) = 0.0;                              /* init total mass... */
    for (i = 0; i < NSUB; ++i)                  /* loop over subnodes */
    {
        if ((q = Subp(p)[i]) != NULL)           /* does subnode exist? */
        {
            if (isCell(q))                     /* and is it a cell? */
            {
                hackCofM(ctx, tree, (NBodyCell*) q, 0.5 * psize); /* find subcell cm */
            }

            Mass(p) += Mass(q);                       /* sum total mass */
                                                      /* weight pos by mass */
            mw_incaddv_s(cmpos, Pos(q), Mass(q));     /* sum c-of-m position */
        }
    }

    if (Mass(p) > 0.0)                          /* usually, cell has mass   */
    {
        mw_incdivs(cmpos, Mass(p));            /* so find c-of-m position  */
    }
    else                                        /* but if no mass inside    */
    {
        cmpos = Pos(p);                /* use geo. center for now  */
    }

    nbCheckTreeStructure(tree, Pos(p), cmpos, psize);

    Rcrit2(p) = findRCrit(ctx, p, tree->rsize, cmpos, psize);            /* set critical radius */
    Pos(p) = cmpos;             /* and center-of-mass pos */
}

/* nbMakeTree: initialize tree structure for hierarchical force calculation
 * from body array btab, which contains ctx.nbody bodies.
 */
NBodyStatus nbMakeTree(const NBodyCtx* ctx, NBodyState* st)
{
    Body* p;
    const Body* endp = st->bodytab + st->nbody;
    NBodyTree* t = &st->tree;

    nbNewTree(st, t);                                /* flush existing tree, etc */

    expandBox(t, st->bodytab, st->nbody);            /* and expand cell to fit */
    for (p = st->bodytab; p < endp; p++)             /* loop over bodies... */
    {
        if (Mass(p) != 0.0)                  /* exclude test particles */
            nbLoadBody(st, t, p);              /* and insert into tree */
    }

    /* Check if tree structure error occured */
    if (st->tree.structureError)
        return NBODY_TREE_STRUCTURE_ERROR;
    
    int scatters = 0;
    FILE *fdata = fopen("data.txt", "a");
    FILE *fvel = fopen("velocity.txt","a");
    FILE *fdens = fopen("density.txt","a");
    FILE *fhalo = fopen("halo_prob.txt","a");
    for (int i = 0; i < st->nbody; i++) {
        //if (Mass(&st->bodytab[i]) != 0.0)
	    scatters = scatters + nbPrintNearestNeighborsOptimized(&st->tree, &st->bodytab[i], i, fdata, fvel, fdens,fhalo,ctx);
    }
    
    fclose(fvel);
    fclose(fdata);
    fclose(fdens);
    fclose(fhalo);
    hackCofM(ctx, &st->tree, t->root, t->rsize);   /* find c-of-m coordinates */

    /* Check if tree structure error occured */
    if (st->tree.structureError)
        return NBODY_TREE_STRUCTURE_ERROR;

    threadTree((NBodyNode*) t->root, NULL);        /* add Next and More links */
    if (ctx->useQuad)                           /* including quad moments? */
        hackQuad(t->root);                      /* assign Quad moments */

    return NBODY_SUCCESS;
}
 

#if 0
/* For testing */
static int luaFindRCrit(lua_State* luaSt)
{
    const NBodyCtx* ctx;
    NBodyCell p;  /* Test cell, just need a set position */
    real rSize, pSize;
    mwvector cmPos;

    ctx = checkNBodyCtx(luaSt, 1);
    Pos(&p) = *checkVector(luaSt, 2);
    rSize = luaL_checknumber(luaSt, 3);
    cmPos = *checkVector(luaSt, 4);
    pSize = luaL_checknumber(luaSt, 5);

    lua_pushnumber(luaSt, findRCrit(ctx, &p, rSize, cmPos, pSize));

    return 1;
}

void registerFindRCrit(lua_State* luaSt)
{
    lua_register(luaSt, "findRCrit", luaFindRCrit);
}
#endif
