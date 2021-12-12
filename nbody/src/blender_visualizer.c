#include "blender_visualizer.h"
#include "milkyway_math.h"

static int appendFileNum(char* prefix, const unsigned int filenumber)
{
    char zeros[5];
    if(!prefix)
    {
        return 0;
    }
    if(filenumber < 10)
    {
        sprintf(zeros, "0000");
    }
    else if(filenumber < 100)
    {
        sprintf(zeros, "000");
    }
    else if(filenumber < 1000)
    {
        sprintf(zeros, "00");
    }
    else if(filenumber < 10000)
    {
        sprintf(zeros, "0");
    }
    else
    {
        zeros[0] = '\0';
    }
    sprintf(prefix, "%s%s%u", prefix, zeros, filenumber);
    return 1;
}


/* Use the center of mass if we have it already in some form,
 * otherwise we can find it */
int nbFindCenterOfMass(mwvector* cmPos, const NBodyState* st) /* While the printCOM stuff is no longer being used, this is still used for the camera. Do not remove. */
{
    if (st->tree.root)
    {
        *cmPos = Pos(st->tree.root);
    }
    else if (st->usesCL)
    {
      #if NBODY_OPENCL
        if (nbDisplayUpdateMarshalBodies(st, cmPos))
        {
            return 1;
        }
      #endif
    }
    else
    {
        /* If we are using exact nbody or haven't constructed the tree
         * yet we need to calculate the center of mass on our own */
        *cmPos = nbCenterOfMass(st);
    }

    return 0;
}

NBodyStatus deleteOldFiles(const NBodyState* st)
{
    int nbody = st->nbody;
    if (nbody != 1) /* Particle sim */
    {
        FILE *f = fopen("blender_event_record.txt", "w");
        if(!f)
        {
            return NBODY_ERROR;
        }
        fclose(f);
        f = fopen("blender_misc_record.txt", "w");
        if(!f)
        {
            return NBODY_ERROR;
        }
        fclose(f);
    }
    else            /* Orbit sim */
    {
        FILE *f = fopen("blender_orbit_record.txt", "w");
        if(!f)
        {
            return NBODY_ERROR;
        }
        fclose(f);
    }
    return NBODY_SUCCESS;
}

/* Particle positions info */
NBodyStatus blenderPrintBodies(const NBodyState* st, const NBodyCtx* ctx)
{
    const Body* b;
    int nbody = st->nbody;
    FILE *f;
    
    if (st->step > 99999)
    {
        printf("Tried to write more than 99999 frames.\n");
        return NBODY_ERROR;
    }
    if (nbody != 1) 
    {
        char filename[19] = "frames/frame_";
        appendFileNum(filename, st->step);
        f = fopen(filename, "w"); /* We are simulating the particles, put this in the right file */
        if(!f)
        {
            return NBODY_ERROR;
        }
    }
    else
    {
        f = fopen("blender_orbit_record.txt", "a"); /* We are simulating the orbit instead, which is a one-particle simulation */
        if(!f)
        {
            return NBODY_ERROR;
        }
    }
    for (int i = 0; i < nbody; i++)
    {
        b = &st->bodytab[i];
        float x = (float) showRealValue(X(Pos(b)));
        float y = (float) showRealValue(Y(Pos(b)));
        float z = (float) showRealValue(Z(Pos(b)));
        fprintf(f, "%0.4f %0.4f %0.4f\n",x,y,z);
    }
    fclose(f);
    return NBODY_SUCCESS;
}

/* Center of mass position info */
NBodyStatus blenderPrintCOM(const NBodyState* st)
{
    mwvector cmPos;
    scene_t* scene = st->scene;

    if (!scene)
    {
        return NBODY_SUCCESS;
    }

    if (nbFindCenterOfMass(&cmPos, st)) 
    {
        return NBODY_ERROR;
    }
    FILE *f = fopen("blender_event_record.txt", "a");
    if(!f)
    {
        return NBODY_ERROR;
    }
    float x = (float) showRealValue(X(cmPos));
    float y = (float) showRealValue(Y(cmPos));
    float z = (float) showRealValue(Z(cmPos));
    fprintf(f, "%0.4f %0.4f %0.4f\n", x,y,z);
    fclose(f);
    return NBODY_SUCCESS;
}

/* Dark vs light matter, total time, and vectors on tidal plane (for camera to look head-on at the plane) */
NBodyStatus blenderPrintMisc(const NBodyState* st, const NBodyCtx* ctx, mwvector cmPos1, mwvector cmPos2)
{
    /* This is the orbit sim (one particle), so the particle sim will have the misc data covered. */
    if (st->nbody == 1) 
    {
        return NBODY_SUCCESS;
    }

    FILE *f = fopen("blender_misc_record.txt", "a");
    if(!f)
    {
        return NBODY_ERROR;
    }

    /* Number of particles, number of frames, evolve time */
    fprintf(f,"%d\n%d\n%f\n", st->nbody, st->step, ctx->timeEvolve);

    /* To find a plane normal to the orbit for the camera to point at */
    fprintf(f,"%f\n%f\n%f\n", showRealValue(X(cmPos1)), showRealValue(Y(cmPos1)), showRealValue(Z(cmPos1)));
    fprintf(f,"%f\n%f\n%f\n", showRealValue(X(cmPos2)), showRealValue(Y(cmPos2)), showRealValue(Z(cmPos2)));

    /* Light or dark particles */
    const Body* b;
    int nbody = st->nbody;
    for (int i = 0; i < nbody; i++)
    {
        b = &st->bodytab[i];
        fprintf(f, "%d\n",ignoreBody(b)); /* 0 if light, 1 if dark */
    }
    fclose(f);
    printf("*Total frames: %u\n", st->step);
    return NBODY_SUCCESS;
}

/*We are looking for three points (two vectors) to approximate the plane of the orbit. One point is (0,0,0), and another is the COM's start position.
* A vector from (0,0,0) to the third point should be close to perpedicular to a vector from (0,0,0) to the COM's start position.
* We find the best third point to be that which is most perpedicular to that vector. */
NBodyStatus blenderPossiblyChangePerpendicularCmPos(mwvector* next, mwvector* perp, mwvector* start)
{

    /* First get a normalizer for the potential third points (no need to normalize "start", since the division would fall out at the if statement) */
    float oldNormer = showRealValue(mw_hypot(mw_hypot(perp->x, perp->y), perp->z));
    float newNormer = showRealValue(mw_hypot(mw_hypot(next->x, next->y), next->z));

    /* Dot product goes to zero as vectors become more perpendicular. Do not favor third points that happen to be closer to zero with "normer". */
    float oldDotAbs = mw_fabs_0(showRealValue(mw_dotv(*perp, *start)))/oldNormer;
    float newDotAbs = mw_fabs_0(showRealValue(mw_dotv(*next, *start)))/oldNormer;

    /* If the new dot product is closer to zero than the old one, drop our old third point and keep this better one for now. */
    if (newDotAbs < oldDotAbs)
    {
        perp->x=next->x;
        perp->y=next->y;
        perp->z=next->z;
    }
    return NBODY_SUCCESS;
}
