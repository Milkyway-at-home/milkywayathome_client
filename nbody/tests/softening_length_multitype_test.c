/* This test checks that the softening length machinery works
 * correctly with more than two particle types. 
 *   1. First it reads in a manual bodies file in the new "#type" format.
 *   2. Then it runs nbCacheEps2Indices() to check that each body 
 *      has its own softening length.
 *   3. Finally it writes an output file to compare to an example file
 *      so it can check the new output format.
 *
 * The values in the eps2[] array are junk numbers so the check works.
 */

#include <lua.h>
#include <string.h>
#include <stdlib.h>

#include "nbody_types.h"
#include "nbody.h"
#include "nbody_lua.h"
#include "nbody_manual_bodies.h"
#include "nbody_lua_type_marshal.h"
#include "nbody_io.h"

#define MULTITYPE_NBODY   12
#define MULTITYPE_NTYPES  4

static const char* inputFile    = "manual_bodies_multitype.in";
static const char* outFile      = "softening_length_multitype_test.out";
static const char* expectedFile = "softening_length_multitype_test_expected.out";

/* Expected contents of manual_bodies_multitype.in, in file order */
static const body_t expectedType[MULTITYPE_NBODY] =
    { -2, -2, -2, -1, -1, -1, 1, 1, 1, 2, 2, 2 };
static const real expectedMass[MULTITYPE_NBODY] =
    { 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7 };

/* Read the manual bodies file. */
static Body* readManualBodiesFile(const char* filename, int* nbodyOut)
{
    lua_State* luaSt;
    int nResults;
    Body* bodytab;

    luaSt = nbLuaOpen(FALSE);
    if (!luaSt)
    {
        mw_printf("Failed to open Lua state\n");
        return NULL;
    }

    lua_newtable(luaSt);
    lua_pushstring(luaSt, filename);
    lua_setfield(luaSt, -2, "body_file");

    nResults = nbGenerateManualBodies(luaSt);
    if (nResults != 1)
    {
        mw_printf("generatemanualbodies() did not return a body table\n");
        lua_close(luaSt);
        return NULL;
    }

    lua_remove(luaSt, 1); /* drop the argument table; only the returned body table is needed now */

    bodytab = readModels(luaSt, 1, nbodyOut);
    lua_close(luaSt); /* bodytab is a plain malloc'd array and does not depend on the Lua state */

    return bodytab;
}

int main(void)
{
    int failed = 0;
    int i;
    int nbody = 0;
    Body* bodytab;
    NBodyCtx ctx;
    NBodyState st;
    NBodyFlags nbf;
    real eps2[MULTITYPE_NTYPES * MULTITYPE_NTYPES];
    int eps2_index[MULTITYPE_NTYPES] = { -2, -1, 1, 2 };
    real expectedEps2Min;
    int wrc;

    /* 1. Read the manual bodies file (new #type format, 4 particle types) */
    bodytab = readManualBodiesFile(inputFile, &nbody);
    if (!bodytab || nbody != MULTITYPE_NBODY)
    {
        mw_printf("Expected to read %d bodies from '%s', got %d\n",
                   MULTITYPE_NBODY, inputFile, nbody);
        return 1; /* Nothing else can be checked without the bodies */
    }
    mw_printf("Read %d bodies from '%s'\n", nbody, inputFile);

    for (i = 0; i < MULTITYPE_NBODY; ++i)
    {
        if (Type(&bodytab[i]) != expectedType[i] || idBody(&bodytab[i]) != (unsigned int) i
            || Mass(&bodytab[i]) != expectedMass[i])
        {
            mw_printf("Body %d read incorrectly: expected type %d, id %d, mass %f; "
                       "got type %d, id %u, mass %f\n",
                       i, (int) expectedType[i], i, expectedMass[i],
                       (int) Type(&bodytab[i]), idBody(&bodytab[i]), Mass(&bodytab[i]));
            failed = 1;
        }
    }

    /* 2. Build a minimal ctx/st and run the nbCacheEps2Indices() function. */
    for (i = 0; i < MULTITYPE_NTYPES * MULTITYPE_NTYPES; ++i)
    {
        eps2[i] = (real) (i + 1);
    }

    memset(&ctx, 0, sizeof(ctx));
    ctx.eps2 = eps2;
    ctx.eps2_index = eps2_index;
    ctx.eps2_size = MULTITYPE_NTYPES;
    ctx.SimpleOutput = TRUE; /* keep the output file simple/deterministic for this test */

    memset(&st, 0, sizeof(st));
    st.bodytab = bodytab;
    st.nbody = nbody;

    nbCacheEps2Indices(&ctx, &st);

    expectedEps2Min = eps2[0];
    for (i = 1; i < MULTITYPE_NTYPES * MULTITYPE_NTYPES; ++i)
    {
        if (eps2[i] < expectedEps2Min)
        {
            expectedEps2Min = eps2[i];
        }
    }
    if (ctx.eps2_min != expectedEps2Min)
    {
        mw_printf("ctx.eps2_min incorrect: expected %f, got %f\n",
                   expectedEps2Min, ctx.eps2_min);
        failed = 1;
    }

    for (i = 0; i < MULTITYPE_NBODY; ++i)
    {
        int expectedIdx = -1;
        int j;
        for (j = 0; j < MULTITYPE_NTYPES; ++j)
        {
            if (eps2_index[j] == Type(&bodytab[i]))
            {
                expectedIdx = j;
                break;
            }
        }

        if (expectedIdx < 0)
        {
            mw_printf("Test bug: body %d's type %d has no entry in eps2_index[]\n",
                       i, (int) Type(&bodytab[i]));
            failed = 1;
            continue;
        }

        if (bodytab[i].eps2Index != expectedIdx)
        {
            mw_printf("Body %d (type %d) got eps2Index %d, expected %d\n",
                       i, (int) Type(&bodytab[i]), bodytab[i].eps2Index, expectedIdx);
            failed = 1;
        }
    }

    if (!failed)
    {
        mw_printf("nbCacheEps2Indices() correctly matched all %d particle types\n",
                   MULTITYPE_NTYPES);
    }

    /* 3. Write the bodies back out to check that the output file
     * can be written for a >2-type system. */
    memset(&nbf, 0, sizeof(nbf));
    nbf.outFileName = outFile;
    /* Deliberately nonexistent: this test never runs a real workunit
     * script, and nbWriteBodies()'s internal histogram-parameter lookup
     * needs a path it can safely fail to open (as opposed to NULL, which
     * isn't handled the same way) rather than crash. */
    nbf.inputFile = "no_such_lua_script_for_this_test.lua";
    nbf.outputBinary = FALSE;

    wrc = nbWriteBodies(&ctx, &st, &nbf);
    if (wrc != 0)
    {
        mw_printf("nbWriteBodies() failed to write '%s'\n", outFile);
        failed = 1;
    }
    else
    {
        mw_printf("Wrote output bodies to '%s'\n", outFile);

        /* Compare against the expected output file.*/
        FILE* ef = fopen(expectedFile, "r");
        if (!ef)
        {
            mw_printf("Expected output file '%s' not found -- cannot verify '%s'\n",
                       expectedFile, outFile);
            failed = 1;
        }
        else
        {
            FILE* af = fopen(outFile, "r");
            int mismatch = 0;
            int ec, ac;

            if (!af)
            {
                mw_printf("Could not reopen '%s' for comparison\n", outFile);
                mismatch = 1;
            }
            else
            {
                do
                {
                    ec = fgetc(ef);
                    ac = fgetc(af);
                    if (ec != ac)
                    {
                        mismatch = 1;
                        break;
                    }
                } while (ec != EOF);
                fclose(af);
            }
            fclose(ef);

            if (mismatch)
            {
                mw_printf("Output file '%s' does not match expected file '%s'\n",
                           outFile, expectedFile);
                failed = 1;
            }
            else
            {
                mw_printf("Output file matches expected file '%s'\n", expectedFile);
            }
        }
    }

    free(bodytab);
    return failed;
}
