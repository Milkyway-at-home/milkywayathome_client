/* This test checks that the softening length (eps2) machinery works
 * correctly with more than two particle types. It does not run a
 * simulation -- it only exercises three pieces of production code against
 * a manual bodies file containing four particle types (two light, two
 * dark):
 *
 *   1. nbGenerateManualBodies()/readModels() -- reading a manual bodies
 *      file in the new "#type" format into a flat Body array.
 *   2. nbCacheEps2Indices() -- matching each body's type to its row/column
 *      in ctx->eps2[], which is the function that would fail if the
 *      eps2/eps2_index machinery didn't generalize past two types.
 *   3. nbWriteBodies() -- writing the resulting bodies back out, so the
 *      output path is also exercised for a >2-type system.
 *
 * The values in the eps2[] array itself are irrelevant here (they are
 * never used by a real force calculation in this test) -- only the
 * eps2_index[] -> body type matching matters.
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

/* Read the manual bodies file through the same Lua-facing entry point
 * (generatemanualbodies) and Lua->C body marshalling (readModels()) that
 * a real run uses -- just invoked directly instead of from a workunit
 * script, since no simulation is being run. */
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

    /* Build the { body_file = filename } named-argument table that
     * nbGenerateManualBodies() (registered as the Lua global
     * generatemanualbodies()) expects as its one argument. */
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

    /* 2. Build a minimal ctx/st and run the function under test:
     * nbCacheEps2Indices() matches each body's type to its eps2_index[]
     * row/column. The eps2[] values themselves are never read by a force
     * calculation in this test, so their exact values don't matter -- only
     * eps2_index[] (the set of types with an entry) matters. */
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

    /* 3. Write the bodies back out. No simulation is run -- this just
     * checks that an output file can be written for a >2-type system. */
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

        /* If a golden output file has been checked in alongside this test,
         * compare against it byte-for-byte. Until then, there's nothing to
         * compare against -- see the project doc for how to generate one. */
        FILE* ef = fopen(expectedFile, "r");
        if (!ef)
        {
            mw_printf("No golden output fixture at '%s' yet -- skipping byte-for-byte "
                       "comparison. After checking that '%s' looks correct, copy it to "
                       "'%s' to enable the comparison on future runs.\n",
                       expectedFile, outFile, expectedFile);
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
                mw_printf("Output file '%s' does not match golden fixture '%s'\n",
                           outFile, expectedFile);
                failed = 1;
            }
            else
            {
                mw_printf("Output file matches golden fixture '%s'\n", expectedFile);
            }
        }
    }

    free(bodytab);
    return failed;
}
