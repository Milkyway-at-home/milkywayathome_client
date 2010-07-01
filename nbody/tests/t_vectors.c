/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#include "nbody_tests.h"
#include "nbody_priv.h"

/* When the result is supposed to be a vector */
static void vectorError(const char* testfunc, const char* op, vector reference, vector answer, vector a, vector b)
{
    char* buf1;
    char* buf2;
    char* buf3;
    char* buf4;
    char* buf5;

    vector diff;

    diff[0] = reference[0] - answer[0];
    diff[1] = reference[1] - answer[1];
    diff[2] = reference[2] - answer[2];

    buf1 = showVector(a);
    buf2 = showVector(b);
    buf3 = showVector(answer);
    buf4 = showVector(reference);
    buf5 = showVector(diff);

    fprintf(stderr,
            "Failed vector test %s: \n"
            "\t%s %s %s,\n"
            "\tgot %s,\n"
            "\texpected %s.\n"
            "\tDifference = %s\n",
            testfunc,
            buf1,
            op,
            buf2,
            buf3,
            buf4,
            buf5);

    free(buf1);
    free(buf2);
    free(buf3);
    free(buf4);
    free(buf5);
}


static void vectorIncError(const char* testfunc, const char* op, vector reference, vector answer, vector aOrig, vector b)
{
    char* buf1;
    char* buf2;
    char* buf3;
    char* buf4;
    char* buf5;

    vector diff;

    diff[0] = reference[0] - answer[0];
    diff[1] = reference[1] - answer[1];
    diff[2] = reference[2] - answer[2];

    buf1 = showVector(aOrig);
    buf2 = showVector(b);
    buf3 = showVector(answer);
    buf4 = showVector(reference);
    buf5 = showVector(diff);

    fprintf(stderr,
            "Failed vectorinc test %s: \n"
            "\t%s(%s, %s),\n"
            "\tgot %s,\n"
            "\texpected %s.\n"
            "\tDifference = %s\n",
            testfunc,
            op,
            buf1,
            buf2,
            buf3,
            buf4,
            buf5);

    free(buf1);
    free(buf2);
    free(buf3);
    free(buf4);
    free(buf5);
}

static void vectorErrorOne(const char* testfunc, const char* op, vector reference, vector answer, vector a)
{
    char* buf1;
    char* buf2;
    char* buf3;
    char* buf4;

    vector diff;

    diff[0] = reference[0] - answer[0];
    diff[1] = reference[1] - answer[1];
    diff[2] = reference[2] - answer[2];

    buf1 = showVector(a);
    buf2 = showVector(answer);
    buf3 = showVector(reference);
    buf4 = showVector(diff);

    fprintf(stderr,
            "Failed vector test %s: %s(%s),\n"
            "\tgot %s, expected %s.\n\tDifference = %s\n",
            testfunc,
            op,
            buf1,
            buf2,
            buf3,
            buf4);

    free(buf1);
    free(buf2);
    free(buf3);
    free(buf4);
}

/* When the result is supposed to be a scalar */
static void scalarError(const char* testfunc, const char* op, real reference, real answer, vector a, vector b)
{
    char* buf1;
    char* buf2;

    buf1 = showVector(a);
    buf2 = showVector(b);

    fprintf(stderr,
            "Failed scalar test %s : ( %s %s %s), got %g, "
            "expected %g. Difference = %g, error = %g%%\n",
            testfunc,
            buf1,
            op,
            buf2,
            answer,
            reference,
            reference - answer,
            REL_ERROR(reference, answer));

    free(buf1);
    free(buf2);
}

static void scalarIncError(const char* testfunc, const char* op, vector reference, vector answer, vector aOrig, real b)
{
    char* buf1;
    char* buf2;
    char* buf3;
    char* buf4;

    vector diff;

    diff[0] = reference[0] - answer[0];
    diff[1] = reference[1] - answer[1];
    diff[2] = reference[2] - answer[2];

    buf1 = showVector(aOrig);
    buf2 = showVector(answer);
    buf3 = showVector(reference);
    buf4 = showVector(diff);

    fprintf(stderr,
            "Failed scalar inc test %s : %s(%s, %g),\n"
            "\tgot %s,\n"
            "\texpected %s.\n"
            "\tDifference = %s\n",
            testfunc,
            op,
            buf1,
            b,

            buf2,
            buf3,
            buf4);

    free(buf1);
    free(buf2);
    free(buf3);
    free(buf4);
}


/* for :: Vector -> Scalar -> Vector ones */
static void vectorScalarError(const char* testfunc, const char* op, vector reference, vector answer, vector a, real b)
{
    char* buf1;
    char* buf2;
    char* buf3;
    char* buf4;
    vector diff;

    diff[0] = reference[0] - answer[0];
    diff[1] = reference[1] - answer[1];
    diff[2] = reference[2] - answer[2];

    buf1 = showVector(a);
    buf2 = showVector(answer);
    buf3 = showVector(reference);
    buf4 = showVector(diff);

    fprintf(stderr,
            "Failed vector-scalar test %s:\n"
            "\t%s(%s, %g),\n"
            "\tgot %s,\n"
            "\texpected %s.\n"
            "\tDifference = %s\n",
            testfunc,
            op,
            buf1,
            b,
            buf2,
            buf3,
            buf4);

    free(buf1);
    free(buf2);
    free(buf3);
    free(buf4);
}


static void scalarErrorOne(const char* testfunc, const char* op, real reference, real answer, vector a)
{
    char* buf1;

    buf1 = showVector(a);

    fprintf(stderr,
            "Failed test %s : %s(%s), got %g, "
            "expected %g. Difference = %g, error = %g%%\n",
            testfunc,
            op,
            buf1,
            reference,
            answer,
            reference - answer,
            REL_ERROR(reference, answer));

    free(buf1);
}

/*****************************************************************************/

/* Test the dot product function against a simple reference */
inline static int test_dotvp()
{
    int failed = 0;

    vector a = RANDOM_VECTOR;
    vector b = RANDOM_VECTOR;

    real answer;

    real reference = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);

    DOTVP(answer, a, b);

    if (answer != reference)
    {
        failed = 1;
        scalarError(__func__, ".", reference, answer, a, b);
    }

    return failed;
}

inline static int test_sqrv()
{
    int failed = 0;
    vector a = RANDOM_VECTOR;
    real answer;
    real reference = (a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]);

    SQRV(answer, a);

    if (answer != reference)
    {
        failed = 1;
        scalarErrorOne(__func__, "`sqrv`", reference, answer, a);
    }

    return failed;
}

inline static int test_subv()
{
    int failed = 0;
    vector answer;
    vector a = RANDOM_VECTOR;
    vector b = RANDOM_VECTOR;
    vector reference = { a[0] - b[0],  a[1] - b[1], a[2] - b[2] };

    SUBV(answer, a, b);

    if (!VECEQ(answer, reference))
    {
        failed = 1;
        vectorError(__func__, "-", reference, answer, a, b);
    }

    return failed;
}

inline static int test_mulvs()
{
    int failed = 0;
    vector answer;
    vector a = RANDOM_VECTOR;
    real s = RANDOM_REAL;
    vector reference = { s * a[0],  s * a[1], s * a[2] };

    MULVS(answer, a, s);

    if (!VECEQ(answer, reference))
    {
        failed = 1;
        vectorScalarError(__func__, "mulvs", reference, answer, a, s);
    }

    return failed;
}

inline static int test_divvs()
{
    int failed = 0;
    vector answer;
    vector a = RANDOM_VECTOR;
    real s = RANDOM_REAL;
    vector reference = { a[0] / s, a[1] / s, a[2] /s };

    DIVVS(answer, a, s);

    if (!VECEQ(answer, reference))
    {
        failed = 1;
        vectorScalarError(__func__, "divvs", reference, answer, a, s);
    }

    return failed;
}

inline static int test_addvs()
{
    int failed = 0;
    vector answer;
    vector a = RANDOM_VECTOR;
    real s = RANDOM_REAL;
    vector reference = { s + a[0],  s + a[1], s + a[2] };

    ADDVS(answer, a, s);

    if (!VECEQ(answer, reference))
    {
        failed = 1;
        vectorScalarError(__func__, "addvs", reference, answer, a, s);
    }

    return failed;
}


inline static int test_setv()
{
    int failed = 0;
    vector answer;
    vector a = RANDOM_VECTOR;
    vector b = RANDOM_VECTOR;
    vector reference;

    vector aOrig;

    COPYVECTOR(aOrig, a);

    COPYVECTOR(reference, b);

    SETV(answer, b);

    if (!VECEQ(answer, reference))
    {
        failed = 1;
        vectorError(__func__, "`setv`", reference, answer, aOrig, b);
    }

    return failed;
}

inline static int test_addv()
{
    int failed = 0;
    vector answer;
    vector a = RANDOM_VECTOR;
    vector b = RANDOM_VECTOR;
    vector reference = { a[0] + b[0],  a[1] + b[1], a[2] + b[2] };

    ADDV(answer, a, b);

    if (!VECEQ(answer, reference))
    {
        failed = 1;
        vectorError(__func__, "+", reference, answer, a, b);
    }

    return failed;
}

inline static int test_clrv()
{
    int failed = 0;
    vector a = RANDOM_VECTOR;
    vector reference = ZERO_VECTOR;
    vector aOrig;

    COPYVECTOR(aOrig, a);

    CLRV(a);

    if (!VECEQ(a, reference))
    {
        failed = 1;
        vectorErrorOne(__func__, "clear", reference, a, aOrig);
    }

    return failed;
}

inline static int test_incaddv()
{
    int failed = 0;
    vector a = RANDOM_VECTOR;
    vector b = RANDOM_VECTOR;
    vector reference = { a[0] + b[0],  a[1] + b[1], a[2] + b[2] };
    vector aOrig;

    COPYVECTOR(aOrig, a);

    INCADDV(a, b);

    if (!VECEQ(a, reference))
    {
        failed = 1;
        vectorIncError(__func__, "incaddv", reference, a, aOrig, b);
    }

    return failed;
}


inline static int test_incmulvs()
{
    int failed = 0;
    vector a = RANDOM_VECTOR;
    real s = RANDOM_REAL;
    vector reference = { s * a[0], s * a[1], s * a[2] };
    vector aOrig;

    COPYVECTOR(aOrig, a);

    INCMULVS(a, s);

    if (!VECEQ(a, reference))
    {
        failed = 1;
        scalarIncError(__func__, "incmulvs", reference, a, aOrig, s);
    }

    return failed;
}

inline static int test_incdivvs()
{
    int failed = 0;
    vector a = RANDOM_VECTOR;
    real s = RANDOM_REAL;
    vector reference = { a[0] / s, a[1] / s, a[2] / s };
    vector aOrig;

    COPYVECTOR(aOrig, a);

    INCDIVVS(a, s);

    if (!VECEQ(a, reference))
    {
        failed = 1;
        scalarIncError(__func__, "incdivvs", reference, a, aOrig, s);
    }

    return failed;
}

inline static int test_incsubv()
{
    int failed = 0;
    vector a = RANDOM_VECTOR;
    vector b = RANDOM_VECTOR;
    vector reference = { a[0] - b[0],  a[1] - b[1], a[2] - b[2] };
    vector aOrig;

    COPYVECTOR(aOrig, a);

    INCSUBV(a, b);

    if (!VECEQ(a, reference))
    {
        failed = 1;
        vectorIncError(__func__, "incsubv", reference, a, aOrig, b);
    }

    return failed;
}

inline static int test_absv()
{
    int failed = 0;
    real answer;
    vector a = RANDOM_VECTOR;
    real reference = sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]));

    ABSV(answer, a);

    if (answer != reference)
    {
        failed = 1;
        scalarErrorOne(__func__, "abs", reference, answer, a);
    }

    return failed;
}


inline static int test_unitv()
{
    vector a = RANDOM_VECTOR;
    real absval;
    unsigned int dir = RANDOM_INT(0, 2);
    char* buf;

    UNITV(a, dir);
    ABSV(absval, a);

    if ( absval != 1.0 || a[dir] != 1.0)
    {
        buf = showVector(a);
        fprintf(stderr,
                "Failed test_unitv: got %s, %u, abs /= 1.0\n",
                buf,
                dir);
        free(buf);
        return 1;
    }

    return 0;
}

inline static int test_incnegv()
{
    char* origStr;
    char* aStr;
    char* refStr;

    vector reference, aOrig;
    vector a = RANDOM_VECTOR;


    COPYVECTOR(aOrig, a);
    COPYVECTOR(reference, a);

    reference[0] = -reference[0];
    reference[1] = -reference[1];
    reference[2] = -reference[2];

    INCNEGV(a);


    if (!VECEQ(a, reference))
    {
        aStr    = showVector(a);
        origStr = showVector(aOrig);
        refStr  = showVector(reference);
        fprintf(stderr,
                "Failed test_incnegv:\n"
                "\tvector = %s\n"
                "\tgot %s, expected %s\n",
                origStr,
                aStr,
                refStr);
        free(aStr);
        free(origStr);
        free(refStr);
        return 1;
    }

    return 0;
}

inline static int test_negv()
{
    char* origStr;
    char* aStr;
    char* refStr;

    vector answer;
    vector reference;
    vector a = RANDOM_VECTOR;


    reference[0] = -a[0];
    reference[1] = -a[1];
    reference[2] = -a[2];

    NEGV(answer, a);

    if (!VECEQ(answer, reference))
    {
        aStr    = showVector(answer);
        origStr = showVector(a);
        refStr  = showVector(reference);
        fprintf(stderr,
                "Failed test_negvv:\n"
                "\tvector = %s\n"
                "\tgot %s, expected %s\n",
                origStr,
                aStr,
                refStr);
        free(aStr);
        free(origStr);
        free(refStr);
        return 1;
    }

    return 0;
}


/* One argument as number of tests to run */
int main(int argc, char** argv)
{
    size_t i;
    unsigned int fails = 0;
    unsigned int numTests;

    if (argc != 2)
        numTests = 1000000;
    else
        numTests = strtol(argv[1], NULL, 10);

    for ( i = 0; i < numTests; ++i)
    {
        fails += test_clrv();
        fails += test_unitv();
        fails += test_setv();
        fails += test_addv();
        fails += test_subv();
        fails += test_mulvs();
        fails += test_divvs();
        fails += test_dotvp();
        fails += test_sqrv();
        fails += test_absv();

        fails += test_incaddv();
        fails += test_incsubv();
        fails += test_incmulvs();
        fails += test_incdivvs();
        fails += test_addvs();

        fails += test_negv();
        fails += test_incnegv();


        /*
        fails += test_distv();
        fails += test_crossvp();
        */

        /* matrix ops */
        /*
        fails += test_clrm();
        fails += test_setmi();
        fails += test_setm();
        fails += test_tranm();
        fails += test_subm();
        fails += test_mulm();
        fails += test_mulms();
        fails += test_divms();
        fails += test_mulmv();
        fails += test_outvp();
        fails += test_tracem();

        fails += test_setms();
        */

    }

    if (fails)
        fprintf(stderr, "Vector tests: %u out of %u tests failed\n", fails, 17 * numTests);

    return fails;
}


