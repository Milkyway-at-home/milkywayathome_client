/* ************************************************************************** */
/* GETPARAM.C: export version prompts user for values. Public routines: */
/* initparam(), getparam(), getiparam(), getbparam(), getrparam(). */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "stdinc.h"
#include "real.h"
#include <string.h>
#include <stdlib.h>

void* allocate(int);
void error(char*, ...);

static char** defaults = NULL;          /* "name=value" char*s */

/*  * INITPARAM: ignore arg vector, remember defaults.
 */

void initparam(char** argv, char** defv)
{
    defaults = defv;
}

/*  * GETPARAM: export version prompts user for value of name.
 */

static int scanbind(char**, char*);
static char* extrvalue(char*);

char* getparam(char* name)
{
    int i, len;
    char* def;
    char buf[128];

    if (defaults == NULL)           /* check initialization */
        error("getparam: called before initparam\n");
    i = scanbind(defaults, name);       /* find name in defaults */
    if (i < 0)
        error("getparam: %s unknown\n", name);
    def = extrvalue(defaults[i]);       /* extract default value */
    if (*def == NULL)
        fprintf(stderr, "enter %s: ", name);    /* prompt user for value */
    else
        fprintf(stderr, "enter %s [%s]: ", name, def);
    gets(buf);                  /* read users response */
    len = strlen(buf);
    if (len > 0)                /* if user gave a value... */
        return (strcpy((char*) allocate(len + 1), buf));
    else                    /* else return default */
        return (def);
}

/*  * GETIPARAM, ..., GETDPARAM: get int, long, bool, or double parameters.
 */

int getiparam(char* name)
{
    char* val;

    for (val = ""; *val == NULL; )      /* while nothing input */
        val = getparam(name);                   /* obtain value from user */
    return (atoi(val));                         /* convert to an integer */
}

bool getbparam(char* name)
{
    char* val;

    for (val = ""; *val == NULL; )
        val = getparam(name);
    if (strchr("tTyY1", *val) != NULL)      /* is value true? */
        return (TRUE);
    if (strchr("fFnN0", *val) != NULL)      /* is value false? */
        return (FALSE);
    error("getbparam: %s=%s not bool\n", name, val);
    return (FALSE);
}

real getrparam(char* name)
{
    char* val;

    for (val = ""; *val == NULL; )
        val = getparam(name);
    return ((real) atof(val));          /* convert & return real */
}

/*  * SCANBIND: scan binding vector for name, return index.
 */

static bool matchname(char*, char*);

static int scanbind(char* bvec[], char* name)
{
    int i;

    for (i = 0; bvec[i] != NULL; i++)
        if (matchname(bvec[i], name))
            return (i);
    return (-1);
}

/*  * MATCHNAME: determine if "name=value" matches "name".
 */

static bool matchname(char* bind, char* name)
{
    char* bp, *np;

    bp = bind;
    np = name;
    while (*bp == *np)
    {
        bp++;
        np++;
    }
    return (*bp == '=' && *np == NULL);
}

/*  * EXTRVALUE: extract value from name=value char*.
 */

static char* extrvalue(char* arg)
{
    char* ap;

    ap = (char*) arg;
    while (*ap != NULL)
        if (*ap++ == '=')
            return ((char*) ap);
    return (NULL);
}
