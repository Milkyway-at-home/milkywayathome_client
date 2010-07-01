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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>

static int runUninterrupted(const char* binPath, const char* config, const char* outfile, const char* checkpoint)
{
    pid_t pid;
    int status;

    const char* args[] = { binPath, "-f", config, "-o", outfile, "-c", checkpoint, NULL };

    pid = fork();

    if (pid == 0) /* child */
    {
        if (execv(binPath, args) < 0)
        {
            perror("execv");
            exit(EXIT_FAILURE);
        }
    }
    else if (pid < 0)
    {
        fprintf(stderr, "Failed to fork\n");
        exit(EXIT_FAILURE);
    }
    else  /* parent */
    {
        if (waitpid(pid, &status, 0) < 0)
        {
            perror("waitpid");
            exit(EXIT_FAILURE);
        }
    }
    return 0;
}

static int runResume(const char* binPath, const char* config, const char* outfile, const char* checkpoint)
{
    pid_t pid;
    int status;

    const char* args[] = { binPath, "-f", config, "-o", outfile, "-c", checkpoint, NULL };

    pid = fork();

    if (pid == 0) /* child */
    {
        if (execv(binPath, args) < 0)
        {
            perror("execv");
            exit(EXIT_FAILURE);
        }
    }
    else if (pid < 0)
    {
        fprintf(stderr, "Failed to fork\n");
        exit(EXIT_FAILURE);
    }
    else  /* parent */
    {
        if (waitpid(pid, &status, 0) < 0)
        {
            perror("waitpid");
            exit(EXIT_FAILURE);
        }
    }

    return 0;
}

static int runInterrupt(const char* binPath, const char* config, const char* outfile, const char* checkpoint, unsigned int stime)
{
    pid_t pid;

    const char* args[] = { binPath, "-f", config, "-o", outfile, "-c", checkpoint, NULL };

    pid = fork();

    if (pid == 0) /* child */
    {
        if (execv(binPath, args) < 0)
        {
            perror("execv");
            exit(EXIT_FAILURE);
        }
    }
    else if (pid < 0)
    {
        fprintf(stderr, "Failed to fork\n");
        exit(EXIT_FAILURE);
    }
    else  /* parent */
    {
       /* Let it run a bit, then kill it some random number of times
         * at random intervals. Then we restart it. The result should
         * be the same as the uninterrupted one. Also assumes the
         * boinc checkpoint interval is set lower */

        printf("Sleeping for %u\n", stime);
        sleep(stime);
        printf("Awake\n");

        if ( kill(pid, SIGQUIT) < 0 )
        {
            perror("kill");
            exit(EXIT_FAILURE);
        }
    }

    return 0;
}

int main(int argc, char** argv)
{
    int rc;
    unsigned int n, numInt = 2;

    const char* binPath       = "bin/milkyway_nbody";
    const char* check         = "test_check";
    const char* testFile      = "medium_test.js";
    const char* testOutNormal = "test_out_normal";
    const char* testOutInt    = "test_out_interrupted";

    remove(check);

    printf("Running normally\n");
    runUninterrupted(binPath, testFile, testOutNormal, check);
    printf("Uninterrupted run complete\n");

    for ( n = 0; n < numInt; ++n )
    {
        printf("Beginning interrupted run %u\n", n);
        runInterrupt(binPath, testFile, testOutInt, check, 10);
        printf("Interrupted run %u complete\n", n);
    }

    printf("Begin resume run\n");
    runResume(binPath, testFile, testOutInt, check);
    printf("Resumed run complete\n");

    /* Check if result files differ */
    if ((rc = system("diff test_out_normal test_out_interrupted")) < 0)
    {
        perror("system");
        exit(EXIT_FAILURE);
    }

    return (rc != 0); /* Failure if they differ */
}


