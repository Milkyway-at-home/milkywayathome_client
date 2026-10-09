/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
     Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright (c) 2016-2018 Siddhartha Shelton
This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.    If not, see <http://www.gnu.org/licenses/>.

The minimum bracketing method is based on results from "Numerical Recipes,
3rd ed." and conforms to the authors' defintion of intellectual
property under their license, and as such does not conflict with
their copyright to their programs which execute similar algorithms.
*/
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "milkyway_math.h"
#include "milkyway_lua.h"
#include "nbody_lua_types.h"
#include "nbody_manual_bodies.h"
#include "nbody_types.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*      DWARF GENERATION        */
static int nbGenerateManualBodiescore(lua_State* luaSt, const char* body_file)
{   
    
    /*initializing particles:*/
    int table;
    Body b;
    FILE* body_inputs;

    unsigned int lineNum = 0;
    char lineBuf[1024];
    int rc __attribute__((unused)) = 0;
    mwbool error __attribute__((unused)) = FALSE;
    
    
    body_inputs = mwOpenResolved(body_file, "r");
    size_t fsize = mwCountLinesInFile(body_inputs);//get the number of lines in the file
    
    
    if (body_inputs == NULL)//make sure the file is available
    {
        mw_printf("Error opening data file '%s'\n", body_file);
        return 0;
    }
    
    if (fsize == 0)//if the file is empty then throw error
    {
        mw_printf("Data file line count = 0\n");
        return 0;
    }
    
    char headerLine[sizeof(lineBuf)] = "";
    mwbool sawFirstLine = FALSE;

    while (fgets(lineBuf, (int) sizeof(lineBuf), body_inputs))
    {
        if (!sawFirstLine)
        {
            sawFirstLine = TRUE;
            if (lineBuf[0] == '#')
            {
                strncpy(headerLine, lineBuf, sizeof(headerLine) - 1);
            }
        }

        /* Skip comments and blank lines */
        if (lineBuf[0] == '#' || lineBuf[0] == '\n')
        {
            fsize -= 1;//removing the commented and blank lines from the total line count
        }
        continue;
    }
    fclose(body_inputs);
    body_inputs = mwOpenResolved(body_file, "r");

    /* The column header comment line tells us which convention the file's
     * first column uses: a "#type" header is the new format, where the
     * column is the body's signed type directly (negative = dark matter,
     * positive = light matter). Anything else (like "#ignore") is the legacy
     * format, where the column is a boolean (0 = baryon, 1 = dark matter). */
    mwbool newBodyTypeFormat = (strstr(headerLine, "type") != NULL);
   
    real * ty = mwCalloc(fsize, sizeof(real));
    real * id = mwCalloc(fsize, sizeof(real));
    real * x  = mwCalloc(fsize, sizeof(real));
    real * y  = mwCalloc(fsize, sizeof(real));
    real * z  = mwCalloc(fsize, sizeof(real));
    real * vx = mwCalloc(fsize, sizeof(real));
    real * vy = mwCalloc(fsize, sizeof(real));
    real * vz = mwCalloc(fsize, sizeof(real));
    real * masses = mwCalloc(fsize, sizeof(real));
    
    unsigned int nbody = fsize;
    memset(&b, 0, sizeof(b));
    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);      
    

    int counter = 0;
    while (fgets(lineBuf, (int) sizeof(lineBuf), body_inputs))
    {
        ++lineNum;

        if (strlen(lineBuf) + 1 >= sizeof(lineBuf))//check if we can store the data from that line
        {
            mw_printf("Error reading data file line %d (Line buffer too small): %s", lineNum, lineBuf);
            error = TRUE;
            break;
        }

        /* Skip comments and blank lines */
        if (lineBuf[0] == '#' || lineBuf[0] == '\n')
            continue;

        rc = sscanf(lineBuf,
                    "%lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
                    &ty[counter],
                    &id[counter],
                    &x[counter],
                    &y[counter],
                    &z[counter],
                    &vx[counter],
                    &vy[counter],
                    &vz[counter],
                    &masses[counter]
                   );
        
        counter++;
    }
    
    fclose(body_inputs);

    /* pushing the bodies */
    for (int i = 0; i < (int) nbody; i++)
    {
        if (newBodyTypeFormat)
        {
            /* New format: the column is the body's particle type.
             * Light matter is positive and dark matter is negative.
             * Bodies must never be 0 since the tree code uses type 0.*/
            if (ty[i] == 0.0)
            {
                mw_fail("Body %d in manual bodies file has type 0, which is "
                        "reserved for tree code; use a nonzero value "
                        "(negative for dark matter, positive for light "
                        "matter).\n", (int) id[i]);
            }
            b.bodynode.type = (body_t) ty[i];
        }
        else
        {
            /* Legacy format: 0 = baryon, 1 = dark matter; tree uses 1 / -1 */
            b.bodynode.type = BODY(ty[i]);
        }
        b.bodynode.id = id[i];
        b.bodynode.mass = masses[i];

        /*this actually gets the position and velocity vectors and pushes table of bodies*/
        /*They are meant to give the dwarf an initial position and vel*/
        /* you have to work for your bodynode */
        b.bodynode.pos.x = x[i];
        b.bodynode.pos.y = y[i];
        b.bodynode.pos.z = z[i];
        
        b.vel.x = vx[i];
        b.vel.y = vy[i];
        b.vel.z = vz[i];
        
//         mw_printf("%f %f %f %f %f %f %f\n", b.bodynode.pos.x, b.bodynode.pos.y, b.bodynode.pos.z, b.vel.x, b.vel.y, b.vel.z, b.bodynode.mass);
        assert(nbPositionValid(b.bodynode.pos));
        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
    }
    
    
    /* go now and be free!*/
    free(ty);
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
    free(masses);
    free(id);
    return 1;             
        
}

int nbGenerateManualBodies(lua_State* luaSt)
{
        static const char* body_file = NULL;
        static const MWNamedArg argTable[] =
        {
            { "body_file",            LUA_TSTRING,     NULL,                    TRUE,    &body_file,           1 },
            END_MW_NAMED_ARG
            
        };

        if (lua_gettop(luaSt) != 1)
            return luaL_argerror(luaSt, 1, "Expected 1 arguments");
        
        handleNamedArgumentTable(luaSt, argTable, 1);
        
        return nbGenerateManualBodiescore(luaSt, body_file);
}


void registerGenerateManualBodies(lua_State* luaSt)
{
    lua_register(luaSt, "generatemanualbodies", nbGenerateManualBodies);
}


// As this code runs, know that it is running on the rotting corpses of a thousand bugs.
