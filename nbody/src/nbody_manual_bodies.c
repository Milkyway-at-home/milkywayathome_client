/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
     Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright (c) 2016 Siddhartha Shelton
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

/*      DWARF GENERATION        */
static int nbGenerateManualBodiescore(lua_State* luaSt, unsigned int nbody, const char* body_file)
{
    /* generatePlummer: generate Plummer model initial conditions for test
    * runs, scaled to units such that M = -4E = G = 1 (Henon, Heggie,
    * etc).    See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
    * 183.
    */
        
//         fscanf(body_inputs, "%f %f %f %f %f %f %f %f", xx, yy, zz, hh, rr, tt, ee, ww);
//         mw_printf("%f\n", zz);
        
        unsigned int i;
        int table;
        Body b;
        real r, v;
 
        real * x  = mwCalloc(nbody, sizeof(real));
        real * y  = mwCalloc(nbody, sizeof(real));
        real * z  = mwCalloc(nbody, sizeof(real));
        real * vx = mwCalloc(nbody, sizeof(real));
        real * vy = mwCalloc(nbody, sizeof(real));
        real * vz = mwCalloc(nbody, sizeof(real));
        real * masses = mwCalloc(nbody, sizeof(real));
        
        FILE* body_inputs;
        body_inputs = fopen(body_file, "r");
//         real tty, xx, yy, zz, vxx, vyy, vzz, mm;
        char buf1[100];
        char buf2[100];
        char buf3[100];
        char buf4[100];
        char buf5[100];
        char buf6[100];
        char buf7[100];
        char buf8[100];
//         double buf1;
//         double buf2;
//         double buf3;
//         double buf4;
//         double buf5;
//         double buf6;
//         double buf7;
//         double buf8;
        
//         char buf1;
//         char buf2;
//         char buf3;
//         char buf4;
//         char buf5;
//         char buf6;
//         char buf7;
//         char buf8;
        unsigned int counter = 0;
        while (fscanf(body_inputs,"%s %s %s %s %s %s %s %s ",buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8)==8)
//         while (fscanf(body_inputs,"%f %f %f %f %f %f %f %f ",buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8)==8)
        {
//             x[counter]      = strtod(buf2, NULL);
//             y[counter]      = strtod(buf3, NULL);
//             z[counter]      = strtod(buf4, NULL);
//             vx[counter]     = strtod(buf5, NULL);
//             vy[counter]     = strtod(buf6, NULL);
//             vz[counter]     = strtod(buf7, NULL);
//             masses[counter] = strtod(buf8, NULL);
            
            x[counter]      = atof(buf2);
            y[counter]      = atof(buf3);
            z[counter]      = atof(buf4);
            vx[counter]     = atof(buf5);
            vy[counter]     = atof(buf6);
            vz[counter]     = atof(buf7);
            masses[counter] = atof(buf8);
            mw_printf("%s %s %s %s %s %s %s %s\n", buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8);
//             mw_printf("%f %f %f %f %f %f %f %f\n", buf1, buf2, buf3, buf4, buf5, buf6, buf7, buf8);
            counter++;
//             mw_printf("%f %f %f %f %f %f %f\n", x[counter], y[counter], z[counter], vx[counter], vy[counter], vz[counter], masses[counter]);
        }
        
        fclose(body_inputs);
        mwvector vec;
        mwbool isdark = TRUE;
        mwbool islight = FALSE;
       
        
     /*initializing particles:*/
        memset(&b, 0, sizeof(b));
        lua_createtable(luaSt, nbody, 0);
        table = lua_gettop(luaSt);      
//         int counter = 0;
        

        /*getting the radii and velocities for the bodies*/

        /* pushing the bodies */
        for (i = 0; i < nbody; i++)
        {
            b.bodynode.type = BODY(islight);
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
            
            assert(nbPositionValid(b.bodynode.pos));
            pushBody(luaSt, &b);
            lua_rawseti(luaSt, table, i + 1);
        }
        
        /* go now and be free!*/
        free(x);
        free(y);
        free(z);
        free(vx);
        free(vy);
        free(vz);
        free(masses);
        
        return 1;             
        
}

int nbGenerateManualBodies(lua_State* luaSt)
{
        static const char* body_file = NULL;
        static real nbodyf = 0.0;
        static const MWNamedArg argTable[] =
        {
            { "nbody",                LUA_TNUMBER,     NULL,                    TRUE,    &nbodyf            },
            { "body_file",            LUA_TSTRING,     NULL,                    TRUE,    &body_file         },
            END_MW_NAMED_ARG
            
        };

        if (lua_gettop(luaSt) != 1)
            return luaL_argerror(luaSt, 1, "Expected 1 arguments");
        
        handleNamedArgumentTable(luaSt, argTable, 1);
        
        return nbGenerateManualBodiescore(luaSt, (unsigned int) nbodyf, body_file);
}


void registerGenerateManualBodies(lua_State* luaSt)
{
    lua_register(luaSt, "generatemanualbodies", nbGenerateManualBodies);
}


// As this code runs, know that it is running on the rotting corpses of a thousand bugs.