/*
Copyright (C) 2011  Matthew Arsenault
Copyright (c) 2016-2018 Siddhartha Shelton

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

#include "nbody_types.h"
#include "nbody_defaults.h"

const NBodyCtx defaultNBodyCtx =
{
    /* Grr lack of C99 named struct initializers in MSVC */
    /* .eps2            */  0.0,
    /* .theta           */  -1.0,  /* Invalid */

    /* .timestep        */  0.0,
    /* .timeEvolve      */  0.0,
    /* .timeBack        */  0.0,
    /* .treeRSize       */  DEFAULT_TREE_ROOT_SIZE,
    /* .sunGCDist       */  DEFAULT_SUN_GC_DISTANCE,
    /* .sunVelx         */  DEFAULT_SUN_VEL_X,
    /* .sunVely         */  DEFAULT_SUN_VEL_Y,
    /* .sunVelz         */  DEFAULT_SUN_VEL_Z,
    /* .NGPdec          */  DEFAULT_NGP_DEC,
    /* .NGPra           */  DEFAULT_NGP_RA,
    /* .lNCP            */  DEFAULT_L_NCP,

    /* .b               */  DEFAULT_B_START_COORD,
    /* .r               */  DEFAULT_R_START_COORD,
    /* .vx              */  DEFAULT_VX_START_COORD,
    /* .vy              */  DEFAULT_VY_START_COORD,
    /* .vz              */  DEFAULT_VZ_START_COORD,

    /* .criterion       */  DEFAULT_CRITERION,
    /* .potentialType   */  EXTERNAL_POTENTIAL_DEFAULT,

    /* .Nstep_control   */  FALSE,
    /* .useBestLike     */  DEFAULT_USE_BEST_LIKELIHOOD,
    /* .useVelDisp      */  DEFAULT_USE_VEL_DISP,
    /* .useBetaDisp     */  DEFAULT_USE_BETA_DISP,
    /* .useBetaComp     */  DEFAULT_USE_BETA_COMP,
    /* .useVlos         */  DEFAULT_USE_VLOS,
    /* .useDist         */  DEFAULT_USE_DIST, 
    /* .usePropMot      */  DEFAULT_USE_PROP_MOT,
    /* .MultiOutput     */  FALSE,
    /* .InitialOutput   */  FALSE,
    /* .SimpleOutput    */  TRUE,

    /* .useQuad         */  DEFAULT_USE_QUADRUPOLE_MOMENTS,
    /* .allowIncest     */  DEFAULT_ALLOW_INCEST,
    /* .quietErrors     */  DEFAULT_QUIET_ERRORS,

    /* .BestLikeStart   */  DEFAULT_BEST_LIKELIHOOD_START,
    /* .OutputFreq      */  DEFAULT_OUTPUT_FREQUENCY,
    
    /* .BetaSigma       */  DEFAULT_SIGMA_CUTOFF,
    /* .VelSigma        */  DEFAULT_SIGMA_CUTOFF,
    /* .DistSigma       */  DEFAULT_SIGMA_CUTOFF,
    /* .PMSigma         */  DEFAULT_SIGMA_CUTOFF,
    /* .IterMax         */  DEFAULT_SIGMA_ITER,
    /* .BetaCorrect     */  DEFAULT_DISP_CORRECTION,
    /* .VelCorrect      */  DEFAULT_DISP_CORRECTION,
    /* .DistCorrect     */  DEFAULT_DISP_CORRECTION,
    /* .PMCorrect       */  DEFAULT_DISP_CORRECTION,

    /* .LMC             */  FALSE,
    /* .LMCfunction     */  1,
    /* .LMCmass         */  0.0,
    /* .LMCscale        */  0.0,
    /* .LMCscale2       */  0.0,
    /* .LMCDynaFric     */  FALSE,
    /* .coulomb_log     */  0.0,

    /* .Ntsteps         */  0,
    /* .checkpointT     */  NOBOINC_DEFAULT_CHECKPOINT_PERIOD,
    /* .nStep           */  0,

    /* .calibrationRuns */   0,

    /* .pot             */  EMPTY_POTENTIAL
};

const HistogramParams defaultHistogramParams =
{
    /* .phi         */  histogramPhi,
    /* .theta       */  histogramTheta,
    /* .psi         */  histogramPsi,
    /* .lambdaStart */  histogramlambdaStart,
    /* .lambdaEnd   */  histogramlambdaEnd,
    /* .lambdaBins  */  histogramlambdaBins,
    /* .lambdaStart */  histogrambetaStart,
    /* .lambdaEnd   */  histogrambetaEnd,
    /* .lambdaBins  */  histogrambetaBins
};



