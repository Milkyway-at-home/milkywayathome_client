/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "nbody_priv.h"
#include "nbody_orbit_integrator.h"
#include "nbody_potential.h"
#include "nbody_io.h"
#include "nbody_coordinates.h"
#include "nbody_defaults.h"
#include "milkyway_util.h"
#include "nbody.h"
#include "nbody_friction.h"
/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
/* Altered to be consistent with nbody integrator.
 * shelton June 25 2018 */

mwvector* shiftByLMC = NULL; //Ptr to LMC Shift Array (default is NULL)
size_t nShiftLMC = 0;

mwvector LMCpos = ZERO_VECTOR; //Ptr to LMC position (default is NULL)
mwvector LMCvel = ZERO_VECTOR; //Ptr to LMC velocity (default is NULL)

#if AUTODIFF
  static inline void getTimeDerivativeInfo(const Potential* pot, mwvector* dwarfPos, mwvector* dwarfVel, mwvector* LMCPos, mwvector* LMCVel, real* LMCMass, real* LMCScale, mwbool dynaFric, real_0 time)
  {
      mwvector mw_origin = ZERO_VECTOR;

      mwvector accel_dwarf_by_MW = nbExtAcceleration(pot, dwarfPos, time);

      /* Calculate jerk of dwarf by MW potential */
      mwvector shift_dwarf_x = *dwarfPos;
      mwvector shift_dwarf_y = *dwarfPos;
      mwvector shift_dwarf_z = *dwarfPos;
      shift_dwarf_x.x = mw_mul(&shift_dwarf_x.x, 1.01);
      shift_dwarf_y.y = mw_mul(&shift_dwarf_y.y, 1.01);
      shift_dwarf_z.z = mw_mul(&shift_dwarf_z.z, 1.01);
      mwvector accel_dwarf_by_MW_x = nbExtAcceleration(pot, &shift_dwarf_x, time);
      mwvector accel_dwarf_by_MW_y = nbExtAcceleration(pot, &shift_dwarf_y, time);
      mwvector accel_dwarf_by_MW_z = nbExtAcceleration(pot, &shift_dwarf_z, time);

      real_0 dax_dx = (accel_dwarf_by_MW_x.x.value - accel_dwarf_by_MW.x.value)/0.01/showRealValue(&dwarfPos->x);
      real_0 day_dx = (accel_dwarf_by_MW_x.y.value - accel_dwarf_by_MW.y.value)/0.01/showRealValue(&dwarfPos->x);
      real_0 daz_dx = (accel_dwarf_by_MW_x.z.value - accel_dwarf_by_MW.z.value)/0.01/showRealValue(&dwarfPos->x);

      real_0 dax_dy = (accel_dwarf_by_MW_y.x.value - accel_dwarf_by_MW.x.value)/0.01/showRealValue(&dwarfPos->y);
      real_0 day_dy = (accel_dwarf_by_MW_y.y.value - accel_dwarf_by_MW.y.value)/0.01/showRealValue(&dwarfPos->y);
      real_0 daz_dy = (accel_dwarf_by_MW_y.z.value - accel_dwarf_by_MW.z.value)/0.01/showRealValue(&dwarfPos->y);

      real_0 dax_dz = (accel_dwarf_by_MW_z.x.value - accel_dwarf_by_MW.x.value)/0.01/showRealValue(&dwarfPos->z);
      real_0 day_dz = (accel_dwarf_by_MW_z.y.value - accel_dwarf_by_MW.y.value)/0.01/showRealValue(&dwarfPos->z);
      real_0 daz_dz = (accel_dwarf_by_MW_z.z.value - accel_dwarf_by_MW.z.value)/0.01/showRealValue(&dwarfPos->z);

      real_0 jerk_dwarf_x = dax_dx * (dwarfVel->x.value) + dax_dy * (dwarfVel->y.value) + dax_dz * (dwarfVel->z.value);
      real_0 jerk_dwarf_y = day_dx * (dwarfVel->x.value) + day_dy * (dwarfVel->y.value) + day_dz * (dwarfVel->z.value);
      real_0 jerk_dwarf_z = daz_dx * (dwarfVel->x.value) + daz_dy * (dwarfVel->y.value) + daz_dz * (dwarfVel->z.value);

      /* Modify first order time derivative for dwarf position */
      setRealGradient(&dwarfPos->x, -(dwarfVel->x.value), BACKWARDS_TIME_POS);
      setRealGradient(&dwarfPos->y, -(dwarfVel->y.value), BACKWARDS_TIME_POS);
      setRealGradient(&dwarfPos->z, -(dwarfVel->z.value), BACKWARDS_TIME_POS);

      if(LMCPos)
      {
          /* Modify dwarf jerk with jerk from LMC */
          mwvector accel_dwarf_by_LMC = plummerAccel(dwarfPos, LMCPos, LMCmass, LMCscale);
          mwvector accel_dwarf_by_LMC_x = plummerAccel(&shift_dwarf_x, LMCPos, LMCmass, LMCscale);
          mwvector accel_dwarf_by_LMC_y = plummerAccel(&shift_dwarf_y, LMCPos, LMCmass, LMCscale);
          mwvector accel_dwarf_by_LMC_z = plummerAccel(&shift_dwarf_z, LMCPos, LMCmass, LMCscale);

          dax_dx = (accel_dwarf_by_LMC_x.x.value - accel_dwarf_by_LMC.x.value)/0.01/showRealValue(&dwarfPos->x);
          day_dx = (accel_dwarf_by_LMC_x.y.value - accel_dwarf_by_LMC.y.value)/0.01/showRealValue(&dwarfPos->x);
          daz_dx = (accel_dwarf_by_LMC_x.z.value - accel_dwarf_by_LMC.z.value)/0.01/showRealValue(&dwarfPos->x);

          dax_dy = (accel_dwarf_by_LMC_y.x.value - accel_dwarf_by_LMC.x.value)/0.01/showRealValue(&dwarfPos->y);
          day_dy = (accel_dwarf_by_LMC_y.y.value - accel_dwarf_by_LMC.y.value)/0.01/showRealValue(&dwarfPos->y);
          daz_dy = (accel_dwarf_by_LMC_y.z.value - accel_dwarf_by_LMC.z.value)/0.01/showRealValue(&dwarfPos->y);

          dax_dz = (accel_dwarf_by_LMC_z.x.value - accel_dwarf_by_LMC.x.value)/0.01/showRealValue(&dwarfPos->z);
          day_dz = (accel_dwarf_by_LMC_z.y.value - accel_dwarf_by_LMC.y.value)/0.01/showRealValue(&dwarfPos->z);
          daz_dz = (accel_dwarf_by_LMC_z.z.value - accel_dwarf_by_LMC.z.value)/0.01/showRealValue(&dwarfPos->z);

          jerk_dwarf_x += dax_dx * (dwarfVel->x.value - LMCVel->x.value) + dax_dy * (dwarfVel->y.value - LMCVel->y.value) + dax_dz * (dwarfVel->z.value - LMCVel->z.value);
          jerk_dwarf_y += day_dx * (dwarfVel->x.value - LMCVel->x.value) + day_dy * (dwarfVel->y.value - LMCVel->y.value) + day_dz * (dwarfVel->z.value - LMCVel->z.value);
          jerk_dwarf_z += daz_dx * (dwarfVel->x.value - LMCVel->x.value) + daz_dy * (dwarfVel->y.value - LMCVel->y.value) + daz_dz * (dwarfVel->z.value - LMCVel->z.value);

          /* Modify dwarf jerk with jerk from MW shift */
          mwvector accel_MW_by_LMC = plummerAccel(&mw_origin, LMCPos, LMCmass, LMCscale);
          mwvector shift_LMC_x = *LMCPos;
          mwvector shift_LMC_y = *LMCPos;
          mwvector shift_LMC_z = *LMCPos;
          shift_LMC_x.x = mw_mul(&shift_LMC_x.x, 1.01);
          shift_LMC_y.y = mw_mul(&shift_LMC_y.y, 1.01);
          shift_LMC_z.z = mw_mul(&shift_LMC_z.z, 1.01);
          mwvector accel_MW_by_LMC_x = plummerAccel(&mw_origin, &shift_LMC_x, LMCmass, LMCscale);
          mwvector accel_MW_by_LMC_y = plummerAccel(&mw_origin, &shift_LMC_y, LMCmass, LMCscale);
          mwvector accel_MW_by_LMC_z = plummerAccel(&mw_origin, &shift_LMC_z, LMCmass, LMCscale);

          dax_dx = (accel_MW_by_LMC_x.x.value - accel_MW_by_LMC.x.value)/0.01/showRealValue(&LMCPos->x);
          day_dx = (accel_MW_by_LMC_x.y.value - accel_MW_by_LMC.y.value)/0.01/showRealValue(&LMCPos->x);
          daz_dx = (accel_MW_by_LMC_x.z.value - accel_MW_by_LMC.z.value)/0.01/showRealValue(&LMCPos->x);

          dax_dy = (accel_MW_by_LMC_y.x.value - accel_MW_by_LMC.x.value)/0.01/showRealValue(&LMCPos->y);
          day_dy = (accel_MW_by_LMC_y.y.value - accel_MW_by_LMC.y.value)/0.01/showRealValue(&LMCPos->y);
          daz_dy = (accel_MW_by_LMC_y.z.value - accel_MW_by_LMC.z.value)/0.01/showRealValue(&LMCPos->y);

          dax_dz = (accel_MW_by_LMC_z.x.value - accel_MW_by_LMC.x.value)/0.01/showRealValue(&LMCPos->z);
          day_dz = (accel_MW_by_LMC_z.y.value - accel_MW_by_LMC.y.value)/0.01/showRealValue(&LMCPos->z);
          daz_dz = (accel_MW_by_LMC_z.z.value - accel_MW_by_LMC.z.value)/0.01/showRealValue(&LMCPos->z);

          jerk_dwarf_x -= dax_dx * (LMCVel->x.value) + dax_dy * (LMCVel->y.value) + dax_dz * (LMCVel->z.value);
          jerk_dwarf_y -= day_dx * (LMCVel->x.value) + day_dy * (LMCVel->y.value) + day_dz * (LMCVel->z.value);
          jerk_dwarf_z -= daz_dx * (LMCVel->x.value) + daz_dy * (LMCVel->y.value) + daz_dz * (LMCVel->z.value);

          /* This term is also needed for LMC jerk */
          real_0 jerk_LMC_x = -(dax_dx * (LMCVel->x.value) + dax_dy * (LMCVel->y.value) + dax_dz * (LMCVel->z.value));
          real_0 jerk_LMC_y = -(day_dx * (LMCVel->x.value) + day_dy * (LMCVel->y.value) + day_dz * (LMCVel->z.value));
          real_0 jerk_LMC_z = -(daz_dx * (LMCVel->x.value) + daz_dy * (LMCVel->y.value) + daz_dz * (LMCVel->z.value));

          mwvector total_accel_dwarf = mw_addv(&accel_dwarf_by_MW, &accel_dwarf_by_LMC);
          total_accel_dwarf = mw_subv(&total_accel_dwarf, &accel_MW_by_LMC);

          /* Calculate total LMC acceleration */
          mwvector accel_LMC_by_MW = nbExtAcceleration(pot, LMCPos, time);
          mwvector DF_acc = dynamicalFriction_LMC(pot, LMCPos, LMCVel, LMCMass, LMCScale, dynaFric, time);
          accel_LMC_by_MW = mw_subv(&accel_LMC_by_MW, &DF_acc); /* Inverting drag force because of time reversal */
          mwvector total_accel_LMC = mw_subv(&accel_LMC_by_MW, &accel_MW_by_LMC);

          /* Calculate jerk of LMC by MW */
          mwvector accel_LMC_by_MW_x = nbExtAcceleration(pot, &shift_LMC_x, time);
          mwvector accel_LMC_by_MW_y = nbExtAcceleration(pot, &shift_LMC_y, time);
          mwvector accel_LMC_by_MW_z = nbExtAcceleration(pot, &shift_LMC_z, time);

          dax_dx = (accel_LMC_by_MW_x.x.value - accel_LMC_by_MW.x.value)/0.01/showRealValue(&LMCPos->x);
          day_dx = (accel_LMC_by_MW_x.y.value - accel_LMC_by_MW.y.value)/0.01/showRealValue(&LMCPos->x);
          daz_dx = (accel_LMC_by_MW_x.z.value - accel_LMC_by_MW.z.value)/0.01/showRealValue(&LMCPos->x);

          dax_dy = (accel_LMC_by_MW_y.x.value - accel_LMC_by_MW.x.value)/0.01/showRealValue(&LMCPos->y);
          day_dy = (accel_LMC_by_MW_y.y.value - accel_LMC_by_MW.y.value)/0.01/showRealValue(&LMCPos->y);
          daz_dy = (accel_LMC_by_MW_y.z.value - accel_LMC_by_MW.z.value)/0.01/showRealValue(&LMCPos->y);

          dax_dz = (accel_LMC_by_MW_z.x.value - accel_LMC_by_MW.x.value)/0.01/showRealValue(&LMCPos->z);
          day_dz = (accel_LMC_by_MW_z.y.value - accel_LMC_by_MW.y.value)/0.01/showRealValue(&LMCPos->z);
          daz_dz = (accel_LMC_by_MW_z.z.value - accel_LMC_by_MW.z.value)/0.01/showRealValue(&LMCPos->z);

          jerk_LMC_x += dax_dx * (LMCVel->x.value) + dax_dy * (LMCVel->y.value) + dax_dz * (LMCVel->z.value);
          jerk_LMC_y += day_dx * (LMCVel->x.value) + day_dy * (LMCVel->y.value) + day_dz * (LMCVel->z.value);
          jerk_LMC_z += daz_dx * (LMCVel->x.value) + daz_dy * (LMCVel->y.value) + daz_dz * (LMCVel->z.value);

          /* Calculate jerk of LMC by Dynamical Friction (position term) */
          mwvector DF_acc_x = dynamicalFriction_LMC(pot, &shift_LMC_x, LMCVel, LMCMass, LMCScale, dynaFric, time);
          mwvector DF_acc_y = dynamicalFriction_LMC(pot, &shift_LMC_y, LMCVel, LMCMass, LMCScale, dynaFric, time);
          mwvector DF_acc_z = dynamicalFriction_LMC(pot, &shift_LMC_z, LMCVel, LMCMass, LMCScale, dynaFric, time);

          dax_dx = (DF_acc_x.x.value - DF_acc.x.value)/0.01/showRealValue(&LMCPos->x);
          day_dx = (DF_acc_x.y.value - DF_acc.y.value)/0.01/showRealValue(&LMCPos->x);
          daz_dx = (DF_acc_x.z.value - DF_acc.z.value)/0.01/showRealValue(&LMCPos->x);

          dax_dy = (DF_acc_y.x.value - DF_acc.x.value)/0.01/showRealValue(&LMCPos->y);
          day_dy = (DF_acc_y.y.value - DF_acc.y.value)/0.01/showRealValue(&LMCPos->y);
          daz_dy = (DF_acc_y.z.value - DF_acc.z.value)/0.01/showRealValue(&LMCPos->y);

          dax_dz = (DF_acc_z.x.value - DF_acc.x.value)/0.01/showRealValue(&LMCPos->z);
          day_dz = (DF_acc_z.y.value - DF_acc.y.value)/0.01/showRealValue(&LMCPos->z);
          daz_dz = (DF_acc_z.z.value - DF_acc.z.value)/0.01/showRealValue(&LMCPos->z);

          jerk_LMC_x -= dax_dx * (LMCVel->x.value) + dax_dy * (LMCVel->y.value) + dax_dz * (LMCVel->z.value);
          jerk_LMC_y -= day_dx * (LMCVel->x.value) + day_dy * (LMCVel->y.value) + day_dz * (LMCVel->z.value);
          jerk_LMC_z -= daz_dx * (LMCVel->x.value) + daz_dy * (LMCVel->y.value) + daz_dz * (LMCVel->z.value);

          /* Calculate jerk of LMC by Dynamical Friction (velocity term) */
          mwvector shift_LMC_vx = *LMCVel;
          mwvector shift_LMC_vy = *LMCVel;
          mwvector shift_LMC_vz = *LMCVel;
          shift_LMC_vx.x = mw_mul(&shift_LMC_vx.x, 1.01);
          shift_LMC_vy.y = mw_mul(&shift_LMC_vy.y, 1.01);
          shift_LMC_vz.z = mw_mul(&shift_LMC_vz.z, 1.01);

          DF_acc_x = dynamicalFriction_LMC(pot, LMCPos, &shift_LMC_x, LMCMass, LMCScale, dynaFric, time);
          DF_acc_y = dynamicalFriction_LMC(pot, LMCPos, &shift_LMC_y, LMCMass, LMCScale, dynaFric, time);
          DF_acc_z = dynamicalFriction_LMC(pot, LMCPos, &shift_LMC_z, LMCMass, LMCScale, dynaFric, time);

          dax_dx = (DF_acc_x.x.value - DF_acc.x.value)/0.01/showRealValue(&LMCVel->x);
          day_dx = (DF_acc_x.y.value - DF_acc.y.value)/0.01/showRealValue(&LMCVel->x);
          daz_dx = (DF_acc_x.z.value - DF_acc.z.value)/0.01/showRealValue(&LMCVel->x);

          dax_dy = (DF_acc_y.x.value - DF_acc.x.value)/0.01/showRealValue(&LMCVel->y);
          day_dy = (DF_acc_y.y.value - DF_acc.y.value)/0.01/showRealValue(&LMCVel->y);
          daz_dy = (DF_acc_y.z.value - DF_acc.z.value)/0.01/showRealValue(&LMCVel->y);

          dax_dz = (DF_acc_z.x.value - DF_acc.x.value)/0.01/showRealValue(&LMCVel->z);
          day_dz = (DF_acc_z.y.value - DF_acc.y.value)/0.01/showRealValue(&LMCVel->z);
          daz_dz = (DF_acc_z.z.value - DF_acc.z.value)/0.01/showRealValue(&LMCVel->z);

          jerk_LMC_x -= dax_dx * (total_accel_LMC.x.value) + dax_dy * (total_accel_LMC.y.value) + dax_dz * (total_accel_LMC.z.value);
          jerk_LMC_y -= day_dx * (total_accel_LMC.x.value) + day_dy * (total_accel_LMC.y.value) + day_dz * (total_accel_LMC.z.value);
          jerk_LMC_z -= daz_dx * (total_accel_LMC.x.value) + daz_dy * (total_accel_LMC.y.value) + daz_dz * (total_accel_LMC.z.value);

          /* Modify first order time derivative for dwarf velocity */
          setRealGradient(&dwarfVel->x, -(total_accel_dwarf.x.value), BACKWARDS_TIME_POS);
          setRealGradient(&dwarfVel->y, -(total_accel_dwarf.y.value), BACKWARDS_TIME_POS);
          setRealGradient(&dwarfVel->z, -(total_accel_dwarf.z.value), BACKWARDS_TIME_POS);

          /* Modify first order time derivative for LMC position */
          setRealGradient(&LMCPos->x, -(LMCVel->x.value), BACKWARDS_TIME_POS);
          setRealGradient(&LMCPos->y, -(LMCVel->y.value), BACKWARDS_TIME_POS);
          setRealGradient(&LMCPos->z, -(LMCVel->z.value), BACKWARDS_TIME_POS);

          /* Modify first order time derivative for LMC velocity */
          setRealGradient(&LMCPos->x, -(total_accel_LMC.x.value), BACKWARDS_TIME_POS);
          setRealGradient(&LMCPos->y, -(total_accel_LMC.y.value), BACKWARDS_TIME_POS);
          setRealGradient(&LMCPos->z, -(total_accel_LMC.z.value), BACKWARDS_TIME_POS);

          /* Set all other second order time derivatives */
          for(int i = 0; i < NumberOfModelParameters; i++)
          {
              if(i == BACKWARDS_TIME_POS)
              {
                  setRealHessian(&dwarfPos->x, total_accel_dwarf.x.value, BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfPos->y, total_accel_dwarf.y.value, BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfPos->z, total_accel_dwarf.z.value, BACKWARDS_TIME_POS, i);

                  setRealHessian(&dwarfVel->x, jerk_dwarf_x, BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfVel->y, jerk_dwarf_y, BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfVel->z, jerk_dwarf_z, BACKWARDS_TIME_POS, i);

                  setRealHessian(&LMCPos->x, total_accel_LMC.x.value, BACKWARDS_TIME_POS, i);
                  setRealHessian(&LMCPos->y, total_accel_LMC.y.value, BACKWARDS_TIME_POS, i);
                  setRealHessian(&LMCPos->z, total_accel_LMC.z.value, BACKWARDS_TIME_POS, i);

                  setRealHessian(&LMCVel->x, jerk_LMC_x, BACKWARDS_TIME_POS, i);
                  setRealHessian(&LMCVel->y, jerk_LMC_y, BACKWARDS_TIME_POS, i);
                  setRealHessian(&LMCVel->z, jerk_LMC_z, BACKWARDS_TIME_POS, i);
              }
              else
              {
                  setRealHessian(&dwarfPos->x, -(dwarfVel->x.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfPos->y, -(dwarfVel->y.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfPos->z, -(dwarfVel->z.gradient[i]), BACKWARDS_TIME_POS, i);

                  setRealHessian(&dwarfVel->x, -(total_accel_dwarf.x.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfVel->y, -(total_accel_dwarf.y.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfVel->z, -(total_accel_dwarf.z.gradient[i]), BACKWARDS_TIME_POS, i);

                  setRealHessian(&LMCPos->x, -(LMCVel->x.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&LMCPos->y, -(LMCVel->y.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&LMCPos->z, -(LMCVel->z.gradient[i]), BACKWARDS_TIME_POS, i);

                  setRealHessian(&LMCVel->x, -(total_accel_LMC.x.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&LMCVel->y, -(total_accel_LMC.y.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&LMCVel->z, -(total_accel_LMC.z.gradient[i]), BACKWARDS_TIME_POS, i);
              }
          }
      }
      else  //FOR NO LMC
      {
          /* Modify first order time derivative for dwarf velocity */
          setRealGradient(&dwarfVel->x, -(accel_dwarf_by_MW.x.value), BACKWARDS_TIME_POS);
          setRealGradient(&dwarfVel->y, -(accel_dwarf_by_MW.y.value), BACKWARDS_TIME_POS);
          setRealGradient(&dwarfVel->z, -(accel_dwarf_by_MW.z.value), BACKWARDS_TIME_POS);

          /* Set all other second order time derivatives */
          for(int i = 0; i < NumberOfModelParameters; i++)
          {
              if(i == BACKWARDS_TIME_POS)
              {
                  setRealHessian(&dwarfPos->x, accel_dwarf_by_MW.x.value, BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfPos->y, accel_dwarf_by_MW.y.value, BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfPos->z, accel_dwarf_by_MW.z.value, BACKWARDS_TIME_POS, i);

                  setRealHessian(&dwarfVel->x, jerk_dwarf_x, BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfVel->y, jerk_dwarf_y, BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfVel->z, jerk_dwarf_z, BACKWARDS_TIME_POS, i);
              }
              else
              {
                  setRealHessian(&dwarfPos->x, -(dwarfVel->x.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfPos->y, -(dwarfVel->y.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfPos->z, -(dwarfVel->z.gradient[i]), BACKWARDS_TIME_POS, i);

                  setRealHessian(&dwarfVel->x, -(accel_dwarf_by_MW.x.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfVel->y, -(accel_dwarf_by_MW.y.gradient[i]), BACKWARDS_TIME_POS, i);
                  setRealHessian(&dwarfVel->z, -(accel_dwarf_by_MW.z.gradient[i]), BACKWARDS_TIME_POS, i);
              }
          }
      }
  }
#endif //AUTODIFF

void nbReverseOrbit(mwvector* finalPos,
                    mwvector* finalVel,
                    const Potential* pot,
                    mwvector* pos,
                    mwvector* vel,
                    real_0 tstop,
                    real_0 dt,
                    real_0 sun_dist)
{
    mw_printf("Performing Reverse Orbit Calculation...\n");
    mwvector v_var = ZERO_VECTOR;
    mwvector x_var = ZERO_VECTOR;
    mwvector x_lbr = ZERO_VECTOR;
    mwvector dv2 = ZERO_VECTOR;
    mwvector dx = ZERO_VECTOR;
    mwvector acc = ZERO_VECTOR;
    mwvector v = ZERO_VECTOR;
    mwvector x = ZERO_VECTOR;
    real_0 t = 0.0;
    real_0 dt_half = dt / 2.0;
    int initialLArrayIndex = tstop/dt;

    // Set derivative information for AUTODIFF
    //mw_printf("OLD X = [%.15f, %.15f, %.15f]\n", showRealValue(&pos->x), showRealValue(&pos->y), showRealValue(&pos->z));
    x_lbr = cartesianToLbr(pos, sun_dist);
    B(&x_lbr) = mw_real_var(showRealValue(&B(&x_lbr)), B_COORD_POS);
    R(&x_lbr) = mw_real_var(showRealValue(&R(&x_lbr)), R_COORD_POS);
    x_var = lbrToCartesian(&x_lbr, sun_dist);
    //mw_printf("NEW X = [%.15f, %.15f, %.15f]\n", showRealValue(&x_var.x), showRealValue(&x_var.y), showRealValue(&x_var.z));

    //mw_printf("OLD V = [%.15f, %.15f, %.15f]\n", showRealValue(&vel->x), showRealValue(&vel->y), showRealValue(&vel->z));
    v_var.x = mw_real_var(showRealValue(&X(vel)), VX_COORD_POS);
    v_var.y = mw_real_var(showRealValue(&Y(vel)), VY_COORD_POS);
    v_var.z = mw_real_var(showRealValue(&Z(vel)), VZ_COORD_POS);
    //mw_printf("NEW V = [%.15f, %.15f, %.15f]\n", showRealValue(&v_var.x), showRealValue(&v_var.y), showRealValue(&v_var.z));

    // Set the initial conditions
    x = x_var;
    v = mw_negv(&v_var);

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, &x, 0);

    //do this loop backward in order to get an accurate time for time-dependent potentials
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        //mw_printf("ACC  = [%.15f, %.15f, %.15f]\n", showRealValue(&acc.x), showRealValue(&acc.y), showRealValue(&acc.z));
        //mw_printf("POS  = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
        //mw_printf("VEL  = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));
        // Update the velocities and positions
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        //mw_printf("DV/2 = [%.15f, %.15f, %.15f]\n", showRealValue(&dv2.x), showRealValue(&dv2.y), showRealValue(&dv2.z));
        v = mw_addv(&v, &dv2);

        dx.x = mw_mul_s(&v.x, dt);
        dx.y = mw_mul_s(&v.y, dt);
        dx.z = mw_mul_s(&v.z, dt);
        //mw_printf("DX   = [%.15f, %.15f, %.15f]\n", showRealValue(&dx.x), showRealValue(&dx.y), showRealValue(&dx.z));
        x = mw_addv(&x, &dx); 

        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, &x, t);
        
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        //mw_printf("DV/2 = [%.15f, %.15f, %.15f]\n", showRealValue(&dv2.x), showRealValue(&dv2.y), showRealValue(&dv2.z));
        v = mw_addv(&v, &dv2);
    }
    
    /* Report the final values (don't forget to reverse the velocities) */
    v = mw_negv(&v);
    //mw_printf("FINAL POS  = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
    //mw_printf("FINAL VEL  = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));

  #if AUTODIFF
    /* Insert reverse time derivative information into positions and velocities */
    getTimeDerivativeInfo(pot, &x, &v, NULL, NULL, NULL, NULL, FALSE, t);
  #endif
    *finalPos = x;
    *finalVel = v;

    mw_printf("Initial Dwarf Position = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
    mw_printf("Initial Dwarf Velocity = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));
}

/*WARNING: This function is very tempermental in AUTODIFF.
  If you make any changes here, be sure to check that the
  outputs of the AUTODIFF version match that of the one
  without AUTODIFF. It could be an alignment problem, but
  I'm not entirely sure what causes the differences. The
  code so far does match, so only make changes here if you
  absolutely need to. I found that using declaring multiple
  vectors to hold intermediate values seems to fix it.
  (See <<< code below)*/
void nbReverseOrbit_LMC(mwvector* finalPos,
                    mwvector* finalVel,
                    mwvector* LMCfinalPos,
                    mwvector* LMCfinalVel,
                    const Potential* pot,
                    mwvector* pos,
                    mwvector* vel,
                    mwvector* LMCposition,
                    mwvector* LMCvelocity,
                    mwbool LMCDynaFric,
                    real_0 ftime,
                    real_0 tstop,
                    real_0 dt,
                    real* LMCmass,
                    real* LMCscale,
                    real_0 sun_dist)
{
    mw_printf("Performing Reverse Orbit Calculation with LMC...\n");	
    unsigned int steps = mw_ceil_0((tstop)/(10*dt)) + 1;
    unsigned int exSteps = mw_abs_0(mw_ceil_0((ftime-tstop)/(10*dt)) + 1);
    unsigned int maxSteps = MAX(steps + 1, exSteps + 3);
    unsigned int i = 0, j = 0, k = 0;
    real_0 dt_half = dt / 2.0;
    mwvector v_var, x_var, x_lbr;
    mwvector acc, v, x, v_for, x_for, mw_acc, LMC_acc, DF_acc, LMCv, LMCx, LMCv_for, LMCx_for, tmp, dv2, dx, dLMCv2, dLMCx, dv2_for, dx_for, dLMCv2_for, dLMCx_for; //<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    mwvector mw_x = ZERO_VECTOR;
    mwvector* ArrayPlaceholder = NULL;

    unsigned int vectorSize = sizeof(mwvector);

    //Placeholder arrays for LMC acceleration corrections
    ArrayPlaceholder = (mwvector*)mwCallocA(maxSteps, vectorSize);
    //Allocate memory for the shift array equal to (x,y,z) i times with extra wiggle room dependent on evolve time
    shiftByLMC = (mwvector*)mwCallocA(steps + exSteps + 4, vectorSize); 

    real_0 t;

    // Set derivative information for AUTODIFF
    x_lbr = cartesianToLbr(pos, sun_dist);
    B(&x_lbr) = mw_real_var(showRealValue(&B(&x_lbr)), B_COORD_POS);
    R(&x_lbr) = mw_real_var(showRealValue(&R(&x_lbr)), R_COORD_POS);
    x_var = lbrToCartesian(&x_lbr, sun_dist);

    v_var.x = mw_real_var(showRealValue(&X(vel)), VX_COORD_POS);
    v_var.y = mw_real_var(showRealValue(&Y(vel)), VY_COORD_POS);
    v_var.z = mw_real_var(showRealValue(&Z(vel)), VZ_COORD_POS);

    // Set the initial conditions for reverse orbit
    mw_printf("    Calculating backward orbit...\n");
    x = x_var;
    v = mw_negv(&v_var);
    LMCv = mw_negv(LMCvelocity);
    LMCx = *LMCposition;


    // Get the initial acceleration
    LMC_acc = nbExtAcceleration(pot, &LMCx, 0);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, TRUE, 0);
        LMC_acc = mw_subv(&LMC_acc, &DF_acc); /* Inverting drag force for reverse orbit */
     }
    acc = nbExtAcceleration(pot, &x, 0);
    tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
    acc = mw_addv(&acc, &tmp);

    // Shift the body
    mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
    LMC_acc = mw_subv(&LMC_acc, &mw_acc);
    acc = mw_subv(&acc, &mw_acc);

    real_0 negT = 0;
    for (t = 0; t <= tstop; t += dt)
    {   
        //negate this time for use in time-dependent potentials
        negT = t*(-1);
    	steps = (int) mw_round_0(t/dt);
    	if( steps % 10 == 0){ 
    		ArrayPlaceholder[i] = mw_negv(&mw_acc);
        	i++;
    	}

        //mw_printf("ACC  = [%.15f, %.15f, %.15f]\n", showRealValue(&acc.x), showRealValue(&acc.y), showRealValue(&acc.z));
        //mw_printf("POS  = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
        //mw_printf("VEL  = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));

        // Update the velocities and positions
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v = mw_addv(&v, &dv2);

        dx.x = mw_mul_s(&v.x, dt);
        dx.y = mw_mul_s(&v.y, dt);
        dx.z = mw_mul_s(&v.z, dt);
        x = mw_addv(&x, &dx);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv = mw_addv(&LMCv, &dLMCv2);

        dLMCx.x = mw_mul_s(&LMCv.x, dt);
        dLMCx.y = mw_mul_s(&LMCv.y, dt);
        dLMCx.z = mw_mul_s(&LMCv.z, dt);
        LMCx = mw_addv(&LMCx, &dLMCx);
        
        // Compute the new acceleration
        LMC_acc = nbExtAcceleration(pot, &LMCx, negT);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, TRUE, negT);
            LMC_acc = mw_subv(&LMC_acc, &DF_acc); /* Inverting drag force for reverse orbit */
        }
        acc = nbExtAcceleration(pot, &x, negT);
        tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
        acc = mw_addv(&acc, &tmp);

    	// Shift the body
        mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
        LMC_acc = mw_subv(&LMC_acc, &mw_acc);
        acc = mw_subv(&acc, &mw_acc);
        
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v = mw_addv(&v, &dv2);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv = mw_addv(&LMCv, &dLMCv2);

    }
    ArrayPlaceholder[i] = mw_negv(&mw_acc); //set the last index after the loop ends

    //Fill reverse orbit of shift array
    for(j = 0; j < i+1; j++) {
        tmp = ArrayPlaceholder[i-j];
        shiftByLMC[j] = tmp;
    }

    // Check if forward time is larger than backward time. We will need to manually compute additional LMC accelerations in that case.
    if (ftime >= tstop)
    {
        // Set the initial conditions for forward orbit
        mw_printf("    Calculating forward orbit...\n");
        x_for = x_var;
        v_for = v_var;
        LMCv_for = *LMCvelocity;
        LMCx_for = *LMCposition;

        // Get the initial acceleration
        DF_acc = dynamicalFriction_LMC(pot, &LMCx_for, &LMCv_for, LMCmass, LMCscale, LMCDynaFric, 0.0);
        LMC_acc = nbExtAcceleration(pot, &LMCx_for, 0.0);
        LMC_acc = mw_addv(&LMC_acc, &DF_acc);
        acc = nbExtAcceleration(pot, &x_for, 0.0);
        tmp = plummerAccel(&x_for, &LMCx_for, LMCmass, LMCscale);
        acc = mw_addv(&acc, &tmp);

        // Shift the body
        mw_acc = plummerAccel(&mw_x, &LMCx_for, LMCmass, LMCscale);
        LMC_acc = mw_subv(&LMC_acc, &mw_acc);
        acc = mw_subv(&acc, &mw_acc);

        for (t = 0; t < (ftime-tstop+(21*dt)); t += dt)
        {   
    	    exSteps = (int) mw_round_0(t/dt);
    	    if ((exSteps % 10 == 0)&&(t!=0)) { 
    	        ArrayPlaceholder[k] = mw_negv(&mw_acc);
                k++;
    	    }

            // Update the velocities and positions
            dv2_for.x = mw_mul_s(&acc.x, dt_half);
            dv2_for.y = mw_mul_s(&acc.y, dt_half);
            dv2_for.z = mw_mul_s(&acc.z, dt_half);
            v_for = mw_addv(&v_for, &dv2_for);

            dx_for.x = mw_mul_s(&v_for.x, dt);
            dx_for.y = mw_mul_s(&v_for.y, dt);
            dx_for.z = mw_mul_s(&v_for.z, dt);
            x_for = mw_addv(&x_for, &dx);

            dLMCv2_for.x = mw_mul_s(&LMC_acc.x, dt_half);
            dLMCv2_for.y = mw_mul_s(&LMC_acc.y, dt_half);
            dLMCv2_for.z = mw_mul_s(&LMC_acc.z, dt_half);
            LMCv_for = mw_addv(&LMCv_for, &dLMCv2_for);

            dLMCx_for.x = mw_mul_s(&LMCv_for.x, dt);
            dLMCx_for.y = mw_mul_s(&LMCv_for.y, dt);
            dLMCx_for.z = mw_mul_s(&LMCv_for.z, dt);
            LMCx_for = mw_addv(&LMCx_for, &dLMCx_for);
        
            // Compute the new acceleration
            DF_acc = dynamicalFriction_LMC(pot, &LMCx_for, &LMCv_for, LMCmass, LMCscale, LMCDynaFric, t);
            LMC_acc = nbExtAcceleration(pot, &LMCx_for, t);
            LMC_acc = mw_addv(&LMC_acc, &DF_acc);
            acc = nbExtAcceleration(pot, &x_for, t);
            tmp = plummerAccel(&x_for, &LMCx_for, LMCmass, LMCscale);
    	    acc = mw_addv(&acc, &tmp);

    	    // Shift the body
            mw_acc = plummerAccel(&mw_x, &LMCx_for, LMCmass, LMCscale);
            LMC_acc = mw_subv(&LMC_acc, &mw_acc);
            acc = mw_subv(&acc, &mw_acc);
        
            dv2_for.x = mw_mul_s(&acc.x, dt_half);
            dv2_for.y = mw_mul_s(&acc.y, dt_half);
            dv2_for.z = mw_mul_s(&acc.z, dt_half);
            v_for = mw_addv(&v_for, &dv2_for);

            dLMCv2_for.x = mw_mul_s(&LMC_acc.x, dt_half);
            dLMCv2_for.y = mw_mul_s(&LMC_acc.y, dt_half);
            dLMCv2_for.z = mw_mul_s(&LMC_acc.z, dt_half);
            LMCv_for = mw_addv(&LMCv_for, &dLMCv2_for);

        }
        ArrayPlaceholder[k] = mw_negv(&mw_acc); //set the last index after the loop ends

        //Fill forward orbit of shift array
        for(j = 0; j < k+1; j++) {
            tmp = ArrayPlaceholder[j];
            shiftByLMC[i+1+j] = tmp;
        }
    }

    unsigned int size = i + k + 2;
//    for(int l = 0; l < size; l++) {
//        mw_printf("shiftByLMC[%u] = [%.15f, %.15f, %.15f]\n", l, showRealValue(&shiftByLMC[l].x), showRealValue(&shiftByLMC[l].y), showRealValue(&shiftByLMC[l].z));
//    }

    //Free placeholder arrays
    mwFreeA(ArrayPlaceholder);
    nShiftLMC = size;

    /* Report the final values (don't forget to reverse the velocities) */
    v = mw_negv(&v);
    LMCv = mw_negv(&LMCv);

  #if AUTODIFF
    /* Insert reverse time derivative information into positions and velocities */
    getTimeDerivativeInfo(pot, &x, &v, &LMCx, &LMCv, LMCmass, LMCscale, LMCDynaFric, negT);
  #endif

    *finalPos = x;
    *finalVel = v;
    *LMCfinalPos = LMCx;
    *LMCfinalVel = LMCv;

    mw_printf("Initial Dwarf Position = [%.15f, %.15f, %.15f]\n", showRealValue(&x.x), showRealValue(&x.y), showRealValue(&x.z));
    mw_printf("Initial Dwarf Velocity = [%.15f, %.15f, %.15f]\n", showRealValue(&v.x), showRealValue(&v.y), showRealValue(&v.z));
    mw_printf("Initial LMC Position = [%.15f, %.15f, %.15f]\n", showRealValue(&LMCx.x), showRealValue(&LMCx.y), showRealValue(&LMCx.z));
    mw_printf("Initial LMC Velocity = [%.15f, %.15f, %.15f]\n", showRealValue(&LMCv.x), showRealValue(&LMCv.y), showRealValue(&LMCv.z));

    //Store LMC position and velocity
    LMCpos = LMCx;
    LMCvel = LMCv;
    mw_printf("Reverse Orbit Calculated!\n");
}

void getLMCArray(mwvector ** shiftArrayPtr, size_t * shiftSizePtr) {
    //Allows access to shift array
    *shiftArrayPtr = shiftByLMC;
    *shiftSizePtr = nShiftLMC;
}

void getLMCPosVel(mwvector * LMCposPtr, mwvector * LMCvelPtr) {
    //Allows access to LMC position and velocity
    *LMCposPtr = LMCpos;
    *LMCvelPtr = LMCvel;
}


void nbPrintReverseOrbit(mwvector* finalPos,
                         mwvector* finalVel,
                         const Potential* pot,
                         mwvector* pos,
                         mwvector* vel,
                         real_0 tstop,
                         real_0 tstopforward,
                         real_0 dt)
{
    mwvector acc, v, x;
    mwvector v_for, x_for, tmp;
    mwvector lbr;
    real_0 t;
    real_0 dt_half = dt / 2.0;

    // Set the initial conditions
    x = *pos;
    v = mw_negv(vel);
    x_for = *pos;
    v_for = *vel;

    // Get the initial acceleration
    acc = nbExtAcceleration(pot, &x, 0);

    FILE * fp;
    fp = fopen("reverse_orbit.out", "w");
    // Loop through time
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v = mw_addv(&v, &tmp);

        tmp.x = mw_mul_s(&X(&v), dt);
        tmp.y = mw_mul_s(&Y(&v), dt);
        tmp.z = mw_mul_s(&Z(&v), dt);
        x = mw_addv(&x, &tmp);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, &x, t);
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v = mw_addv(&v, &tmp);
        
        lbr = cartesianToLbr(&x, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(&X(&x)), showRealValue(&Y(&x)), showRealValue(&Z(&x)), showRealValue(&X(&lbr)), showRealValue(&Y(&lbr)), showRealValue(&Z(&lbr)), showRealValue(&X(&v)), showRealValue(&Y(&v)), showRealValue(&Z(&v)));
    }

    fclose(fp);
    fp = fopen("forward_orbit.out", "w");
    acc = nbExtAcceleration(pot, &x_for, 0);
    for (t = 0; t <= tstopforward; t += dt)
    {
        // Update the velocities and positions
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v_for = mw_addv(&v_for, &tmp);

        tmp.x = mw_mul_s(&X(&v_for), dt);
        tmp.y = mw_mul_s(&Y(&v_for), dt);
        tmp.z = mw_mul_s(&Z(&v_for), dt);
        x_for = mw_addv(&x_for, &tmp);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, &x_for, t);
        tmp.x = mw_mul_s(&X(&acc), dt_half);
        tmp.y = mw_mul_s(&Y(&acc), dt_half);
        tmp.z = mw_mul_s(&Z(&acc), dt_half);
        v_for = mw_addv(&v_for, &tmp);
        
        lbr = cartesianToLbr(&x_for, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(&X(&x_for)), showRealValue(&Y(&x_for)), showRealValue(&Z(&x_for)), showRealValue(&X(&lbr)), showRealValue(&Y(&lbr)), showRealValue(&Z(&lbr)), showRealValue(&X(&v_for)), showRealValue(&Y(&v_for)), showRealValue(&Z(&v_for)));
    }
    fclose(fp);
    
    /* Report the final values (don't forget to reverse the velocities) */
    v = mw_negv(&v);

    *finalPos = x;
    *finalVel = v;
}

void nbPrintReverseOrbit_LMC(mwvector* finalPos,
                         mwvector* finalVel,
                         mwvector* LMCfinalPos,
                         mwvector* LMCfinalVel,
                         const Potential* pot,
                         mwvector* pos,
                         mwvector* vel,
                         mwvector* LMCposition,
                         mwvector* LMCvelocity,
                         mwbool LMCDynaFric,
                         real_0 tstop,
                         real_0 tstopforward,
                         real_0 dt,
                         real* LMCmass,
                         real* LMCscale)
{
    mwvector acc, mw_acc, LMC_acc, DF_acc, tmp, dv2, dx, dLMCv2, dLMCx;
    mwvector v, x, LMCv, LMCx;
    mwvector v_for, x_for, LMCv_for, LMCx_for;
    mwvector mw_x = mw_vec(ZERO_REAL, ZERO_REAL, ZERO_REAL);
    mwvector lbr;
    real_0 t;
    real_0 dt_half = dt/2.0;

    // Set the initial conditions
    x = *pos;
    v = *vel;
    LMCx = *LMCposition;
    LMCv = *LMCvelocity;
    x_for = *pos;
    v_for = *vel;
    LMCx_for = *LMCposition;
    LMCv_for = *LMCvelocity;
    
    v = mw_negv(&v);
    LMCv = mw_negv(&LMCv);

    // Get the initial acceleration (reverse)
    mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
    LMC_acc = nbExtAcceleration(pot, &LMCx, 0);
    if (LMCDynaFric) {
        DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, TRUE, 0);
        DF_acc = mw_negv(&DF_acc); /* Inverting drag force for reverse orbit */
        LMC_acc = mw_addv(&LMC_acc, &DF_acc);
     }
    acc = nbExtAcceleration(pot, &x, 0);
    tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
    acc = mw_addv(&acc, &tmp);
    LMC_acc = mw_subv(&LMC_acc,& mw_acc);
    acc = mw_subv(&acc, &mw_acc);

    FILE * fp;
    fp = fopen("reverse_orbit.out", "w");

    // Loop through time (reverse)
    for (t = 0; t >= tstop*(-1); t -= dt)
    {
        // Update the velocities and positions
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v = mw_addv(&v, &dv2);

        dx.x = mw_mul_s(&v.x, dt);
        dx.y = mw_mul_s(&v.y, dt);
        dx.z = mw_mul_s(&v.z, dt);
        x = mw_addv(&x, &dx);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv = mw_addv(&LMCv, &dLMCv2);

        dLMCx.x = mw_mul_s(&LMCv.x, dt);
        dLMCx.y = mw_mul_s(&LMCv.y, dt);
        dLMCx.z = mw_mul_s(&LMCv.z, dt);
        LMCx = mw_addv(&LMCx, &dLMCx);
        
        // Compute the new acceleration
        mw_acc = plummerAccel(&mw_x, &LMCx, LMCmass, LMCscale);
        LMC_acc = nbExtAcceleration(pot, &LMCx, t);
        if (LMCDynaFric) {
            DF_acc = dynamicalFriction_LMC(pot, &LMCx, &LMCv, LMCmass, LMCscale, TRUE, t);
            LMC_acc = mw_subv(&LMC_acc, &DF_acc); /* Inverting drag force for reverse orbit */
        }
        acc = nbExtAcceleration(pot, &x, t);
        tmp = plummerAccel(&x, &LMCx, LMCmass, LMCscale);
    	acc = mw_addv(&acc, &tmp);

    	// Shift the body
    	mw_acc = mw_negv(&mw_acc);
        LMC_acc = mw_addv(&LMC_acc, &mw_acc);
        acc = mw_addv(&acc, &mw_acc);
        
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v = mw_addv(&v, &dv2);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv = mw_addv(&LMCv, &dLMCv2);

        
        lbr = cartesianToLbr(&x, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(&X(&x)), showRealValue(&Y(&x)), showRealValue(&Z(&x)), showRealValue(&X(&lbr)), showRealValue(&Y(&lbr)), showRealValue(&Z(&lbr)), showRealValue(&X(&v)), showRealValue(&Y(&v)), showRealValue(&Z(&v)));
    }
    fclose(fp);

    fp = fopen("forward_orbit.out", "w");

    //Get the initial acceleration (forward)
    mw_acc = plummerAccel(&mw_x, &LMCx_for, LMCmass, LMCscale);
    DF_acc = dynamicalFriction_LMC(pot, &LMCx_for, &LMCv_for, LMCmass, LMCscale, LMCDynaFric, 0);
    LMC_acc = nbExtAcceleration(pot, &LMCx_for, 0);
    LMC_acc = mw_addv(&LMC_acc, &DF_acc);
    acc = nbExtAcceleration(pot, &x_for, 0);
    tmp = plummerAccel(&x_for, &LMCx_for, LMCmass, LMCscale);
    acc = mw_addv(&acc, &tmp);
    LMC_acc = mw_subv(&LMC_acc, &mw_acc);
    acc = mw_subv(&acc, &mw_acc);

    //Loop through time (forward)
    for (t = 0; t <= tstopforward; t += dt)
    {
        // Update the velocities and positions
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v_for = mw_addv(&v_for, &dv2);

        dx.x = mw_mul_s(&v_for.x, dt);
        dx.y = mw_mul_s(&v_for.y, dt);
        dx.z = mw_mul_s(&v_for.z, dt);
        x_for = mw_addv(&x_for, &dx);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv_for = mw_addv(&LMCv_for, &dLMCv2);

        dLMCx.x = mw_mul_s(&LMCv_for.x, dt);
        dLMCx.y = mw_mul_s(&LMCv_for.y, dt);
        dLMCx.z = mw_mul_s(&LMCv_for.z, dt);
        LMCx_for = mw_addv(&LMCx_for, &dLMCx);
    
        // Compute the new acceleration
        mw_acc = plummerAccel(&mw_x, &LMCx_for, LMCmass, LMCscale);
        DF_acc = dynamicalFriction_LMC(pot, &LMCx_for, &LMCv_for, LMCmass, LMCscale, LMCDynaFric, t);
        LMC_acc = nbExtAcceleration(pot, &LMCx_for, t);
        LMC_acc = mw_addv(&LMC_acc, &DF_acc);
        acc = nbExtAcceleration(pot, &x_for, t);
        tmp = plummerAccel(&x_for, &LMCx_for, LMCmass, LMCscale);
        acc = mw_addv(&acc, &tmp);
        // Shift the body
        mw_acc = mw_negv(&mw_acc);
        LMC_acc = mw_addv(&LMC_acc, &mw_acc);
        acc = mw_addv(&acc, &mw_acc);
        
        dv2.x = mw_mul_s(&acc.x, dt_half);
        dv2.y = mw_mul_s(&acc.y, dt_half);
        dv2.z = mw_mul_s(&acc.z, dt_half);
        v_for = mw_addv(&v_for, &dv2);

        dLMCv2.x = mw_mul_s(&LMC_acc.x, dt_half);
        dLMCv2.y = mw_mul_s(&LMC_acc.y, dt_half);
        dLMCv2.z = mw_mul_s(&LMC_acc.z, dt_half);
        LMCv_for = mw_addv(&LMCv_for, &dLMCv2);
        
        lbr = cartesianToLbr(&x_for, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", showRealValue(&X(&x_for)), showRealValue(&Y(&x_for)), showRealValue(&Z(&x_for)), showRealValue(&X(&lbr)), showRealValue(&Y(&lbr)), showRealValue(&Z(&lbr)), showRealValue(&X(&v_for)), showRealValue(&Y(&v_for)), showRealValue(&Z(&v_for)));
    }
    fclose(fp);
    
    /* Report the final values (don't forget to reverse the velocities) */
    *finalPos = x;
    *finalVel = mw_negv(&v);
    *LMCfinalPos = LMCx;
    *LMCfinalVel = mw_negv(&LMCv);
}
