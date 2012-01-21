/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute.
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "cl_compile_flags.h"
#include "setup_cl.h"

#include <sstream>

/* Return flag for Nvidia compiler for maximum registers to use. */
static const char* getNvidiaRegCount(const DevInfo* di)
{
    const char* regCount32 = "-cl-nv-maxrregcount=32 ";
    const char* regDefault = "";

    if (mwComputeCapabilityIs(di, 1, 3)) /* 1.3 == GT200 */
    {
        /* 32 allows for greatest number of threads at a time */
        mw_printf("Found a compute capability 1.3 device. Using %s\n", regCount32);
        return regCount32;
    }

    /* Higher or other is Fermi or unknown, */
    return regDefault;
}

/* Get string of options to pass to the CL compiler. Result must be freed */
char* getCompilerFlags(const CLInfo* ci, const AstronomyParameters* ap, cl_bool useILKernel)
{
    const DevInfo* di = &ci->di;
    std::string flagStr;
    const char* str;
    std::stringstream flags(std::stringstream::out);
    flags.precision(15);

    if (DOUBLEPREC)
    {
        flags << "-D DOUBLEPREC=1 ";
    }
    else
    {
        flags << "-D DOUBLEPREC=0 -cl-single-precision-constant ";
    }

    /* Math options for CL compiler */
    flags << "-cl-mad-enable ";
    flags << "-cl-no-signed-zeros ";
    flags << "-cl-strict-aliasing ";
    flags << "-cl-finite-math-only ";

    /* Get constant definitions */
    flags << "-D FAST_H_PROB="    << ap->fast_h_prob    << " ";
    flags << "-D AUX_BG_PROFILE=" << ap->aux_bg_profile << " ";
    flags << "-D NSTREAM="        << ap->number_streams << " ";
    flags << "-D CONVOLVE="       << ap->convolve       << " ";
    flags << "-D R0="             << ap->r0             << " ";
    flags << "-D SUN_R0="         << ap->sun_r0         << " ";
    flags << "-D Q_INV_SQR="      << ap->q_inv_sqr      << " ";
    flags << "-D BG_A="           << ap->bg_a           << " ";
    flags << "-D BG_B="           << ap->bg_b           << " ";
    flags << "-D BG_C="           << ap->bg_c           << " ";

    flags << "-D ALPHA="          << ap->alpha          << " ";
    flags << "-D ALPHA_DELTA_3="  << ap->alpha_delta3   << " ";


    /* FIXME: Device vendor not necessarily the platform vendor */
    if (mwHasNvidiaCompilerFlags(di))
    {
        flags << "-cl-nv-verbose ";
        flags << getNvidiaRegCount(di);
    }

    if (useILKernel)
    {
        /* Options to only emit AMD IL in the binary */
        flags << "-fno-bin-exe -fno-bin-llvmir -fno-bin-source -fbin-amdil ";
    }

    flagStr = flags.str();
    str = flagStr.c_str();
    return str ? strdup(str) : NULL;
}

