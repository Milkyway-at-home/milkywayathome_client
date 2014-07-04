/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <stdint.h>

#include "undvc_common/parse_xml.hxx"

using std::cerr;
using std::endl;
using std::vector;
using std::string;

/**
 *  Need to pass in nbodies somehow
 */

void calculate_fpops(const vector<double> &parameters, double &rsc_fpops_est, double &rsc_fpops_bound, string workunit_extra_xml) {
    double simulation_time  = parameters[0];
    //    double orbit_time       = parameters[1];
    double radius_1         = parameters[2];
    double radius_2         = parameters[3];
    double mass_1           = parameters[4];
    double mass_2           = parameters[5];
    
    double max_radius;
    if (radius_2 > radius_1)    max_radius = radius_2;
    else                        max_radius = radius_1;
    
    uint64_t n_bodies = 0;
    try {
        n_bodies = parse_xml<uint64_t>(workunit_extra_xml, "n_bodies");
    } catch (string ex_msg) {
        cerr << "ERROR parsing workunit_extra_xml on " << __FILE__ << ":" << __LINE__ << endl;
        cerr << "\t" << ex_msg << endl;
        cerr << "xml is:" << endl;
        cerr << workunit_extra_xml << endl;
        exit(1);
    }
    
    
    double timestep = (0.1 * 0.1) * sqrt(M_PI * (4.0/3.0) * max_radius * max_radius * max_radius / (mass_1 + mass_2));
    double step_fpops = (6 + 3 + (7 * 5) + (2 * 10) + 20) * (n_bodies * n_bodies);
    double fpops = step_fpops * (simulation_time / timestep);
    
    //    cerr << "workunit_extra_xml: " << workunit_extra_xml << endl;
    //    cerr << "n_bodies: " << n_bodies << endl;
    //    cerr << "parameters: " << vector_to_string(parameters) << endl;
    //    cerr << "timestep: " << timestep << endl;
    //    cerr << "step_fpops: " << step_fpops << endl;
    //    cerr << "fpops: " << fpops << endl;
    
    rsc_fpops_est = fpops * 10;
    rsc_fpops_bound = fpops * 100000;
}
