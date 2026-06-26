//This test is designed to test the softening length array calculations for various different component profiles
//If a new dwarf profile is added, please add a corresponding test case below
//This test will also run a two component Nbody to make sure the array is correctly passed 
//(single component is tested in model tests)
//Like with the model tests, the expected result may need to be updated if changes are made to the simulation
#include "test_env_util.h"
#include "nbody_lua_models.h"
#include "nbody_potential_types.h"

int check_result(real* eps_array, real eps[3]) {
    for (int i = 0; i < 3; i++) {
        if (eps_array[i] != eps[i]) {
            return 0;
        }
    }
    return 1;
}

int main() {

    //Check softening lengths for various models. These are the same params as used in mixeddwarf test

    real eps_l;
    real eps_cross;
    real eps_d;
    real eps[3];
    real* eps_array;
    unsigned int lm_nbody = 10000;
    unsigned int nbody = 40000;
    Dwarf* light_comp = mwMalloc(sizeof(Dwarf));
    *light_comp = (Dwarf)EMPTY_DWARF;
    Dwarf* dark_comp = mwMalloc(sizeof(Dwarf));
    *dark_comp = (Dwarf)EMPTY_DWARF;
    int failed = 0;

    light_comp->mass = 12;
    light_comp->scaleLength = 0.2;
    dark_comp->mass = 48;
    dark_comp->scaleLength = 0.8;

    //Plummer-Plummer
    eps_l = 0.00000209086716398582624413214625658241629935218952596187591552734375;
    eps_cross = 0.00000247278851242636896762684510509000546107927220873534679412841796875;
    eps_d = 0.00000685814317670222965375075377392199982296006055548787117004394531250;

    eps[0] = eps_l;
    eps[1] = eps_cross;
    eps[2] = eps_d;

    light_comp->type = 0;
    dark_comp->type = 0;

    eps_array = nbCalculateEps2_NEW(light_comp, dark_comp, lm_nbody, nbody); 

    if (!check_result(eps_array, eps)) {
        mw_printf("Test failed for Plummer-Plummer model\n");
        mw_printf("Expected: %.80f, %.80f, %.80f\n", eps[0], eps[1], eps[2]);
        mw_printf("Got: %.80f, %.80f, %.80f\n", eps_array[0], eps_array[1], eps_array[2]);
        failed = 1;
    }
    else {
        mw_printf("Plummer-Plummer model passed\n");
    }

    //Plummer-Hernquist
    eps_l = 0.000001847711702106992331048093518297559256780004943720996379852294921875;
    eps_cross = 0.0000051708933699171807162422477566199319198858574964106082916259765625;
    eps_d = 0.00000831266825054630985390551056735120027951779775321483612060546875;

    eps[0] = eps_l;
    eps[1] = eps_cross;
    eps[2] = eps_d;

    dark_comp->type = 2;

    eps_array = nbCalculateEps2_NEW(light_comp, dark_comp, lm_nbody, nbody); 

    if (!check_result(eps_array, eps)) {
        mw_printf("Test failed for Plummer-Hernquist model\n");
        mw_printf("Expected: %.80f, %.80f, %.80f\n", eps[0], eps[1], eps[2]);
        mw_printf("Got: %.80f, %.80f, %.80f\n", eps_array[0], eps_array[1], eps_array[2]);
        failed = 1;
    }
    else {
        mw_printf("Plummer-Hernquist model passed\n");
    }

    //Plummer-NFW
    eps_l = 0.00000204024281439374797800473219921979506352727185003459453582763671875;
    eps_cross = 0.000003070486122463703874789497927366710428032092750072479248046875;
    eps_d = 0.00000607115034132603119810915603959955433310824446380138397216796875;

    eps[0] = eps_l;
    eps[1] = eps_cross;
    eps[2] = eps_d;

    dark_comp->type = 1;

    eps_array = nbCalculateEps2_NEW(light_comp, dark_comp, lm_nbody, nbody); 

    if (!check_result(eps_array, eps)) {
        mw_printf("Test failed for Plummer-NFW model\n");
        mw_printf("Expected: %.80f, %.80f, %.80f\n", eps[0], eps[1], eps[2]);
        mw_printf("Got: %.80f, %.80f, %.80f\n", eps_array[0], eps_array[1], eps_array[2]);
        failed = 1;
    }
    else {
        mw_printf("Plummer-NFW model passed\n");
    }

    //Plummer-Cored
    eps_l = 0.000002281581053752841276613221033198186660229112021625041961669921875;
    eps_cross = 0.000008748023194028815212146266144799255926045589148998260498046875;
    eps_d = 0.000035177338064575728606396542996748166842735372483730316162109375;

    eps[0] = eps_l;
    eps[1] = eps_cross;
    eps[2] = eps_d;

    dark_comp->type = 4;
    dark_comp->r1 = 0.7;
    dark_comp->rc = 0.6;

    eps_array = nbCalculateEps2_NEW(light_comp, dark_comp, lm_nbody, nbody); 

    if (!check_result(eps_array, eps)) {
        mw_printf("Test failed for Plummer-cored model\n");
        mw_printf("Expected: %.80f, %.80f, %.80f\n", eps[0], eps[1], eps[2]);
        mw_printf("Got: %.80f, %.80f, %.80f\n", eps_array[0], eps_array[1], eps_array[2]);
        failed = 1;
    }
    else {
        mw_printf("Plummer-cored model passed\n");
    }

    //Plummer-Cutoff-Cored
    eps_l = 0.000002281581053752841276613221033198186660229112021625041961669921875;
    eps_cross = 0.000008748023194028815212146266144799255926045589148998260498046875;
    eps_d = 0.000035177338064575728606396542996748166842735372483730316162109375;

    eps[0] = eps_l;
    eps[1] = eps_cross;
    eps[2] = eps_d;

    dark_comp->rcut = 4.5;

    eps_array = nbCalculateEps2_NEW(light_comp, dark_comp, lm_nbody, nbody); 

    if (!check_result(eps_array, eps)) {
        mw_printf("Test failed for Plummer-cutoff-cored model\n");
        mw_printf("Expected: %.80f, %.80f, %.80f\n", eps[0], eps[1], eps[2]);
        mw_printf("Got: %.80f, %.80f, %.80f\n", eps_array[0], eps_array[1], eps_array[2]);
        failed = 1;
    }
    else {
        mw_printf("Plummer-cutoff-cored model passed\n");
    }

    //Plummer-Cutoff-NFW
    eps_l = 0.00000204024281439380557624514549164285170945731806568801403045654296875;
    eps_cross = 0.000003070486122463703874789497927366710428032092750072479248046875;
    eps_d = 0.00000607115034132603119810915603959955433310824446380138397216796875;

    eps[0] = eps_l;
    eps[1] = eps_cross;
    eps[2] = eps_d;

    dark_comp->type = 1;
    dark_comp->r1 = 0.0;
    dark_comp->rc = 0.0;

    eps_array = nbCalculateEps2_NEW(light_comp, dark_comp, lm_nbody, nbody); 

    if (!check_result(eps_array, eps)) {
        mw_printf("Test failed for Plummer-cutoff-NFW model\n");
        mw_printf("Expected: %.80f, %.80f, %.80f\n", eps[0], eps[1], eps[2]);
        mw_printf("Got: %.80f, %.80f, %.80f\n", eps_array[0], eps_array[1], eps_array[2]);
        failed = 1;
    }
    else {
        mw_printf("Plummer-Cutoff-NFW model passed\n");
    }

    return failed;
}

