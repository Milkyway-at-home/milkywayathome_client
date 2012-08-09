#include <cfloat>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <unistd.h>
#include <sys/wait.h> 

extern "C" {
#include "lauxlib.h"
}

#include "tao/evolutionary_algorithms/particle_swarm.hxx"
#include "tao/evolutionary_algorithms/differential_evolution.hxx"

#include "tao/synchronous_algorithms/parameter_sweep.hxx"
#include "tao/synchronous_algorithms/synchronous_gradient_descent.hxx"
#include "tao/synchronous_algorithms/synchronous_newton_method.hxx"

//from undvc_common
#include "tao/undvc_common/arguments.hxx"
#include "tao/undvc_common/vector_io.hxx"

using namespace std;

//Run Independent Constants
char* PATH_TO_SEPARATION;
char* PATH_TO_STARS;  //Needs -s in front of path name without any spaces
int number_streams;
int number_cuts;
int wedge;
vector<double> area;

int lua_get_area(lua_State* L) 
{
    int c = 0;
    lua_getglobal(L, "area");
    lua_pushnil(L); 
    while(lua_next(L, -2))
    {
        lua_getfield(L, -1, "r_min");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "r_max");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "r_steps");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "mu_min");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "mu_max");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "mu_steps");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "nu_min");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "nu_max");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "nu_steps");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in area" << endl;
            return 0;
        }
        area.push_back(lua_tonumber(L, -1));
        lua_pop(L, 2);
        c++;
    }
    lua_pop(L, 2);

    number_cuts = c;

    return 1;
}

//  Assumes table containing background is at top of the stack
int lua_get_background(lua_State* L, vector<double> & v) 
{
    lua_getfield(L, -1, "background");
    lua_getfield(L, -1, "q");
    if (!lua_isnumber(L, -1))
    {
        cerr << "invalid component in background" << endl;
        return 0;
    }
    v.push_back(lua_tonumber(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, -1, "r0");
    if (!lua_isnumber(L, -1))
    {
        cerr << "invalid component in background" << endl;
        return 0;
    }
    v.push_back(lua_tonumber(L, -1));
    lua_pop(L, 2);

    return 1;
}

//  Assumes table containing background is at top of the stack
int lua_get_stream(lua_State* L, vector<double> & v) 
{
    int s = 0;
    lua_getfield(L, -1, "streams");
    lua_pushnil(L);
    while(lua_next(L, -2))
    {
        lua_getfield(L, -1, "epsilon");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in stream" << endl;
            return 0;
        }
        v.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "mu");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in stream" << endl;
            return 0;
        }
        v.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "r");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in stream" << endl;
            return 0;
        }
        v.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "theta");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in stream" << endl;
            return 0;
        }
        v.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "phi");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in stream" << endl;
            return 0;
        }
        v.push_back(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, -1, "sigma");
        if (!lua_isnumber(L, -1))
        {
            cerr << "invalid component in stream" << endl;
            return 0;
        }
        v.push_back(lua_tonumber(L, -1));
        lua_pop(L, 2);
        s++;
    }
    lua_pop(L, 2);

    number_streams = s;

    return 1;
}

double objective_function(const vector<double> &A) 
{
    ofstream input ( "input.lua" );

    input << "\nwedge = " << wedge;

    input << setprecision(16);

    input << "\n\nbackground = {\n\tq = " << A[0] << ",\n\tr0 = " << A[1] << "\n}\n\nstreams = {\n";

    for(int s=0; s<number_streams; s++)
        input << "\t{\n\t\tepsilon = " << A[2+s*6] << ",\n\t\tmu = " << A[3+s*6] << ",\n\t\tr = " << A[4+s*6] 
              << ",\n\t\ttheta = " << A[5+s*6] << ",\n\t\tphi = " << A[6+s*6] << ",\n\t\tsigma = " << A[7+s*6] << "\n\t},\n";

    input << "}\n\narea = {\n";

    for(int c=0; c<number_cuts; c++)
        input << "\t{\n\t\tr_min = " << area[c*9] << ",\n\t\tr_max = " << area[1+c*9] << ",\n\t\tr_steps = " << area[2+c*9] 
               << ",\n\n\t\tmu_min = " << area[3+c*9] << ",\n\t\tmu_max = " << area[4+c*9] << ",\n\t\tmu_steps = "
               << area[5+c*9] << ",\n\n\t\tnu_min = " << area[6+c*9] << ",\n\t\tnu_max = " << area[7+c*9] 
               << ",\n\t\tnu_steps = " << area[8+c*9] << "\n\t},\n";

    input << "}" << endl;

    pid_t pID = fork();
    if (pID == 0)
    {
        // Code only executed by child process
        execl ( PATH_TO_SEPARATION, "milkyway_separation", "-c", "-i", "-t","-ainput.lua", PATH_TO_STARS, "-f", NULL);
        _exit(0);
    }
    else if (pID < 0)
    {
        cerr << "Failed to fork" << endl;
        exit(1);
    }

    int childExitStatus;
    waitpid( pID, &childExitStatus, 0);

    if( !WIFEXITED(childExitStatus) )
    {
        cerr << "waitpid() exited with an error: Status= " << WEXITSTATUS(childExitStatus) << endl;
        return 0;
    }
    else if( WIFSIGNALED(childExitStatus) )
    {
        cerr << "waitpid() exited due to a signal: " << WTERMSIG(childExitStatus) << endl;
        return 0;
    }

    ifstream results ( "results.txt" );

    double likelihood;

    results >> likelihood;

    return likelihood;
}

int main(int number_arguments, char **argv) 
{

    vector<string> arguments(argv, argv + number_arguments);

    string stemp;
    get_argument(arguments, "--separation", true, stemp);
    PATH_TO_SEPARATION = (char*)stemp.c_str();

    get_argument(arguments, "--params", true, stemp);
    char *search_params = (char*)stemp.c_str();

    lua_State *L = lua_open();
    if (!L)
    {
        cerr << "Failed to get Lua state" << endl;
        return 0;
    }

    luaL_dofile(L, search_params);

    double temp;
    lua_getglobal(L, "wedge");
    temp = lua_tonumber(L, -1);
    lua_pop(L, 1);
    wedge = (int) temp;

    lua_get_area(L);    //Also sets number of cuts

    vector<double> min_bound;
    lua_getglobal(L, "min");
    lua_get_background(L, min_bound);
    lua_get_stream(L, min_bound);   //Also sets number of streams
    lua_pop(L, 1);

    vector<double> max_bound;
    lua_getglobal(L, "max");
    lua_get_background(L, max_bound);
    lua_get_stream(L, max_bound);
    lua_pop(L, 1);

    // Not sure why this has to be down here, but opening the Lua file somehow makes the variable unreadable if it is done first
    get_argument(arguments, "--stars", true, stemp);
    stemp = "-s" + stemp;
    PATH_TO_STARS = (char*)stemp.c_str();

    string search_type;
    get_argument(arguments, "--search_type", true, search_type);

    if (search_type.compare("ps") == 0) {
        ParticleSwarm ps(min_bound, max_bound, arguments);
        ps.iterate(objective_function);

    } else if (search_type.compare("de") == 0) {
        DifferentialEvolution de(min_bound, max_bound, arguments);
        de.iterate(objective_function);

    } else if (search_type.compare("sweep") == 0) {
        vector<double> step_size;
        lua_getglobal(L, "step");
        lua_get_background(L, step_size);
        lua_get_stream(L, step_size);
        lua_pop(L, 1);

        parameter_sweep(min_bound, max_bound, step_size, objective_function);

    } else if (search_type.compare("snm") == 0 || search_type.compare("gd") == 0 || search_type.compare("cgd") == 0) {
        vector<double> starting_point;
        lua_getglobal(L, "start");
        lua_get_background(L, starting_point);
        lua_get_stream(L, starting_point);
        lua_pop(L, 1);

        double range;
        if (get_argument(arguments, "--rand", false, range)) 
        {
            cout << "Randomizing search paramaters by +- " << range << endl;
            srand48((long)time(0));
            for(int i=0; i<starting_point.size(); i++)
            {
                starting_point[i] = starting_point[i]*(1.0-range) + drand48()*(starting_point[i]*(1.0+range)-starting_point[i]*(1.0-range));
            }
        }

        vector<double> step_size;
        lua_getglobal(L, "step");
        lua_get_background(L, step_size);
        lua_get_stream(L, step_size);
        lua_pop(L, 1);

        if (search_type.compare("snm") == 0) {
            synchronous_newton_method(arguments, objective_function, starting_point, step_size);
        } else if (search_type.compare("gd") == 0) {
            synchronous_gradient_descent(arguments, objective_function, starting_point, step_size);
        } else if (search_type.compare("cgd") == 0) {
            synchronous_conjugate_gradient_descent(arguments, objective_function, starting_point, step_size);
        }

    } else {
        fprintf(stderr, "Improperly specified search type: '%s'\n", search_type.c_str());
        fprintf(stderr, "Possibilities are:\n");
        fprintf(stderr, "    de     -       differential evolution\n");
        fprintf(stderr, "    ps     -       particle swarm optimization\n");
        fprintf(stderr, "    snm    -       synchronous newton method\n");
        fprintf(stderr, "    gd     -       gradient descent\n");
        fprintf(stderr, "    cgd    -       conjugate gradient descent\n");
        exit(0);
    }

    //lua_close(L);

    remove ( "results.txt" );
    remove ( "input.lua" );
    remove ( "separation_checkpoint");

    return 0;
}
