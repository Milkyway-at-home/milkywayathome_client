/********
        *       Includes for BOINC
********/
#ifdef _WIN32
	#include "boinc_win.h"
	#include "str_util.h"
#else
	#include "config.h"
#endif

#ifndef _WIN32
	#include <cstdio>
	#include <cctype>
	#include <ctime>
	#include <cstring>
	#include <cstdlib>
	#include <csignal>
	#include <unistd.h>
#endif

/* #define BOINC_APP_GRAPHICS */

#ifdef BOINC_APP_GRAPHICS
	#include "graphics_api.h"
	#include "graphics_lib.h"
#endif

#include "diagnostics.h"
#include "util.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"
#include "stdio.h"

using std::string;

#define OUTPUT_FILENAME "out"
#define PARAMETER_FILENAME "parameters.txt"
#define STAR_POINTS_FILENAME "stars.txt"

/********
	*	Includes for astronomy
 ********/
#include "boinc_parameters.h"
#include "boinc_star_points.h"
#include "evaluation.h"

#ifdef _WIN32
	void AppInvalidParameterHandler(const wchar_t* expression, const wchar_t* function, const wchar_t* file, unsigned int line, uintptr_t pReserved ) {
		fprintf(stderr, "Invalid parameter detected in function %s. File: %s Line: %d\n", function, file, line);
		fprintf(stderr, "Expression: %s\n", expression);
		// Cause a Debug Breakpoint.
		DebugBreak();
	}
#endif


void worker() {
	int i;
	/********
		*	READ THE ASTRONOMY PARAMETERS
	 ********/
	ASTRONOMY_PARAMETERS *ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	int retval = boinc_read_astronomy_parameters(PARAMETER_FILENAME, ap);
	if (retval) {
		fprintf(stderr, "APP: error reading astronomy parameters: %d\n", retval);
		boinc_finish(1);
	}

	/********
		*	READ THE STAR POINTS
	 ********/
	STAR_POINTS *sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
	retval = boinc_read_star_points(STAR_POINTS_FILENAME, sp);
	if (retval) {
		fprintf(stderr, "APP: error reading star points: %d\n", retval);
		boinc_finish(1);
	}

	/********
		*	INITIALIZE THE EVALUATION STATE
	 ********/
	EVALUATION_STATE *es = (EVALUATION_STATE*)malloc(sizeof(EVALUATION_STATE));
	initialize_state(es, ap->number_streams);

	/********
		*	READ THE SEARCH PARAMETERS
	 ********/



	/********
		*	CALCULATE THE INTEGRALS
	 ********/
	retval = calculate_integrals(ap, es, sp);
	if (retval) {
		fprintf(stderr, "APP: error calculating integrals: %d\n", retval);
		boinc_finish(retval);
	}

	/********
		*	CALCULATE THE LIKELIHOOD
	 ********/
	retval = calculate_likelihood(ap, es, sp);
	if (retval) {
		fprintf(stderr, "APP: error calculating likelihood: %d\n", retval);
		boinc_finish(retval);
	}
	double* double_parameters = get_as_double_array(ap);

	/********
		*	RESOLVE THE OUTPUT FILE & WRITE THE RESULT
	 ********/
	char output_path[512];
	retval = boinc_resolve_filename(OUTPUT_FILENAME, output_path, sizeof(output_path));

	MFILE out;
	retval = out.open(output_path, "w");
	if (retval) {
		fprintf(stderr, "APP: error opening output file\n");
		boinc_finish(retval);
	}

	es->prob_sum /= (sp->number_stars - es->bad_jacobians);

	out.printf("likelihood: %lf\n", es->prob_sum);
	out.printf("parameters [%d]:", ap->number_parameters);
	for (i = 0; i < ap->number_parameters; i++) {
		out.printf(" %lf", double_parameters[i]);
	}
	out.printf("\nvolume size: %d\n", (ap->r_steps * ap->mu_steps * ap->nu_steps) );
	out.printf("info: %s\n", ap->parameter_info);

	if (out.flush()) {
		fprintf(stderr, "APP: astronomy flush failed.\n");
		boinc_finish(1);
	}
	out.close();

	free(double_parameters);
	free_state(es);
	free(es);
	free_parameters(ap);
	free(ap);
	free_star_points(sp);
	free(sp);

#ifdef _WIN32
	_set_printf_count_output( 1 );
	_get_printf_count_output();
#endif

	boinc_finish(0);
}

int main(int argc, char **argv){
        int retval = 0;

#ifdef BOINC_APP_GRAPHICS
#if defined(_WIN32) || defined(__APPLE__)
	        retval = boinc_init_graphics(worker);
#else
                retval = boinc_init_graphics_lib(worker, argv[0]);
#endif
                if (retval) exit(retval);
#endif

	retval = boinc_init();
        if (retval) exit(retval);
        worker();
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode){
        LPSTR command_line;
        char* argv[100];
	int argc;

	command_line = GetCommandLine();
	argc = parse_command_line( command_line, argv );
	return main(argc, argv);
}
#endif

const char *BOINC_RCSID_33ac47a071 = "$Id: boinc_astronomy.C,v 1.3 2008/09/03 21:33:27 deselt Exp $";
