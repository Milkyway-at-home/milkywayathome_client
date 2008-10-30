/********
	*	FGDO includes
 ********/
#include "search_manager.h"

int main(int argc, char** argv) {
	char working_directory[512], search_name[512];
	/********
		*	Arguments
		*		-d <working_directory>
		*		-s <search_name>
		*		-wus <wus_to_generate>
	 ********/

	search_manager__register_search("nm", start_newton_method);
	search_manager__register_search("gs", start_genetic_search);
	search_manager__register_search("de", start_differential_evolution);
	search_manager__register_search("pso", start_particle_swarm);

	init_search_manager(working_directory, wu_init_func, wu_generate_func, wu_process_func);
	load_asynchronous_search(search_name);
	generate_workunits(search_name, wus_to_generate);
}
