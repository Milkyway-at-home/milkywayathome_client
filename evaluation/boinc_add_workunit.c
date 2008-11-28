#include <stdio.h>
#include "boinc_add_workunit.h"


void fwrite_workunit_info(FILE *file, WORKUNIT_INFO *info) {
	fprintf(file, "rsc_fpops_est: %lf, rsc_fpops_bound: %lf, rsc_memory_bound: %lf, rsc_disk_bound: %lf\n",
		info->rsc_fpops_est, info->rsc_fpops_bound, info->rsc_memory_bound, info->rsc_disk_bound);
	fprintf(file, "delay_bound: %d, min_quorum: %d, target_nresults: %d, max_error_results: %d, max_total_results: %d, max_success_results: %d\n",
		info->delay_bound, info->min_quorum, info->target_nresults, info->max_error_results, info->max_total_results, info->max_success_results);
	fprintf(file, "credit_str: %s, result_xml_path: %s\n", info->credit_str, info->result_xml_path);
	fprintf(file, "required_files[0]: %s, required_files[1]: %s\n", info->required_files[0], info->required_files[1]);
}
