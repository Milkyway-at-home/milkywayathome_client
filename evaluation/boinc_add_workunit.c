#include <stdio.h>
#include "boinc_add_workunit.h"
#include "../util/settings.h"


int fwrite_workunit_info(FILE *file, WORKUNIT_INFO *info) {
	int i;
	fprintf(file, "number_parameters: %d\n", info->number_parameters);
	fprintf(file, "rsc_fpops_est: %lf, rsc_fpops_bound: %lf, rsc_memory_bound: %lf, rsc_disk_bound: %lf\n",
		info->rsc_fpops_est, info->rsc_fpops_bound, info->rsc_memory_bound, info->rsc_disk_bound);
	fprintf(file, "delay_bound: %d, min_quorum: %d, target_nresults: %d, max_error_results: %d, max_total_results: %d, max_success_results: %d\n",
		info->delay_bound, info->min_quorum, info->target_nresults, info->max_error_results, info->max_total_results, info->max_success_results);
	fprintf(file, "credit_str: %s\n", info->credit_str);
	fprintf(file, "result_xml_path: %s\n", info->result_xml_path);
	fprintf(file, "template_filename: %s\n", info->template_filename);
	fprintf(file, "required_files[%d]:\n", info->number_required_files);
	for (i = 0; i < info->number_required_files; i++) {
		fprintf(file, "%s\n", info->required_files[i]);
	}
	return 0;
}

int write_workunit_info(char *filename, WORKUNIT_INFO *info) {
	int retval;
	FILE* data_file = fopen(filename, "w");
	if (!data_file) {
		fprintf(stderr, "Couldn't find output file [%s] to write workunit info.\n", filename);
		return 1;
	}
	retval = fwrite_workunit_info(data_file, info);
	if (retval) {
		fprintf(stderr, "Couldn't write to output file [%s] to write workunit info, error: %d\n", filename, retval);
		return retval;
	}
	fclose(data_file);
	return 0;
}

int fread_workunit_info(FILE *file, WORKUNIT_INFO **info) {
	int i;
	(*info) = (WORKUNIT_INFO*)malloc(sizeof(WORKUNIT_INFO));
	fscanf(file, "number_parameters: %d\n", &((*info)->number_parameters));
	fscanf(file, "rsc_fpops_est: %lf, rsc_fpops_bound: %lf, rsc_memory_bound: %lf, rsc_disk_bound: %lf\n",
		&((*info)->rsc_fpops_est), &((*info)->rsc_fpops_bound), &((*info)->rsc_memory_bound), &((*info)->rsc_disk_bound));
	fscanf(file, "delay_bound: %d, min_quorum: %d, target_nresults: %d, max_error_results: %d, max_total_results: %d, max_success_results: %d\n",
		&((*info)->delay_bound), &((*info)->min_quorum), &((*info)->target_nresults), &((*info)->max_error_results), &((*info)->max_total_results), &((*info)->max_success_results));

	(*info)->credit_str = (char*)malloc(sizeof(char) * 1024);
	(*info)->result_xml_path = (char*)malloc(sizeof(char) * FILENAME_SIZE);
	fscanf(file, "credit_str: %s\n", (*info)->credit_str);
	fscanf(file, "result_xml_path: %s\n", (*info)->result_xml_path);

	(*info)->template_filename = (char*)malloc(sizeof(char) * FILENAME_SIZE);
	fscanf(file, "template_filename: %s\n", (*info)->template_filename);
	if (read_file_malloc((*info)->template_filename, (*info)->template_file)) {
		fprintf(stderr, "Could not read workunit result template: %s\n", (*info)->template_filename);
		exit(0);
	}

	fscanf(file, "required_files[%d]:\n", &((*info)->number_required_files));
	(*info)->required_files = (char**)malloc(sizeof(char*) * (*info)->number_required_files);
	for (i = 0; i < (*info)->number_required_files; i++) {
		(*info)->required_files[i] = (char*)malloc(sizeof(char) * FILENAME_SIZE);
		fscanf(file, "%s\n", (*info)->required_files[i]);
	}
	return 0;
}

int read_workunit_info(char *filename, WORKUNIT_INFO **info) {
	int retval;
	FILE* data_file = fopen(filename, "r");
	if (!data_file) {
		fprintf(stderr, "Couldn't find input file [%s] to read workunit info.\n", filename);
		return 1;
	}
	retval = fread_workunit_info(data_file, info);
	if (retval) {
		fprintf(stderr, "Couldn't read to input file [%s] to write workunit info, error: %d\n", filename, retval);
		return retval;
	}
	fclose(data_file);
	return 0;
}
