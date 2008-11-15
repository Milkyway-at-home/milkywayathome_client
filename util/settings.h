#ifndef FGDO_SETTINGS_H
#define FGDO_SETTINGS_H

#define FILENAME_SIZE		2048
#define SEARCH_QUALIFIER_SIZE	64
#define SEARCH_NAME_SIZE	512
#define METADATA_SIZE		2048
#define BOINC_SEARCH_PATH "./boinc_testing/"

char* get_working_directory();
void set_working_directory(char* working_directory);

#endif
