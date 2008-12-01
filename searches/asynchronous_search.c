#include "asynchronous_search.h"

const char *AS_GEN_STR[] = { "success", "search completed", "ERROR" };
const char *AS_INSERT_STR[] = { "success", "search completed", "fitness is NAN", "fitness is invalid", "parameters contain NAN", "parameters out of bounds" };
const char *AS_CP_STR[] = { "success", "search completed", "ERROR" };

char AS_MSG[1024] = "";
