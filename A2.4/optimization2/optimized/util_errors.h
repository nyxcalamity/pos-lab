/**
 * Provides macos for error handling.
 */
#ifndef __UTIL_ERRORS_H__
#define __UTIL_ERRORS_H__

#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "posl_definitions.h"

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_err(M, ...) fprintf(stdout, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define log_warn(M, ...) fprintf(stdout, "[WARN] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define log_info(M, ...) fprintf(stdout, "[INFO] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#define log_dbg(M, ...) fprintf(stdout, "[DEBUG] (%s:%d): " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#define check_status(R, E, M, ...) if (E != 0) { log_err(M " on process #%d", ##__VA_ARGS__, R); return E; }
#define check_allocation(R, P, M, ...) if (P == NULL) { log_err(M " on process #%d", ##__VA_ARGS__, R); return POSL_ERROR; }
#define check_alloc(P, M, ...) if (P == NULL) { log_err(M, ##__VA_ARGS__); return POSL_ERROR; }

#endif