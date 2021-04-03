#ifndef INC_IO_H
#define INC_IO_H
#include <mpirxx.h>
#include <string>
#include "cruncher.h"

void read_from_file(mpq_class* frac, std::string file_name, bool dec = true);
void output_cf_terms_list(const CFTerms&, std::string file_name, bool verbose = true);

#endif
