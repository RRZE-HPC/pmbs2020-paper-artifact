#ifndef RACE_PARSE_H
#define RACE_PARSE_H

#include <vector>
#include <getopt.h>


struct my_option
{
    option gnu_opt;
    char* desc;
    my_option(const char*name, int has_arg, int *flag, int val, char* desc);
    my_option();
};

enum PinMethod
{
    FILL, SCATTER
};

struct parser
{
        char *mat_file;
        int iter;
        bool validate;
        double tol;
        char *prgname;
        int chunkHeight;
	int sigma;
	int hpcg_size;
	bool RCM_flag;
        int numOptions;
        my_option *long_options;
        option *gnuOptions;
        parser();
        ~parser();
        bool parse_arg(int argc, char **argv);
        int dump_arg();
        void help();
};

#endif
