#include <unistd.h>

#include <string>
#include <iostream>

#include "index.hpp"

int main(int argc, char** argv) {
    std::string outfile = "";

    int ell = 80;
    int mink = 5;
    int n_parts = 1;
    double quantization = 1.;
    char *gap_pattern = (char *) "11111111110001110110010010011101001110001010010100001010011000010111100000001100";

    int c;
    while ((c = getopt(argc, argv, "o:l:s:p:q:v:k:")) >= 0) {
        if (c == 'o') outfile = optarg;
        else if (c == 'l') ell = atoi(optarg);
        else if (c == 'k') mink = atoi(optarg);
        else if (c == 's') gap_pattern = optarg;
        else if (c == 'p') n_parts = atoi(optarg);
        else if (c == 'q') quantization = atof(optarg);
        else if (c == 'v') {
            printf("1.0\n");
            return 0;
        }
    }

    if (argc == optind) {
        fprintf(stderr, "Usage: selkie-index [options] <input>\n");
        fprintf(stderr, "-o FILE     output file [stdout]\n");
        fprintf(stderr, "-l INT      ell [%d]\n", ell);
        fprintf(stderr, "-k INT      k [%d]\n", mink);
        fprintf(stderr, "-q FLOAT    quantization bin size [%.1f]\n", quantization);
        fprintf(stderr, "-p INT      construction partitions [%d]\n", n_parts);
        fprintf(stderr, "-s STR      spacing pattern\n");
        fprintf(stderr, "\"%s\"\n", gap_pattern);
        return 1;
    }

    std::string infile = argv[optind];
    index_t index(infile, ell, mink, quantization, gap_pattern, n_parts);

    std::ofstream os(outfile, std::ios::binary);
    index.save(os);
    os.close();
}
