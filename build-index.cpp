#include <string>
#include <iostream>

#include "index.hpp"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << std::endl;
        return 1;
    }

    int ell = 80;
    int mink = 5;
    char *gap_pattern = (char *) "11111111110001110110010010011101001110001010010100001010011000010111100000001100";

    const std::string outfile = "asd";
    index_t index(argv[1], ell, mink, gap_pattern, 4);

    // emphf::logger() << "Saving index" << std::endl;
    std::ofstream os(outfile, std::ios::binary);
    index.save(os);
    os.close();

    // emphf::logger() << "Done" << std::endl;
}