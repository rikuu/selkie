#include <unistd.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "index.hpp"
#include "rmap.hpp"

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cout << "Invalid arguments." << std::endl;
        return 1;
    }

    char *in_filename = argv[1];
    char *index_filename = argv[2];
    int count_thrs = atoi(argv[3]);
    
    std::ifstream is(std::string(index_filename), std::ios::binary);
    index_t index;
    index.load(is);

    std::vector<std::string> names;
    std::vector<std::vector<double> > forward;
    read_rmaps(in_filename, &names, &forward, NULL);
    for (size_t i = 0; i < forward.size(); i++) {
        std::vector<std::pair<related, unsigned int> > counts;
        size_t counts_size = index.get_related(forward[i], &counts);
        
        size_t filtered = 0;
        for (size_t k = 0; k < counts_size; k++) {
            if (counts[k].second < count_thrs || counts[k].first.related_id/2 <= i)
                continue;
            filtered++;
        }

        if (filtered == 0)
            continue;

        std::cout << names[i] << ":";
        for (size_t k = 0; k < counts_size; k++) {
            if (counts[k].second < count_thrs || counts[k].first.related_id/2 <= i)
                continue;

            const std::string orientation = ((counts[k].first.related_id % 2) == 0) ? "f" : "r";
            std::cout << " " << names[counts[k].first.related_id/2] << "," << orientation;
        }

        std::cout << std::endl;
    }

    return 0;
}