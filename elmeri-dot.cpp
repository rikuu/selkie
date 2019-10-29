#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "lmer_reader.hpp"
#include "rmap.hpp"
#include "index.hpp"

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
    
    // Gathering statistics on the related Rmaps
    int singletons = 0;
    int edges = 0;
    int nodes = 0;

    std::cout << "graph {\n";

    std::vector<std::vector<double> > forward;
    read_rmaps(in_filename, NULL, &forward, NULL);
    for (size_t i = 0; i < forward.size(); i++) {
        // Get the related Rmaps
        std::vector<std::pair<related, unsigned int> > counts;
        size_t counts_size = index.get_related(forward[i], &counts);
        
        // Output the related Rmaps
        nodes++;
        int degree = 0;

        for (size_t k = 0; k < counts_size; k++) {
            if (counts[k].second < count_thrs)
                continue;

            degree++;
            // Only print each edge once, i.e. when the current Rmap has a smaller index than the related one
            if (counts[k].first.related_id/2 > i) {
                // <current Rmap> -- <related Rmap> [weight=<number of shared mers>]
                std::cout << i << " -- " << counts[k].first.related_id/2 << " [weight=" << counts[k].second << "];\n";
                edges++;
            }
        }

        if (degree == 0) {
            singletons++;
        }
    }


    std::cout << "}";
    std::cout << std::endl;

    // Some statistics of the related Rmaps
    fprintf(stderr, "Nodes: %d\n", nodes);  // Number of Rmaps
    fprintf(stderr, "Edges: %d\n", edges);  // Number of related Rmaps pairs
    fprintf(stderr, "Singleton nodes: %d\n", singletons);  // Number of Rmaps with no related Rmaps

    // fprintf(stderr, "lmers: %d\n", index.n);

    return 0;
}
