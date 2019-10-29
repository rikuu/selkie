#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#define MIN_FRAGMENT_LEN 1.0

size_t read_rmaps(const char *filename, std::vector<std::string> *names,
        std::vector<std::vector<double> > *forward, std::vector<std::vector<double> > *reverse) {
    // Open the input file
    std::ifstream file;
    file.open(filename);
    if (!file.good()) {
        std::cerr << "Rmap file invalid." << std::endl;
        exit(1);
    }

    std::string line;
    int i = 0;
    size_t count = 0;
    while (std::getline(file, line)) {
        // Line counting:
        // 0 = header line (rmap name),
        // 1 = fragment line (list of fragments separated by tabs, first fragment is in column 3)
        // 2 = empty line

        // Add the name of the Rmap
        if (i == 0) {
            if (names)
                names->push_back(line);
        }

        // Read the fragment lengths
        if (i == 1) {
            std::vector<double> cols;
            size_t j = 0;
            while (j < line.length()) {
                size_t k = 0;
                while (j+k < line.length() && line[j+k] != ' ' && line[j+k] != '\t') {
                    k++;
                }

                double f = atof(line.substr(j, k).c_str());
                if (f > 0.0) {
                    cols.push_back(f);
                }
                
                j += k+1;
            }

            /*if (cols.size() < 3) {
                std::cout << "Invalid fragment list: " << line << std::endl;
                exit(2);
            }*/

            // Turn the fragment lengths into a sequence of cut sites
            if (forward) {
                std::vector<double> cuts;
                double s = 0.0;
                cuts.push_back(s);
                for (size_t j = 0; j < cols.size(); j++) {
                    s += cols[j];
                    if (cols[j] >= MIN_FRAGMENT_LEN)
                        cuts.push_back(s);
                }
            
                forward->push_back(cuts);
            }
            
            // Figure out the cut sites for the reverse Rmap
            if (reverse) {
                std::vector<double> cuts;
                double s = 0.0;
                cuts.push_back(s);
                for (size_t j = cols.size()-1; j != 0; j--) {
                    s += cols[j];
                    if (cols[j] >= MIN_FRAGMENT_LEN)
                        cuts.push_back(s);
                }
            
                reverse->push_back(cuts);
            }

            count++;
        }
    
        i++;
        if (i > 2)
            i = 0;
    }

    file.close();
    return count;
}
