#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <boost/algorithm/string.hpp>

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
            std::vector<std::string> cols;
            boost::split(cols, line, [](char c) {return c == '\t' || c == ' ';});
            if (cols.size() < 3) {
                std::cout << "Invalid fragment list: " << line << std::endl;
                exit(2);
            }

            // Turn the fragment lengths into a sequence of cut sites
            if (forward) {
                std::vector<double> cuts;
                double s = 0.0;
                cuts.push_back(s);
                for (size_t j = 0; j < cols.size(); j++) {
                    if (atof(cols[j].c_str()) <= 0.0)
                        continue;
                    s += atof(cols[j].c_str());
                    if (atof(cols[j].c_str()) >= MIN_FRAGMENT_LEN)
                        cuts.push_back(s);
                }
            
                forward->push_back(cuts);
            }
            
            // Figure out the cut sites for the reverse Rmap
            if (reverse) {
                std::vector<double> cuts;
                double s = 0.0;
                cuts.push_back(s);
                for (size_t j = cols.size()-1; j >= 3; j--) {
                    if (atof(cols[j].c_str()) <= 0.0)
                        continue;
                    s += atof(cols[j].c_str());
                    if (atof(cols[j].c_str()) >= MIN_FRAGMENT_LEN)
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
