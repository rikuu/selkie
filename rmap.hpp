#ifndef SELKIE_RMAP_HPP
#define SELKIE_RMAP_HPP

#include <vector>
#include <string>

size_t read_rmaps(const char *filename, std::vector<std::string> *names,
        std::vector<std::vector<double> > *forward, std::vector<std::vector<double> > *reverse);

#endif
