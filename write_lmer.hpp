#ifndef SELKIE_WRITE_LMER_HPP
#define SELKIE_WRITE_LMER_HPP

#include <cstring>

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>

#include "lmer_reader.hpp"
#include "rmap.hpp"

static bool cmp_lmer(const std::vector<int> &a, const std::vector<int> &b) {
    if (a.size() == 0) return false;
    if (b.size() == 0) return true;

    size_t minsize = (a.size() < b.size()) ? a.size() : b.size();
    for (size_t i = 0; i < minsize; i++) {
        if (a[i] < b[i])
            return true;
        else if (a[i] > b[i])
            return false;
    }

    return (a.size() < b.size());
}

static bool eq_lmer(const std::vector<int> &a, const std::vector<int> &b) {
    if (a.size() != b.size())
        return false;
    
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != b[i])
            return false;
    }
    
    return true;
}

size_t write_lmers(const std::string &in, const std::string &out,
        const int n_parts, const int ell, const int mink, const double quantization,
        const char* gap_pattern) {
    std::vector<std::vector<double> > forward, reverse;
    size_t rmap_count = read_rmaps(in.c_str(), 0, &forward, &reverse);
    
    size_t all_n = 0;
    size_t rp = rmap_count / n_parts;
    for (int p = 0; p < n_parts; p++) {
        const std::string part_out = out + ((n_parts != 1) ? "." + std::to_string(p) : "");

        std::vector<std::vector<int> > keys;
        for (size_t i = p * rp; i < std::min((p+1) * rp, rmap_count); i++) {
            std::vector<std::vector<int> > lmers_f = extract_lmers(forward[i], ell, mink, quantization, gap_pattern);
        
            for (auto it = lmers_f.begin(); it != lmers_f.end(); ++it) {
                std::vector<int> key = *it;
                keys.push_back(key);
            }

            std::vector<std::vector<int> > lmers_r = extract_lmers(reverse[i], ell, mink, quantization, gap_pattern);
        
            for (auto it = lmers_r.begin(); it != lmers_r.end(); ++it) {
                std::vector<int> key = *it;
                keys.push_back(key);
            }
        }

        std::sort(keys.begin(), keys.end(), cmp_lmer);

        size_t n = 0;
        for (size_t i = 0; i < keys.size(); i++) {
            if (i > 0 && eq_lmer(keys[i], keys[i-1]))
                continue;

            // if (keys[i].size() < (unsigned) mink)
            //    continue;

            n++;
        }

        all_n += n;

        std::ofstream os(part_out, std::ios::binary);
        os.write(reinterpret_cast<const char*>(&n), sizeof(size_t));

        for (size_t i = 0; i < keys.size(); i++) {
            if (i > 0 && eq_lmer(keys[i], keys[i-1]))
                continue;
            
            // if (keys[i].size() < (unsigned) mink)
            //    continue;

            size_t len = keys[i].size();
            os.write(reinterpret_cast<const char*>(&len), sizeof(size_t));
            os.write(reinterpret_cast<const char*>(keys[i].data()), sizeof(int) * keys[i].size());
        }

        os.close();
    }

    return all_n;
}

size_t multiwaymerge(const std::string &out, const int n_parts) {
    size_t n = 0;
    
    std::vector<FILE*> files;
    std::vector<lmer_iterator> readers;
    for (int p = 0; p < n_parts; p++) {
        const std::string part_out = out + "." + std::to_string(p);
        files.push_back(fopen(part_out.c_str(), "rb"));
        readers.push_back(lmer_iterator(files[p]));
    }

    auto end = lmer_iterator();

    std::ofstream os(out, std::ios::binary);
    os.write(reinterpret_cast<const char*>(&n), sizeof(size_t));

    size_t correct_n = 0;
    std::vector<int> prev(0);
    while (readers.size() > 0) {
        auto min_lmer = *(readers[0]);
        size_t min_reader = 0;

        for (size_t p = 1; p < readers.size(); p++) {
            auto lmer = *(readers[p]);
            if (cmp_lmer(lmer, min_lmer)) {
                min_lmer = lmer;
                min_reader = p;
            }
        }
        
        if (readers[min_reader] == end)
            readers.erase(readers.begin() + min_reader);
        else
            ++readers[min_reader];

        if (eq_lmer(min_lmer, prev))
            continue;

        prev = min_lmer;

        correct_n++;
        size_t len = min_lmer.size();
        os.write(reinterpret_cast<const char*>(&len), sizeof(size_t));
        os.write(reinterpret_cast<const char*>(min_lmer.data()), (std::streamsize)(sizeof(int) * min_lmer.size()));
    }

    for (int p = 0; p < n_parts; p++)
        fclose(files[p]);

    os.seekp(0);
    os.write(reinterpret_cast<const char*>(&correct_n), sizeof(size_t));

    os.close();

    return correct_n;
}

#endif
