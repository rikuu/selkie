#ifndef ELMERI_INDEX_HPP
#define ELMERI_INDEX_HPP

#include "streamvbyte.h"

#include <vector>
#include <fstream>
#include <algorithm>

#include "write_lmer.hpp"
#include "mergetable.hpp"

class related {
public:
    int related_id, current_pos, related_pos;
};

static bool sort_related(const related &r1, const related &r2) {
    if (r1.related_id != r2.related_id)
        return r1.related_id < r2.related_id;
    else
        return r1.current_pos < r2.current_pos;
}

struct sort_pairs {
    bool operator()(const std::pair<related, int> &left, const std::pair<related, int> &right) {
        return left.second > right.second;
    }
};

class index_t {
public:
    index_t() {}
    index_t(const std::string &in, int ell, int mink, double quantization, char *gap_pattern, int n_parts) {
        m_ell = ell;
        m_mink = mink;
        m_gap_pattern = gap_pattern;
        m_quantization = quantization;

        m_max_list_length = 0;
        m_sum_of_lengths = 0;
        m_sum_of_sizes = 0;
        
        const std::string key_file = "tmp.xxx";

        emphf::logger() << "Writing " << n_parts << " parts" << std::endl;
        m_n = write_lmers(in, key_file, n_parts, ell, mink, quantization, gap_pattern);

        if (n_parts != 1) {
            emphf::logger() << "Merging parts" << std::endl;
            m_n = multiwaymerge(key_file, n_parts);
        }

        emphf::logger() << "Constructing MPHF" << std::endl;
        mphf_t(key_file).swap(m_mphf);

        emphf::logger() << "Constructing mergetable" << std::endl;
        m_merge = mergetable_t(m_n, key_file, m_mphf);

        m_max_list_length = 8;
        size_t bytes_alloc = 3*m_n/2;
        size_t lengths_alloc = bytes_alloc;
        
        uint32_t *delta_buffer = reinterpret_cast<uint32_t *>(calloc(m_max_list_length, sizeof(uint32_t)));
        m_lists = reinterpret_cast<uint8_t *>(calloc(bytes_alloc, sizeof(uint8_t)));

        m_lengths = sdsl::bit_vector(bytes_alloc, 0);
        m_sizes = sdsl::bit_vector(bytes_alloc, 0);

        m_lengths[0] = 1;
        m_sizes[0] = 1;

        const size_t step = m_n / n_parts;

        size_t max_bytes = 0;
        emphf::logger() << "Gathering lists" << std::endl;
        for (size_t i = 0; i < m_n-step; i += step) {
            std::vector<std::vector<uint32_t> > lists = gather_lists(in, i, i+step);

            size_t length = 0;
            bool max_list_changed = false;
            for (size_t i = 0; i < lists.size(); i++) {
                if (lists[i].size() > m_max_list_length) {
                    m_max_list_length = lists[i].size();
                    max_list_changed = true;
                }
                
                max_bytes += streamvbyte_max_compressedbytes(lists[i].size());
                length += lists[i].size();
            }

            if (max_list_changed)
                delta_buffer = (uint32_t*) realloc(delta_buffer, m_max_list_length * sizeof(uint32_t));
            
            if (max_bytes >= bytes_alloc) {
                // TODO: prealloc more than required
                bytes_alloc = max_bytes;
                m_lists = (uint8_t*) realloc(m_lists, bytes_alloc * sizeof(uint8_t));
                m_sizes.resize(bytes_alloc);
            }
            
            if (m_sum_of_lengths + length >= lengths_alloc) {
                // TODO: for sure prealloc more
                lengths_alloc = m_sum_of_lengths + length;
                m_lengths.resize(lengths_alloc);
            }

            compress_lists(lists, delta_buffer);
        }

        free(delta_buffer);

        m_lengths.resize(m_sum_of_lengths+2);
        m_sizes.resize(m_sum_of_sizes+2);

        m_lengths_select = new sdsl::bit_vector::select_1_type(&m_lengths);
        m_sizes_select = new sdsl::bit_vector::select_1_type(&m_sizes);
    }

    void save(std::ostream &os) const;
    void load(std::istream &is);

    uint8_t *get_values(const std::vector<int> q, size_t *out_size);

    size_t get_related(const std::vector<double> rmap, std::vector<std::pair<related, unsigned int> > *counts);

    std::vector<std::vector<uint32_t> > gather_lists(const std::string &in, size_t start, size_t end);
    void compress_lists(std::vector<std::vector<uint32_t> > lists, uint32_t *delta_buffer);

private:
    sdsl::bit_vector m_lengths, m_sizes;
    sdsl::bit_vector::select_1_type *m_lengths_select, *m_sizes_select;

    uint8_t *m_lists;
    mphf_t m_mphf;
    mergetable_t m_merge;

    int m_ell, m_mink;
    size_t m_n, m_m;
    double m_quantization;
    char* m_gap_pattern;

    size_t m_sum_of_lengths, m_sum_of_sizes, m_max_list_length;
};

uint8_t *index_t::get_values(const std::vector<int> q, size_t *out_size) {
    const uint64_t u = m_mphf.lookup(q);

    uint64_t s, l, l1;
    if (!m_merge.merged(u)) {
        uint64_t r = m_merge.rank(u);
        l = m_lengths_select->select(u+1 - r);
        l1 = m_lengths_select->select(u+2 - r);
        s = m_sizes_select->select(u+1 - r);
    } else {
        uint32_t v = m_merge.end(u);
        uint64_t r = m_merge.rank(v);
        l = m_lengths_select->select(v+1 - r);
        l1 = m_lengths_select->select(v+2 - r);
        s = m_sizes_select->select(v+1 - r);
    }

    *out_size = l1 - l;
    return m_lists + s;
}

size_t index_t::get_related(const std::vector<double> rmap, std::vector<std::pair<related, unsigned int> > *counts) {
    std::vector<std::vector<int> > lmers = extract_lmers(rmap, m_ell, m_mink, m_quantization, m_gap_pattern);

    int exclude_id = counts->size() > 0 ? (*counts)[0].first.related_id : -1;

    // TODO: Get pre-allocated buffer as input
    uint32_t *delta_buffer = (uint32_t *) calloc(m_max_list_length, sizeof(uint32_t));

    std::vector<related> friends;
    std::set<uint64_t> used_lmers;
    for (auto it = lmers.begin(); it != lmers.end(); ++it) {
        uint64_t u = m_mphf.lookup(*it);
        if (used_lmers.find(u) != used_lmers.end())
            continue;

        used_lmers.insert(u);

        size_t out_size = 0;
        uint8_t *v = get_values(*it, &out_size);
        streamvbyte_decode(v, delta_buffer, out_size);
    
        for (size_t i = 1; i < out_size; i++) {
            delta_buffer[i] = delta_buffer[i-1] + delta_buffer[i];
        }

        for (size_t i = 0; i < out_size; i++) {
            related r;
            r.related_id = delta_buffer[i];
            r.current_pos = 0;
            r.related_pos = 0;
            friends.push_back(r);
        }
    }

    free(delta_buffer);

    // Count how many times each related Rmap is contained in the related Rmaps set
    std::sort(friends.begin(), friends.end(), sort_related);

    related r;
    r.related_id = -1;
    r.related_pos = -1;
    r.current_pos = -1;
    unsigned int c = 0;

    for (size_t k = 0; k < friends.size(); k++) {
        if (r.related_id >= 0) {
            if (r.related_id != friends[k].related_id) {
                if (r.related_id != exclude_id) {
                    counts->push_back(std::make_pair(r, c));
                }
                
                r = friends[k];
                c = 1;
            } else {
                c++;
            }
        } else {
            r = friends[k];
            c = 1;
        }
    }

    if (r.related_id != exclude_id) {
        counts->push_back(std::make_pair(r, c));
    }

    // Sort the related Rmaps based on the number of shared mers
    std::sort(counts->begin(), counts->end(), sort_pairs());

    return counts->size();
}

std::vector<std::vector<uint32_t> > index_t::gather_lists(const std::string &in, size_t start, size_t end) {
    std::vector<std::vector<uint32_t> > lists(end - start + 1);

    std::vector<std::vector<double> > forward, reverse;
    size_t rmap_count = read_rmaps(in.c_str(), 0, &forward, &reverse);
    for (size_t i = 0; i < rmap_count; i++) {
        std::vector<std::vector<int> > lmers_f = extract_lmers(forward[i], m_ell, m_mink, m_quantization, m_gap_pattern);

        for (auto it = lmers_f.begin(); it != lmers_f.end(); ++it) {
            const std::vector<int> key = *it;
            const uint64_t u = m_mphf.lookup(key);
            if (u >= m_n)
                continue;

            const uint64_t v = m_merge.merged(u) ? m_merge.end(u) : u;
            
            if (v < start || v > end)
                continue;

            if (lists[v - start].size() != 0 && lists[v-start].back() == (uint32_t) i*2)
                continue;

            lists[v - start].push_back((uint32_t) i*2);
        }

        std::vector<std::vector<int> > lmers_r = extract_lmers(reverse[i], m_ell, m_mink, m_quantization, m_gap_pattern);
    
        for (auto it = lmers_r.begin(); it != lmers_r.end(); ++it) {
            const std::vector<int> key = *it;
            const uint64_t u = m_mphf.lookup(key);
            if (u >= m_n)
                continue;

            const uint64_t v = m_merge.merged(u) ? m_merge.end(u) : u;

            if (v < start || v > end)
                continue;

            if (lists[v - start].size() != 0 && lists[v-start].back() == (uint32_t) (i*2)+1)
                continue;

            lists[v - start].push_back((uint32_t) (i*2)+1);
        }
    }

    return lists;
}

void index_t::compress_lists(std::vector<std::vector<uint32_t> > lists, uint32_t *delta_buffer) {
    for (size_t i = 0; i < lists.size(); i++) {
        if (lists[i].size() == 0)
            continue;

        std::sort(lists[i].begin(), lists[i].end());
        
        uint32_t prev = 0, skipped = 0;
        for (size_t j = 0; j < lists[i].size(); j++) {
            assert(lists[i][j] >= prev);
            uint32_t d = lists[i][j] - prev;

            if (d == 0 && j != 0) {
                skipped++;
                continue;
            }
            
            delta_buffer[j] = d;
            prev = lists[i][j];
        }
        
        size_t size = streamvbyte_encode(delta_buffer, lists[i].size() - skipped, m_lists + m_sum_of_sizes);

        m_sum_of_lengths += lists[i].size() - skipped;
        m_sum_of_sizes += size;

        m_lengths[m_sum_of_lengths] = 1;
        m_sizes[m_sum_of_sizes] = 1;
    }
}

void index_t::save(std::ostream& os) const {
    std::streampos prev = os.tellp();
    m_mphf.save(os);
    
    std::cout << "mphf\t\t" << os.tellp() - prev << std::endl;
    prev = os.tellp();

    os.write(reinterpret_cast<char const*>(&m_n), sizeof(size_t));
    os.write(reinterpret_cast<char const*>(&m_m), sizeof(size_t));
    os.write(reinterpret_cast<char const*>(&m_sum_of_lengths), sizeof(size_t));
    os.write(reinterpret_cast<char const*>(&m_sum_of_sizes), sizeof(size_t));
    os.write(reinterpret_cast<char const*>(&m_max_list_length), sizeof(size_t));

    os.write(reinterpret_cast<char const*>(&m_ell), sizeof(int));
    os.write(reinterpret_cast<char const*>(&m_mink), sizeof(int));
    os.write(reinterpret_cast<char const*>(&m_quantization), sizeof(double));

    const size_t pattern_length = strlen(m_gap_pattern);
    os.write(reinterpret_cast<char const*>(&pattern_length), sizeof(size_t));
    os.write(reinterpret_cast<char const*>(m_gap_pattern), pattern_length);

    std::cout << "constants\t" << os.tellp() - prev << std::endl;
    prev = os.tellp();

    m_lengths.serialize(os);

    std::cout << "lengths\t\t" << os.tellp() - prev << std::endl;
    prev = os.tellp();

    m_sizes.serialize(os);

    std::cout << "sizes\t\t" << os.tellp() - prev << std::endl;
    prev = os.tellp();

    os.write(reinterpret_cast<char const*>(m_lists), (std::streamsize)(sizeof(m_lists[0]) * m_sum_of_sizes));
    
    std::cout << "lists\t\t" << os.tellp() - prev << std::endl;
    prev = os.tellp();

    m_merge.save(os);
    std::cout << "merge\t\t" << os.tellp() - prev << std::endl;
}

void index_t::load(std::istream& is) {
    m_mphf.load(is);
    
    is.read(reinterpret_cast<char*>(&m_n), sizeof(size_t));
    is.read(reinterpret_cast<char*>(&m_m), sizeof(size_t));
    is.read(reinterpret_cast<char*>(&m_sum_of_lengths), sizeof(size_t));
    is.read(reinterpret_cast<char*>(&m_sum_of_sizes), sizeof(size_t));
    is.read(reinterpret_cast<char*>(&m_max_list_length), sizeof(size_t));

    is.read(reinterpret_cast<char*>(&m_ell), sizeof(int));
    is.read(reinterpret_cast<char*>(&m_mink), sizeof(int));
    is.read(reinterpret_cast<char*>(&m_quantization), sizeof(double));

    size_t len_pattern = 0;
    is.read(reinterpret_cast<char*>(&len_pattern), sizeof(size_t));
    m_gap_pattern = reinterpret_cast<char*>(calloc(len_pattern+1, sizeof(char)));
    is.read(reinterpret_cast<char*>(m_gap_pattern), len_pattern);

    m_lengths.load(is);
    m_lengths_select = new sdsl::bit_vector::select_1_type(&m_lengths);

    m_sizes.load(is);
    m_sizes_select = new sdsl::bit_vector::select_1_type(&m_sizes);

    m_lists = reinterpret_cast<uint8_t *>(calloc(m_sum_of_sizes, sizeof(uint8_t)));
    is.read(reinterpret_cast<char*>(m_lists), (std::streamsize)(sizeof(m_lists[0]) * m_sum_of_sizes));

    m_merge.load(is);
}

#endif
