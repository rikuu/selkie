#ifndef ELMERI_MERGETABLE_HPP
#define ELMERI_MERGETABLE_HPP

#include <cstring>

#include <vector>
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <random>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "thirdparty/emphf/common.hpp"
#include "thirdparty/emphf/mphf.hpp"
#include "thirdparty/emphf/base_hash.hpp"
#include "thirdparty/emphf/hypergraph.hpp"
#include "thirdparty/emphf/hypergraph_sorter_seq.hpp"

#include "lmer_reader.hpp"

//#define DEBUG

class mphf_t {
public:
    struct lmer_adaptor {
        emphf::byte_range_t operator()(std::vector<int> const& s);
    };

    mphf_t() {}
    mphf_t(const std::string &key_file);

    uint64_t lookup(const std::vector<int> &q);
    void save(std::ostream &os) const;
    void load(std::istream &is);

    void swap(mphf_t& other) {
        //std::swap(m_adaptor, other.m_adaptor);
        m_mphf.swap(other.m_mphf);
    }

private:
    lmer_adaptor m_adaptor;
    emphf::mphf<emphf::jenkins64_hasher> m_mphf;
};

class mergetable_t {
public:
    mergetable_t() {}
    
    mergetable_t(const size_t n) {
        m_bv = new sdsl::bit_vector(n, 0);
        m_rank = new sdsl::bit_vector::rank_1_type(m_bv);
        m_index = std::vector<uint64_t>();
    }

    mergetable_t(const size_t n, const std::string &key_file, mphf_t mphf);

    inline size_t size() const {
        return m_bv->size();
    }

    inline bool merged(const uint64_t u) const {
      return (*m_bv)[u];
    }

    inline uint64_t end(const uint64_t u) const {
        return m_index[rank(u)];
    }

    inline size_t rank(const uint64_t u) const {
        return m_rank->rank(u) ;
    }

    void save(std::ostream& os) const {
        m_bv->serialize(os);

        size_t n = m_index.size();
        os.write(reinterpret_cast<char const*>(&n), sizeof(size_t));
        os.write(reinterpret_cast<char const*>(m_index.data()), (std::streamsize)(sizeof(m_index[0]) * n));
    }

    void load(std::istream& is) {
        m_bv = new sdsl::bit_vector(100, 0);
        m_bv->load(is);
        m_rank = new sdsl::bit_vector::rank_1_type(m_bv);

        size_t n = 0;
        is.read(reinterpret_cast<char*>(&n), sizeof(size_t));

        m_index.resize(n); // = reinterpret_cast<uint64_t *>(calloc(n, sizeof(uint64_t)));
        is.read(reinterpret_cast<char*>(m_index.data()), (std::streamsize)(sizeof(m_index[0]) * n));
    }

private:
    sdsl::bit_vector *m_bv;
    sdsl::bit_vector::rank_1_type *m_rank;
    std::vector<uint64_t> m_index;
};

#define MER_SIMILARITY_THRS 10

emphf::byte_range_t mphf_t::lmer_adaptor::operator()(std::vector<int> const& s) {
    const uint8_t* start = reinterpret_cast<uint8_t const*>(s.data());
    const uint8_t* end = start + (s.size() * sizeof(int));
    return emphf::byte_range_t(start, end);
}

mphf_t::mphf_t(const std::string &key_file) {
    lmer_reader keys(key_file.c_str());
    size_t n = keys.size() - 1;

    size_t max_nodes = size_t(std::ceil(double(n) * 1.23)) + 2;
    if (max_nodes >= uint64_t(1) << 32) {
        emphf::hypergraph_sorter_seq<emphf::hypergraph<uint64_t>> sorter;
        emphf::mphf<emphf::jenkins64_hasher>(sorter, n, keys, m_adaptor).swap(m_mphf);
    } else {
        emphf::hypergraph_sorter_seq<emphf::hypergraph<uint32_t>> sorter;
        emphf::mphf<emphf::jenkins64_hasher>(sorter, n, keys, m_adaptor).swap(m_mphf);
    }
}

uint64_t mphf_t::lookup(const std::vector<int> &q) {
    return m_mphf.lookup(q, m_adaptor);
}

void mphf_t::save(std::ostream &os) const {
    m_mphf.save(os);
}

void mphf_t::load(std::istream &is) {
    m_mphf.load(is);
}

#define ROOT 0xffffffffffffffff

mergetable_t::mergetable_t(const size_t n, const std::string &key_file, mphf_t mphf) {
    uint64_t *tree = (uint64_t *) calloc(sizeof(uint64_t), n);
    uint32_t *sizes = (uint32_t *) calloc(sizeof(uint32_t), n);
    
    // memset(sizes, 1, n);
    for (size_t i = 0; i < n; i++) {
        tree[i] = ROOT;
        sizes[i] = 1;
    }

    lmer_reader keys(key_file.c_str());
    std::vector<int> *u2lmer = new std::vector<int>[n];
    //std::unordered_map<uint64_t, std::vector<int> > u2lmer;
    for (auto it = keys.begin(); it != keys.end(); ++it) {
      std::vector<int> lmer = *it;
      uint64_t u = mphf.lookup(lmer);
      assert(u < n);
      u2lmer[u] = lmer;
    }


    for (int round = 0; round < 2; round++) {
        for (auto it = keys.begin(); it != keys.end(); ++it) {
            std::vector<int> lmer = *it;
            uint64_t u = mphf.lookup(lmer);
            assert(u < n);

            uint64_t root_u = u;
            while (tree[root_u] != ROOT)
                root_u = tree[root_u];
            
            if (sizes[root_u]+1 > MER_SIMILARITY_THRS)
                continue;

            for (size_t col = 0; col < lmer.size(); col++) {
                {
                    lmer[col]--;
                    uint64_t v = mphf.lookup(lmer);

		    bool ok = true;
		    if (v >= n)
		      ok = false;
		    else {
		      std::vector<int> vlmer = u2lmer[v];
		      if (lmer.size() != vlmer.size()) {
			ok = false;
		      } else {
			for(int i = 0; i < lmer.size(); i++) {
			  if (lmer[i] != vlmer[i]) {
			    ok = false;
			    break;
			  }
			}
		      }
#ifdef DEBUG
		      if (!ok) {
			std::cout << "False merge:" << u << " " << v << std::endl;
			for(int i = 0; i < lmer.size(); i++) {
			  std::cout << "," << lmer[i];
			}
			std::cout << std::endl;
			for(int i = 0; i < vlmer.size(); i++) {
			  std::cout << "," << vlmer[i];
			}
			std::cout << std::endl;
		      }
#endif
		    }
		    lmer[col]++;

                    if (ok && v < n) {
		        uint64_t root_v = v;
                        while (tree[root_v] != ROOT)
                            root_v = tree[root_v];
                        
                        if (root_u == root_v)
                            continue;

                        if (sizes[root_u]+sizes[root_v] > MER_SIMILARITY_THRS)
                            continue;
                        
                        if (sizes[root_u] < sizes[root_v]) {
                            tree[root_u] = root_v;
                            sizes[root_v] += sizes[root_u];
			    root_u = root_v;
                        } else {
                            tree[root_v] = root_u;
                            sizes[root_u] += sizes[root_v];
                        }
                    }
                }

                {
                    lmer[col]++;
                    uint64_t v = mphf.lookup(lmer);
		    bool ok = true;
		    if (v >= n)
		      ok = false;
		    else {
		      std::vector<int> vlmer = u2lmer[v];
		      if (lmer.size() != vlmer.size()) {
			ok = false;
		      } else {
			for(int i = 0; i < lmer.size(); i++) {
			  if (lmer[i] != vlmer[i]) {
			    ok = false;
			    break;
			  }
			}
		      }
#ifdef DEBUG
		      if (!ok) {
			std::cout << "False merge:" << u << " " << v << std::endl;
			for(int i = 0; i < lmer.size(); i++) {
			  std::cout << "," << lmer[i];
			}
			std::cout << std::endl;
			for(int i = 0; i < vlmer.size(); i++) {
			  std::cout << "," << vlmer[i];
			}
			std::cout << std::endl;
		      }
#endif
		    }

                    lmer[col]--;

                    if (ok && v < n) {
		        uint64_t root_v = v;
                        while (tree[root_v] != ROOT)
                            root_v = tree[root_v];
                        
                        if (root_u == root_v)
                            continue;
                        
                        if (sizes[root_u]+sizes[root_v] > MER_SIMILARITY_THRS)
                            continue;
                        
                        if (sizes[root_u] < sizes[root_v]) {
                            tree[root_u] = root_v;
                            sizes[root_v] += sizes[root_u];
			    root_u = root_v;
                        } else {
                            tree[root_v] = root_u;
                            sizes[root_u] += sizes[root_v];
                        }
                    }
                }
            }
        }
    }

    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        if (tree[i] == ROOT)
            continue;
        m++;

        uint64_t root = tree[i];
        while (tree[root] != ROOT)
            root = tree[root];
        
        tree[i] = root;
    }
    
    std::cout << n << " " << m << std::endl;

    m_bv = new sdsl::bit_vector(n, 0);
    for (size_t i = 0; i < n; i++) {
        if (tree[i] == ROOT)
            continue;
        
        (*m_bv)[i] = 1;
        m_index.push_back(tree[i]);
    }

    m_rank = new sdsl::bit_vector::rank_1_type(m_bv);

    // emphf::logger() << "Finished merging" << std::endl;

    free(tree);
    free(sizes);
}

#endif
