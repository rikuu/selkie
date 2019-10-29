#ifndef LMER_READER_HPP
#define LMER_READER_HPP

#include <cstring>
#include <cassert>

#include <vector>
#include <string>
#include <fstream>
#include <iterator>

std::vector<std::vector<int> > extract_lmers(const std::vector<double> &cuts,
    const int ell, const int mink, const double quantization,
    const char *gap_pattern);

class lmer_iterator
    : public std::iterator<std::forward_iterator_tag, const std::vector<int> > {

public:
    lmer_iterator() : m_is(nullptr), m_pos(0) {    
    }

    lmer_iterator(FILE *is) {
        m_is = is;
    
        fseek(m_is, 0, SEEK_SET);
        fread(reinterpret_cast<char *>(&m_n), sizeof(size_t), 1, m_is);
        m_pos = ftell(m_is);

        m_line = std::vector<int>(0);

        m_m = 0;

        advance();
    }

    value_type const& operator*() const {
        return m_line;
    }

    lmer_iterator& operator++() {
        advance();
        return *this;
    }

    friend bool operator==(lmer_iterator const& lhs, lmer_iterator const& rhs) {
        if (!lhs.m_is || !rhs.m_is) {
            if (!lhs.m_is && !rhs.m_is) {
                return true;
            } else {
                return false;
            }
        }

        assert(lhs.m_is == rhs.m_is);

        return rhs.m_pos == lhs.m_pos;
    }

    friend bool operator!=(lmer_iterator const& lhs, lmer_iterator const& rhs) {
        return !(lhs == rhs);
    }

private:
    void advance() {
        // TODO: Experiment with buffering

        assert(m_is);
        fseek(m_is, m_pos, SEEK_SET);

        size_t length = 0;
        fread(reinterpret_cast<char *>(&length), sizeof(size_t), 1, m_is);
        
        m_line.resize(length);
        fread(reinterpret_cast<char *>(m_line.data()), sizeof(int), length, m_is);

        m_pos = ftell(m_is);

        m_m++;
        if (m_m == m_n)
            m_is = nullptr;
    }

    FILE *m_is;
    long m_pos;
    std::vector<int> m_line;
    size_t m_n, m_m;
};

class lmer_reader {
public:
    lmer_reader(const std::string &filename) {
        m_is = fopen(filename.c_str(), "rb");

        if (!m_is) {
            throw std::invalid_argument("Error opening " + filename);
        }

        fread(reinterpret_cast<char *>(&m_n), sizeof(size_t), 1, m_is);
    }

    ~lmer_reader() {
        fclose(m_is);
    }

    lmer_iterator begin() const {
        return lmer_iterator(m_is);
    }

    lmer_iterator end() const { return lmer_iterator(); }

    size_t size() const {
        return m_n;
    }

private:
    // noncopyble
    lmer_reader(lmer_reader &);
    lmer_reader& operator=(lmer_reader const&);
    
    size_t m_n;
    FILE *m_is;
};

#endif