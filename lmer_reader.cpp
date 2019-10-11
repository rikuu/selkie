#include <vector>
#include <iostream>

#include "lmer_reader.hpp"

static inline int quantize(const double d) {
    return (int) (d / 1.);
}

std::vector<std::vector<int> > extract_lmers(
        const std::vector<double> &cuts, int ell, int mink, const char *gap_pattern) {
    std::vector<std::vector<int> > lmers;

    // No el-mer can be extracted because the Rmap is too short
    if (cuts[cuts.size()-1] < ell)
        return lmers;

    // Figure out the last cut contained in the first el-mer
    size_t first_cut = 0, last_cut = 0;
    while (last_cut < cuts.size() && cuts[last_cut] < ell)
        last_cut++;

    const double chunk = (double)ell / (double)strlen(gap_pattern);

    // First cut that is not masked by the gap pattern
    int first_used_cut = -1;

    // Filter cuts based on the gap pattern
    std::vector<double> filtered_cuts;
    
    for (size_t jj = first_cut; jj < last_cut; jj++) {
        if (gap_pattern[(int)(cuts[jj]/chunk)] == '1') {
            filtered_cuts.push_back(cuts[jj]);
            if (first_used_cut < 0)
                first_used_cut = jj;
        }
    }

    // Add additional cuts if there are less than mink+1
    for (size_t jj = last_cut; jj < cuts.size() && filtered_cuts.size() <= (unsigned) mink; jj++) {
        filtered_cuts.push_back(cuts[jj]);
        if (first_used_cut < 0) {
            first_used_cut = jj;
        }
    }

    // Represent the (el,mink)-mer as a string and add it to the list of mers
    std::vector<int> lmer;
    for (size_t jj = 0; jj < filtered_cuts.size()-1; jj++) {
        double frag = filtered_cuts[jj+1]-filtered_cuts[jj];
        lmer.push_back(quantize(frag));
    }

    if (lmer.size() >= mink) {
#ifdef DEBUG
      for (int iii = 0; iii < lmer.size(); iii++) {
	std::cout << "," << lmer[iii];
      }
      std::cout << std::endl;
#endif
      lmers.push_back(lmer);
    }

    // Current offset for the gap pattern in kbp
    double start_pos = 0.0;

    // Add the rest of the (el,mink)-mers
    while (last_cut < cuts.size()) {
        // Find out the minimum shift of the gap pattern that will add or
        // remove a cut site. This can happen at any 1/* or */1 boundary of the
        // gap pattern or at the beginning/end of the gap pattern.

        // Shift needed for the first cut to drop out
        double min_shift = cuts[first_cut]-start_pos;

        // Shift needed for the next cut to enter the gap pattern region
        if (cuts[last_cut] - (start_pos+ell) < min_shift)
            min_shift = cuts[last_cut] - (start_pos+ell);

        // Minimum shift for any 1/* or */1 boundary to hit the next cut
        for (size_t jj = 0; jj < strlen(gap_pattern)-1; jj++) {
            if (gap_pattern[jj+1] != gap_pattern[jj]) {
                // A 1/* or */1 boundary
                // -> Find the next cut after that boundary and compute the shift for that cut to hit the boundary
                for (size_t ii = first_cut; ii < last_cut; ii++) {
                    if (cuts[ii] >= start_pos + (jj+1)*chunk) {
                        if (cuts[ii] - start_pos - (jj+1)*chunk < min_shift) {
                            min_shift = cuts[ii] - start_pos - (jj+1)*chunk;
                            break;
                        }
                    }
                }
            }
        }

        // Add a small constant to the shift to make sure that the cut actually changes status
        start_pos += min_shift + 0.01;

        // Find the first cut of this el-mer
        while (first_cut < cuts.size() && start_pos > cuts[first_cut])
            first_cut++;

        // Find the last cut of this el-mer
        while (last_cut < cuts.size() && start_pos+ell > cuts[last_cut])
            last_cut++;

        // Check for boundary conditions
        if (first_cut >= cuts.size() || last_cut >= cuts.size()) break;
        
        // Figure out the next el-mer
        // Filter cuts based on the gap pattern
        filtered_cuts.clear();
        first_used_cut = -1;
        for (size_t jj = first_cut; jj < last_cut; jj++) {
            if (gap_pattern[(int)((cuts[jj]-start_pos)/chunk)] == '1') {
                filtered_cuts.push_back(cuts[jj]);
                if (first_used_cut < 0)
                    first_used_cut = jj;
                }
        }
        
        // Add additional cuts if there are less than mink+1
        for (size_t jj = last_cut; jj < cuts.size() && filtered_cuts.size() <= (unsigned) mink; jj++) {
            filtered_cuts.push_back(cuts[jj]);
            if (first_used_cut < 0) {
                first_used_cut = jj;
            }
        }

        std::vector<int> lmerr;
        for (size_t jj = 0; jj < filtered_cuts.size()-1; jj++) {
            double frag = filtered_cuts[jj+1]-filtered_cuts[jj];
            lmerr.push_back(quantize(frag));
        }

	if (lmerr.size() >= mink) {
#ifdef DEBUG
	  for (int iii = 0; iii < lmerr.size(); iii++) {
	    std::cout << "," << lmerr[iii];
	  }
	  std::cout << std::endl;
#endif
	  lmers.push_back(lmerr);
	}
    }

    return lmers;
}
