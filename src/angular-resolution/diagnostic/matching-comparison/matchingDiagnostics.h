#ifndef MATCHING_DIAGNOSTICS_H
#define MATCHING_DIAGNOSTICS_H

#include <vector>
#include <set>

using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_f = ROOT::VecOps::RVec<float>;

// Count valid matches (index >= 0)
int countValidMatches(const Vec_i& matches) {
    int n = 0;
    for (auto m : matches) {
        if (m >= 0) n++;
    }
    return n;
}

// Count unique valid targets
int countUniqueMatches(const Vec_i& matches) {
    std::set<int> unique_targets;
    for (auto m : matches) {
        if (m >= 0) unique_targets.insert(m);
    }
    return unique_targets.size();
}

// Check if there are duplicate valid targets
bool hasDuplicateMatches(const Vec_i& matches) {
    std::set<int> unique_targets;
    int n_valid = 0;
    for (auto m : matches) {
        if (m >= 0) {
            n_valid++;
            unique_targets.insert(m);
        }
    }
    return (int)unique_targets.size() < n_valid;
}

// Exact vector equality
bool exactMatchAgreement(const Vec_i& a, const Vec_i& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

// Count per-jet assignment agreement
int countAssignmentAgreement(const Vec_i& a, const Vec_i& b) {
    if (a.size() != b.size()) return 0;
    int n = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == b[i]) n++;
    }
    return n;
}

// Build transitive mapping: reco->gen and gen->parton => reco->parton
Vec_i transitiveMatch(const Vec_i& reco_to_gen, const Vec_i& gen_to_parton) {
    Vec_i out;
    out.reserve(reco_to_gen.size());
    for (auto g : reco_to_gen) {
        if (g >= 0 && g < (int)gen_to_parton.size()) out.push_back(gen_to_parton[g]);
        else out.push_back(-1);
    }
    return out;
}

#endif
