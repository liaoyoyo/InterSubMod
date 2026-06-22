/**
 * @file PhyloLabeler.cpp
 * @brief phylo-v4.1 hierarchical cluster labeling (Phase B1: coarse) — C++ port of phylo_v4.py.
 *
 * Faithful port of the Python reference v4_label (coarse, want_other=False):
 *   peel(upstream) -> UPGMA tree -> recursive { descend-quarantine to balanced node ->
 *   between/within r>=SEP_MIN -> per-subgroup column-permutation BERNOULLI null95 gate ->
 *   pay-for-itself quarantine -> hierarchical labels (larger child = lower number) } ->
 *   post-filter terminal groups < MIN_SZ to "outlier".
 *
 * RNG: std::mt19937_64 with a pinned compile-time seed (Hard-Gate determinism). Cross-language
 * equivalence with numpy PCG64 is statistical (verified by the B3 harness), not bit-identical.
 */
#include "core/PhyloLabeler.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <set>

#include "core/HierarchicalClustering.hpp"

namespace InterSubMod {

PhyloLabeler::PhyloLabeler(const PhyloConfig& cfg) : cfg_(cfg) {
}

// ---------------------------------------------------------------------------
// BERNOULLI distance recompute (matches DistanceMatrix::calculate_bernoulli /
// Python cluster_redesign.bernoulli_dist). Vectorized over the read x CpG matrix.
//   weight(p)=2|p-0.5|, delta=p_i+p_j-2 p_i p_j, Dist=sum(w_i w_j delta)/sum(w_i w_j).
//   invalid (common<min_sz or zero weight) -> -1; diagonal -> 0.
// ---------------------------------------------------------------------------
Eigen::MatrixXd PhyloLabeler::bernoulli_dist(const Eigen::MatrixXd& meth) const {
    const int n = static_cast<int>(meth.rows());
    const int C = static_cast<int>(meth.cols());
    Eigen::MatrixXd Mk(n, C), Pz(n, C), Aw(n, C), Bw(n, C);
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < C; ++c) {
            double p = meth(i, c);
            bool ok = !std::isnan(p);
            Mk(i, c) = ok ? 1.0 : 0.0;
            double pz = ok ? p : 0.0;
            Pz(i, c) = pz;
            double aw = ok ? 2.0 * std::abs(pz - 0.5) : 0.0;  // 0 weight if missing
            Aw(i, c) = aw;
            Bw(i, c) = aw * pz;
        }
    }
    Eigen::MatrixXd sw = Aw * Aw.transpose();
    Eigen::MatrixXd swd = Bw * Aw.transpose() + Aw * Bw.transpose() - 2.0 * (Bw * Bw.transpose());
    Eigen::MatrixXd common = Mk * Mk.transpose();
    Eigen::MatrixXd Dm(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                Dm(i, j) = 0.0;
            } else if (common(i, j) < cfg_.min_common || sw(i, j) < 1e-9) {
                Dm(i, j) = -1.0;
            } else {
                Dm(i, j) = swd(i, j) / sw(i, j);
            }
        }
    }
    return Dm;
}

// ---------------------------------------------------------------------------
// between(S1,S2)/within(S1∪S2) mean ratio over valid (>=0) entries. -1 if undefined.
// ---------------------------------------------------------------------------
double PhyloLabeler::between_within(const Eigen::MatrixXd& s, const std::vector<int>& S1,
                                    const std::vector<int>& S2) const {
    double bet_sum = 0.0;
    int bet_n = 0;
    for (int a : S1)
        for (int b : S2) {
            double d = s(a, b);
            if (d >= 0.0) {
                bet_sum += d;
                ++bet_n;
            }
        }
    double wit_sum = 0.0;
    int wit_n = 0;
    auto add_within = [&](const std::vector<int>& S) {
        for (size_t a = 0; a < S.size(); ++a)
            for (size_t b = a + 1; b < S.size(); ++b) {
                double d = s(S[a], S[b]);
                if (d >= 0.0) {
                    wit_sum += d;
                    ++wit_n;
                }
            }
    };
    add_within(S1);
    add_within(S2);
    if (bet_n == 0 || wit_n == 0) return -1.0;
    double wm = wit_sum / wit_n;
    if (wm <= 1e-6) return -1.0;
    return (bet_sum / bet_n) / wm;
}

// ---------------------------------------------------------------------------
// Build UPGMA tree (HierarchicalClustering) and extract desc (leaf indices per node)
// and ch (children node ids; {-1,-1} for leaves). Node ids: leaves 0..n-1, internal n..2n-2.
// ---------------------------------------------------------------------------
void PhyloLabeler::build_desc_ch(const Eigen::MatrixXd& dist, std::vector<std::vector<int>>& desc,
                                 std::vector<std::pair<int, int>>& ch, int& root, int n) const {
    desc.assign(2 * n - 1, {});
    ch.assign(2 * n - 1, {-1, -1});
    root = -1;
    if (n < 2) {
        if (n == 1) {
            desc[0] = {0};
            root = 0;
        }
        return;
    }
    std::vector<std::string> names(n);
    for (int i = 0; i < n; ++i) names[i] = std::to_string(i);
    HierarchicalClustering clusterer(LinkageMethod::UPGMA);
    Tree tree = clusterer.build_tree(dist, names);
    auto rt = tree.get_root();
    if (!rt) return;
    root = rt->node_id;
    // DFS fill
    std::vector<std::shared_ptr<TreeNode>> stack{rt};
    std::set<int> seen;
    while (!stack.empty()) {
        auto nd = stack.back();
        stack.pop_back();
        if (!nd || nd->node_id < 0 || nd->node_id >= 2 * n - 1) continue;
        if (!seen.insert(nd->node_id).second) continue;
        desc[nd->node_id] = nd->leaf_indices;
        if (nd->left && nd->right) {
            ch[nd->node_id] = {nd->left->node_id, nd->right->node_id};
            stack.push_back(nd->left);
            stack.push_back(nd->right);
        }
    }
}

// ---------------------------------------------------------------------------
// Per-subgroup column-permutation null: permute each CpG column among non-NaN reads of
// the subgroup, recompute BERNOULLI distance, re-cluster, take the root sibling-split
// between/within ratio; return the null_pct percentile over RNULL replicates.
// ---------------------------------------------------------------------------
double PhyloLabeler::subgroup_null_pct(const Eigen::MatrixXd& meth, const std::vector<int>& leaves,
                                       uint64_t seed) const {
    const int m = static_cast<int>(leaves.size());
    const int C = static_cast<int>(meth.cols());
    Eigen::MatrixXd Psub(m, C);
    for (int i = 0; i < m; ++i) Psub.row(i) = meth.row(leaves[i]);

    std::mt19937_64 rng(seed);
    std::vector<double> ratios;
    ratios.reserve(cfg_.rnull);
    for (int rep = 0; rep < cfg_.rnull; ++rep) {
        Eigen::MatrixXd Pn = Psub;
        for (int c = 0; c < C; ++c) {
            std::vector<int> vi;
            for (int i = 0; i < m; ++i)
                if (!std::isnan(Pn(i, c))) vi.push_back(i);
            if (vi.size() > 1) {
                std::vector<double> vals(vi.size());
                for (size_t k = 0; k < vi.size(); ++k) vals[k] = Pn(vi[k], c);
                std::shuffle(vals.begin(), vals.end(), rng);
                for (size_t k = 0; k < vi.size(); ++k) Pn(vi[k], c) = vals[k];
            }
        }
        Eigen::MatrixXd Dn = bernoulli_dist(Pn);
        std::vector<std::vector<int>> ndesc;
        std::vector<std::pair<int, int>> nch;
        int nroot = -1;
        build_desc_ch(Dn, ndesc, nch, nroot, m);
        if (nroot < 0 || nch[nroot].first < 0) continue;
        double r = between_within(Dn, ndesc[nch[nroot].first], ndesc[nch[nroot].second]);
        if (r >= 0.0) ratios.push_back(r);
    }
    if (ratios.empty()) return 0.0;
    std::sort(ratios.begin(), ratios.end());
    // numpy-style linear-interpolation percentile
    double rank = cfg_.null_pct / 100.0 * (ratios.size() - 1);
    int lo = static_cast<int>(std::floor(rank));
    int hi = static_cast<int>(std::ceil(rank));
    double frac = rank - lo;
    return ratios[lo] * (1.0 - frac) + ratios[hi] * frac;
}

// ---------------------------------------------------------------------------
// Main labeling (Phase B1: coarse).
// ---------------------------------------------------------------------------
PhyloResult PhyloLabeler::label(const Eigen::MatrixXd& dist, const Eigen::MatrixXd& meth) {
    PhyloResult res;
    const int n = static_cast<int>(dist.rows());
    if (n < 2 * cfg_.min_sz) return res;  // valid=false

    std::vector<std::vector<int>> desc;
    std::vector<std::pair<int, int>> ch;
    int root = -1;
    build_desc_ch(dist, desc, ch, root, n);
    if (root < 0) return res;

    std::vector<std::string> lab(n, "");

    // descend: peel (tiny,big) caterpillar to the first node whose both children >= min_sz.
    auto descend = [&](int node, std::vector<int>& quar) -> int {
        int cur = node;
        while (ch[cur].first >= 0) {
            int c1 = ch[cur].first, c2 = ch[cur].second;
            int s1 = static_cast<int>(desc[c1].size()), s2 = static_cast<int>(desc[c2].size());
            if (std::min(s1, s2) >= cfg_.min_sz) return cur;
            int small = (s1 < s2) ? c1 : c2;
            int big = (s1 < s2) ? c2 : c1;
            quar.insert(quar.end(), desc[small].begin(), desc[small].end());
            cur = big;
        }
        return -1;  // no balanced split
    };

    // test_split: r >= SEP_MIN and r > per-subgroup null95.
    auto test_split = [&](int bnode) -> bool {
        int c1 = ch[bnode].first, c2 = ch[bnode].second;
        double r = between_within(dist, desc[c1], desc[c2]);
        if (r < 0.0 || r < cfg_.sep_min) return false;
        // distinct seed per node (deterministic) to avoid identical null streams across nodes
        uint64_t s = cfg_.seed + static_cast<uint64_t>(bnode) * 0x9E3779B97F4A7C15ull;
        return r > subgroup_null_pct(meth, desc[bnode], s);
    };

    // recursive labeling
    std::function<void(int, const std::string&)> rec = [&](int node, const std::string& label) {
        const std::vector<int>& leaves = desc[node];
        if (static_cast<int>(leaves.size()) < 2 * cfg_.min_sz) {
            for (int i : leaves) lab[i] = label;
            return;
        }
        std::vector<int> quar;
        int bnode = descend(node, quar);
        if (bnode >= 0 && test_split(bnode)) {
            for (int i : quar) lab[i] = "outlier";  // pay-for-itself: commit only on success
            int c1 = ch[bnode].first, c2 = ch[bnode].second;
            int big = (desc[c1].size() >= desc[c2].size()) ? c1 : c2;
            int small = (desc[c1].size() >= desc[c2].size()) ? c2 : c1;
            rec(big, label + "-1");
            rec(small, label + "-2");
        } else {
            for (int i : leaves) lab[i] = label;  // whole node, no quarantine erosion
        }
    };
    rec(root, "1");

    // empty -> outlier; post-filter terminal labels with < min_sz reads -> outlier
    std::map<std::string, int> cnt;
    for (auto& l : lab) {
        if (l.empty()) l = "outlier";
        if (l != "outlier") cnt[l]++;
    }
    std::set<std::string> small;
    for (auto& kv : cnt)
        if (kv.second < cfg_.min_sz) small.insert(kv.first);
    std::set<std::string> terms;
    for (auto& l : lab) {
        if (l != "outlier" && small.count(l)) l = "outlier";
        if (l != "outlier") terms.insert(l);
    }
    res.labels = lab;
    res.n_groups = static_cast<int>(terms.size());
    res.n_outlier = static_cast<int>(std::count(lab.begin(), lab.end(), std::string("outlier")));
    res.valid = true;
    return res;
}

}  // namespace InterSubMod
