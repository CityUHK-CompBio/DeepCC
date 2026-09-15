#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

static double sparseES(const double* sortedVals,
                       const int* hitRanks,
                       int m,
                       int N,
                       double exponent) {
  if (m == 0 || N <= 0) return 0.0;
  double nr = 0.0;
  for (int i = 0; i < m; i++) nr += std::pow(std::abs(sortedVals[hitRanks[i]]), exponent);
  if (nr == 0.0 || std::isnan(nr) || std::isinf(nr)) return 0.0;

  double step = 1.0 / (double)(N - m);
  double hi = 0.0, lo = 0.0;
  double cum = 0.0;

  // Maxima: value after each hit
  for (int i = 0; i < m; i++) {
    cum += std::pow(std::abs(sortedVals[hitRanks[i]]), exponent) / nr;
    double cand = cum - (double)(hitRanks[i] - i) * step;
    if (cand > hi) hi = cand;
    // Minima: end of the non-member gap between hit i and hit i+1
    if (i < m - 1 && hitRanks[i + 1] > hitRanks[i] + 1) {
      double lcand = cum - (double)(hitRanks[i + 1] - i - 1) * step;
      if (lcand < lo) lo = lcand;
    }
  }
  // Before first hit
  if (hitRanks[0] > 0) {
    double lcand = -(double)hitRanks[0] * step;
    if (lcand < lo) lo = lcand;
  }
  // After last hit
  {
    double lcand = cum - (double)(N - m) * step;
    if (lcand < lo) lo = lcand;
  }

  if (std::abs(hi) > std::abs(lo)) return hi;
  if (std::abs(lo) > std::abs(hi)) return lo;
  return lo;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix calcEnrichmentScoreBatchCPP(Rcpp::NumericMatrix expr,
                                                Rcpp::IntegerVector setStarts,
                                                Rcpp::IntegerVector setMembers,
                                                int nthreads = 1,
                                                double exponent = 1.0) {
  int S = expr.nrow();
  int G = expr.ncol();
  int P = setStarts.size() - 1;
  if (P <= 0) return Rcpp::NumericMatrix(S, 0);

  std::vector<int> ss(setStarts.begin(), setStarts.end());
  std::vector<int> sm(setMembers.begin(), setMembers.end());
  Rcpp::NumericMatrix out(S, P);
  for (int s = 0; s < S; s++) {
    std::vector<int> sortIdx(G);
    std::vector<int> rankOf(G, -1);
    std::vector<double> sv(G);
    std::vector<int> hitBuf;

    int vc = 0;
    for (int j = 0; j < G; j++) {
        double v = expr(s, j);
      if (!R_IsNA(v) && !std::isnan(v)) {
        sortIdx[vc] = j;
        vc++;
      }
    }
    if (vc == 0) {
        for (int p = 0; p < P; p++) out(s, p) = 0.0;
      continue;
    }

    std::stable_sort(sortIdx.begin(), sortIdx.begin() + vc,
                       [&](int a, int b) { return expr(s, a) > expr(s, b); });
      for (int j = 0; j < vc; j++) sv[j] = expr(s, sortIdx[j]);
    for (int j = 0; j < vc; j++) rankOf[sortIdx[j]] = j;

      for (int p = 0; p < P; p++) {
        hitBuf.clear();
      for (int k = ss[p]; k < ss[p + 1]; k++) {
        int gi = sm[k];
        if (gi >= 0 && gi < G && rankOf[gi] >= 0) hitBuf.push_back(rankOf[gi]);
      }
      if (hitBuf.size() > 1) std::sort(hitBuf.begin(), hitBuf.end());
      out(s, p) = sparseES(sv.data(), hitBuf.data(),
                                                 (int)hitBuf.size(), vc, exponent);
    }
  }

  return out;
}
