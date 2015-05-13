// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import (
	"math"
	"sort"
)

// A MannWhitneyUTestResult is the result of a Mann-Whitney U-test.
type MannWhitneyUTestResult struct {
	// N1 and N2 are the sizes of the input samples.
	N1, N2 int

	// U is the value of the Mann-Whitney U statistic for this
	// test, generalized by counting ties as 0.5.
	//
	// Given the Cartesian product of the two samples, this is the
	// number of pairs in which the value from the first sample is
	// less than the value of the second, plus 0.5 times the
	// number of pairs where the values from the two samples are
	// equal. Hence, U is always an integer multiple of 0.5 (it is
	// a whole integer if there are no ties) in the range [0, N1*N2].
	//
	// The value of U given here is always the smaller of the two
	// possible values of U (depending which sample is "first").
	// The other U can be calculated as N1*N2 - U.
	//
	// There are many equivalent statistics with slightly
	// different definitions. The Wilcoxon (1945) W statistic
	// (generalized for ties) is U + (N1(N1+1))/2. It is also
	// common to use 2U to eliminate the half steps and Smid
	// (1956) uses N1*N2 - 2U to additionally center the
	// distribution.
	//
	// TODO: Maybe this should always be the U statistic of the
	// first sample. Currently you can't back out which of the two
	// samples this is the U for.
	U float64

	// P is the two-tailed p-value of the Mann-Whitney test.
	//
	// TODO: One-tailed? It's not just P/2. Should I add
	// one-tailed as another field, or take it as an argument?
	P float64
}

// MannWhitneyExactLimit gives the largest sample size for which the
// exact U distribution will be used for the Mann-Whitney U-test.
//
// Using the exact distribution is necessary for small sample sizes
// because the distribution is highly irregular. However, computing
// the distribution for large sample sizes is both computationally
// expensive and unnecessary because it quickly approaches a normal
// approximation. Computing the distribution for two 50 value samples
// takes a few milliseconds on a 2014 laptop.
var MannWhitneyExactLimit = 50

// MannWhitneyTiesExactLimit gives the largest sample size for which
// the exact U distribution will be used for the Mann-Whitney U-test
// in the presence of ties.
//
// Computing this distribution is (currently) much more expensive than
// computing the distribution without ties, so this is set much lower.
// Computing this distribution for two 9 value samples takes a few
// tens of milliseconds on a 2014 laptop.
var MannWhitneyTiesExactLimit = 9

// MannWhitneyUTest performs a Mann-Whitney U-test [1,2] of the null
// hypothesis that two samples come from the same population against
// the alternative hypothesis that one sample tends to have larger or
// smaller values than the other.
//
// This is similar to a t-test, but unlike the t-test, the
// Mann-Whitney U-test is non-parametric (it does not assume a normal
// distribution). It has very slightly lower efficiency than the
// t-test on normal distributions.
//
// Computing the exact U distribution is expensive for large sample
// sizes, so this uses a normal approximation for sample sizes larger
// than MannWhitneyExactLimit if there are no ties or
// MannWhitneyTiesExactLimit if there are ties. This normal
// approximation uses both the tie correction and the continuity
// correction.
//
// This can fail with ErrSampleSize if either sample is empty or
// ErrSamplesEqual if all sample values are equal.
//
// This is also known as a Mann-Whitney-Wilcoxon test and is
// equivalent to the Wilcoxon rank-sum test, though the Wilcoxon
// rank-sum test differs in nomenclature.
//
// [1] Mann, Henry B.; Whitney, Donald R. (1947). "On a Test of
// Whether one of Two Random Variables is Stochastically Larger than
// the Other". Annals of Mathematical Statistics 18 (1): 50–60.
//
// [2] Klotz, J. H. (1966). "The Wilcoxon, Ties, and the Computer".
// Journal of the American Statistical Association 61 (315): 772-787.
func MannWhitneyUTest(x1, x2 []float64) (*MannWhitneyUTestResult, error) {
	n1, n2 := len(x1), len(x2)
	if n1 == 0 || n2 == 0 {
		return nil, ErrSampleSize
	}

	// Compute the U statistic and tie vector T.
	x1 = append([]float64(nil), x1...)
	x2 = append([]float64(nil), x2...)
	sort.Float64s(x1)
	sort.Float64s(x2)
	merged, labels := labeledMerge(x1, x2)

	R1 := 0.0
	T, hasTies := []int{}, false
	for i := 0; i < len(merged); {
		rank1, nx1, v1 := i+1, 0, merged[i]
		// Consume samples that tie this sample (including itself).
		for ; i < len(merged) && merged[i] == v1; i++ {
			if labels[i] == 1 {
				nx1++
			}
		}
		// Assign all tied samples the average rank of the
		// samples, where merged[0] has rank 1.
		if nx1 != 0 {
			rank := float64(i+rank1) / 2
			R1 += rank * float64(nx1)
		}
		T = append(T, i-rank1+1)
		if i > rank1 {
			hasTies = true
		}
	}
	U1 := R1 - float64(n1*(n1+1))/2

	// Compute the smaller of U1 and U2
	U2 := float64(n1*n2) - U1
	if U2 < U1 {
		U1, U2 = U2, U1
	}

	var p float64
	if !hasTies && n1 <= MannWhitneyExactLimit && n2 <= MannWhitneyExactLimit ||
		hasTies && n1 <= MannWhitneyTiesExactLimit && n2 <= MannWhitneyTiesExactLimit {
		// Use exact U distribution. U1 will be an integer.
		if len(T) == 1 {
			// All values are equal. Test is meaningless.
			return nil, ErrSamplesEqual
		} else if U1 == U2 {
			// The distribution is symmetric about U1.
			// Since the distribution is discrete, the CDF
			// is discontinuous and if simply double
			// CDF(U1), we'll double count the
			// (non-infinitesimal) probability mass at U1.
			// What we want is just the integral of the
			// whole CDF, which is 1.
			p = 1
		} else {
			p = UDist{M: n1, N: n2, T: T}.CDF(U1) * 2
		}
	} else {
		// Use normal approximation (with tie and continuity
		// correction).
		t := tieCorrection(merged)
		N := float64(n1 + n2)
		μ_U := float64(n1*n2) / 2
		σ_U := math.Sqrt(float64(n1*n2) * ((N + 1) - t/(N*(N-1))) / 12)
		numer := U1 - μ_U
		numer -= sign(numer) * 0.5 // Continuity correction
		if σ_U == 0 {
			return nil, ErrSamplesEqual
		} else {
			z := numer / σ_U
			p = 2 * math.Min(StdNormal.CDF(z), 1-StdNormal.CDF(z))
		}
	}

	return &MannWhitneyUTestResult{N1: n1, N2: n2, U: U1, P: p}, nil
}

// labeledMerge merges sorted lists x1 and x2 into sorted list merged.
// labels[i] is 1 or 2 depending on whether merged[i] is a value from
// x1 or x2, respectively.
func labeledMerge(x1, x2 []float64) (merged []float64, labels []byte) {
	merged = make([]float64, len(x1)+len(x2))
	labels = make([]byte, len(x1)+len(x2))

	i, j, o := 0, 0, 0
	for i < len(x1) && j < len(x2) {
		if x1[i] < x2[j] {
			merged[o] = x1[i]
			labels[o] = 1
			i++
		} else {
			merged[o] = x2[j]
			labels[o] = 2
			j++
		}
		o++
	}
	for ; i < len(x1); i++ {
		merged[o] = x1[i]
		labels[o] = 1
		o++
	}
	for ; j < len(x2); j++ {
		merged[o] = x2[j]
		labels[o] = 2
		o++
	}
	return
}

// tieCorrection computes the tie correction factor Σ_j (t_j³ - t_j)
// where t_j is the number of ties in the j'th rank.
func tieCorrection(xs []float64) float64 {
	t := 0
	for i := 0; i < len(xs); {
		i1, v1 := i, xs[i]
		for ; i < len(xs) && xs[i] == v1; i++ {
		}
		run := i - i1
		t += run*run*run - run
	}
	return float64(t)
}
