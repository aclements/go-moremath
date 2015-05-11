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
	// M and N are the sizes of the input samples.
	M, N int

	// U is the value of the Mann-Whitney U statistic for this
	// test, generalized by counting ties as 0.5.
	//
	// Given the Cartesian product of the two samples, this is the
	// number of pairs in which the value from the first sample is
	// less than the value of the second, plus 0.5 times the
	// number of pairs where the values from the two samples are
	// equal. Hence, U is always an integer multiple of 0.5 (it is
	// a whole integer if there are no ties) in the range [0, M*N].
	//
	// The value of U given here is always the smaller of the two
	// possible values of U (depending which sample is "first").
	// The other U can be calculated as M*N - U.
	//
	// There are many equivalent statistics with slightly
	// different definitions. The Wilcoxon (1945) W statistic
	// (generalized for ties) is U + (M(M+1))/2. It is also common
	// to use 2U to eliminate the half steps and Smid (1956) uses
	// M*N - 2U to additionally center the distribution.
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

// MannWhitneyUTest performs a Mann-Whitney U-test [1] of the null
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
// than MannWhitneyExactLimit. This normal approximation uses both the
// tie correction and the continuity correction.
//
// TODO: This also uses the approximation if the samples have ties.
// This should be fixed.
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
func MannWhitneyUTest(x1, x2 []float64) (*MannWhitneyUTestResult, error) {
	n1, n2 := len(x1), len(x2)
	if n1 == 0 || n2 == 0 {
		return nil, ErrSampleSize
	}

	// Compute the U statistic.
	x1 = append([]float64(nil), x1...)
	x2 = append([]float64(nil), x2...)
	sort.Float64s(x1)
	sort.Float64s(x2)
	merged, labels := labeledMerge(x1, x2)

	R1 := 0.0
	ties := false
	for i := 0; i < len(merged); {
		rank1, nx1, v1 := i+1, 0, merged[i]
		// Consume samples that tie this sample (including itself).
		for ; i < len(merged) && merged[i] == v1; i++ {
			if labels[i] == 1 {
				nx1++
			}
		}
		if i > rank1 {
			ties = true
		}
		// Assign all tied samples the average rank of the
		// samples, where merged[0] has rank 1.
		if nx1 != 0 {
			rank := float64(i+rank1) / 2
			R1 += rank * float64(nx1)
		}
	}
	U1 := R1 - float64(n1*(n1+1))/2

	// Compute the smaller of U1 and U2
	U2 := float64(n1*n2) - U1
	if U2 < U1 {
		U1, U2 = U2, U1
	}

	// TODO: Use exact distribution when there are ties.
	var p float64
	if !ties && n1 <= MannWhitneyExactLimit && n2 <= MannWhitneyExactLimit {
		// Use exact U distribution. U1 will be an integer.
		p = UDist{M: n1, N: n2}.CDF(U1) * 2
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

	return &MannWhitneyUTestResult{M: n1, N: n2, U: U1, P: p}, nil
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
