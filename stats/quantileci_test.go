// Copyright 2020 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import (
	"fmt"
	"math"
	"testing"
)

func TestQuantileCI(t *testing.T) {
	var res QuantileCIResult
	check := func(wlo, whi int, wactual float64, wambig bool) {
		t.Helper()
		if wlo != res.LoOrder || whi != res.HiOrder || !aeq(wactual, res.Confidence) || wambig != res.Ambiguous {
			t.Errorf("want [%v,%v]@%v/%v, got [%v,%v]@%v/%v",
				wlo, whi, wactual, wambig,
				res.LoOrder, res.HiOrder, res.Confidence, res.Ambiguous)
		}
	}
	eq := func(a, b float64) bool {
		return a == b ||
			math.IsInf(a, 1) && math.IsInf(b, 1) ||
			math.IsInf(a, -1) && math.IsInf(b, -1)
	}
	checkSample := func(wlo, whi float64) {
		t.Helper()
		var s Sample
		for i := 1; i <= res.N; i++ {
			s.Xs = append(s.Xs, float64(i))
		}
		s.Sorted = true
		lo, hi := res.FromSample(s)
		if !eq(wlo, lo) || !eq(whi, hi) {
			t.Errorf("want [%v,%v], got [%v,%v]", wlo, whi, lo, hi)
		}
	}

	binomBuckets := func(n int, p float64) []float64 {
		t.Helper()
		dist := BinomialDist{N: n, P: p}
		bs := make([]float64, n+1)
		t.Logf("B(%d,%v):", n, p)
		for i := range bs {
			bs[i] = dist.PMF(float64(i))
			t.Logf("  %d | %v", i, bs[i])
		}
		return bs
	}
	normBuckets := func(n int, p float64) []float64 {
		t.Helper()
		norm := BinomialDist{N: n, P: p}.NormalApprox()
		bs := make([]float64, n+1)
		t.Logf("normal approximation to B(%d,%v):", n, p)
		for i := range bs {
			bs[i] = norm.CDF(float64(i)+0.5) - norm.CDF(float64(i)-0.5)
			t.Logf("  %d | %v", i, bs[i])
		}
		return bs
	}

	// Confidence is so low that it has to fall directly around
	// the quantile.
	binomBuckets(4, 0.5) // Just for logging
	res = QuantileCI(4, 0.5, 0.001)
	check(2, 3, 0.375, false)
	checkSample(2, 3)
	res = QuantileCI(4, 0.25, 0.001)
	check(1, 2, 0.421875, false)
	checkSample(1, 2)
	// Quantile near 0.
	res = QuantileCI(4, 0, 0.001)
	check(0, 1, 1, false)
	checkSample(-inf, 1)
	res = QuantileCI(4, 0.0001, 0.001)
	check(0, 1, binomBuckets(4, 0.0001)[0], false)
	// Quantile near 1.
	res = QuantileCI(4, 1, 0.001)
	check(4, 5, 1, false)
	checkSample(4, inf)
	res = QuantileCI(4, 0.999, 0.001)
	check(4, 5, binomBuckets(4, 0.999)[4], false)
	// Confidence is exactly the PMF.
	res = QuantileCI(4, 0.5, 0.375)
	check(2, 3, 0.375, false)
	// And just beyond the PMF. This should be left-biased.
	res = QuantileCI(4, 0.5, 0.3750001)
	check(1, 3, 0.375+0.25, true)
	// Confidence is 1 or nearly 1.
	res = QuantileCI(4, 0.5, 1)
	check(0, 5, 1, false)
	res = QuantileCI(4, 0.5, 0.99)
	check(0, 5, 1, false)
	// Confidence is enough to trim one bucket. This should be
	// left-biased.
	res = QuantileCI(4, 0.5, 0.99-0.0625)
	check(0, 4, 0.375+2*0.25+0.0625, true)

	// Odd sample size with very low confidence. This should be
	// left-biased.
	binomBuckets(5, 0.5) // Just for logging
	res = QuantileCI(5, 0.5, 0.001)
	check(2, 3, 0.3125, true)
	// Confidence is exactly the PMF. This should be left-biased.
	res = QuantileCI(5, 0.5, 0.3125)
	check(2, 3, 0.3125, true)
	// And just beyond the PMF.
	res = QuantileCI(5, 0.5, 0.3125001)
	check(2, 4, 0.3125*2, false)
	// Confidence is 1 or nearly 1.
	res = QuantileCI(5, 0.5, 1)
	check(0, 6, 1, false)
	res = QuantileCI(5, 0.5, 0.99)
	check(0, 6, 1, false)
	// Confidence trims one bucket.
	res = QuantileCI(5, 0.5, 0.99-0.03125)
	check(0, 5, 1-0.03125, true)

	// Test normal approximation with even sample size.
	defer func(x int) { quantileCIApproxThreshold = x }(quantileCIApproxThreshold)
	quantileCIApproxThreshold = 0
	n := normBuckets(4, 0.5)
	// Low confidence directly around the quantile.
	res = QuantileCI(4, 0.5, 0.001)
	check(2, 3, n[2], false)
	// Confidence exactly equal to the center band.
	res = QuantileCI(4, 0.5, n[2])
	check(2, 3, n[2], false)
	// And just above. This should be left-biased.
	res = QuantileCI(4, 0.5, n[2]+0.00001)
	check(1, 3, n[1]+n[2], true)
	// Confidence is 1.
	res = QuantileCI(4, 0.5, 1)
	check(0, 5, 1, false)
	// Confidence is nearly 1. Because of the approximation, we
	// have to drop fairly low before we lose a tail, so this is
	// still the full range.
	res = QuantileCI(4, 0.5, 0.99)
	check(0, 5, 1, false)
	// Confidence is low enough to lose the right-most band. This
	// should be left-biased.
	res = QuantileCI(4, 0.5, 0.90)
	check(0, 4, n[0]+n[1]+n[2]+n[3], true)

	// Test normal approximation with odd sample size.
	n = normBuckets(5, 0.5)
	// Low confidence directly around the quantile. Left-biased.
	res = QuantileCI(5, 0.5, 0.001)
	check(2, 3, n[2], true)
	// Confidence exactly equal to the mode band. Left-biased.
	res = QuantileCI(5, 0.5, n[2])
	check(2, 3, n[2], true)
	// And just above. Symmetric.
	res = QuantileCI(5, 0.5, n[2]+0.00001)
	check(2, 4, n[2]+n[3], false)

	// Test normal approximation degenerate cases.
	res = QuantileCI(5, 0, 0.95) // 0%ile
	check(0, 1, 1, false)
	res = QuantileCI(5, 0.001, 0.95)
	check(0, 1, 1, false)
	res = QuantileCI(5, 1, 0.95) // 100%ile
	check(5, 6, 1, false)
	res = QuantileCI(5, 0.999, 0.95)
	check(5, 6, 1, false)
}

func BenchmarkQuantileCI(b *testing.B) {
	defer func(x int) { quantileCIApproxThreshold = x }(quantileCIApproxThreshold)
	for n := 5; n <= 100; n += 5 {
		for _, approx := range []bool{false, true} {
			if approx {
				quantileCIApproxThreshold = 0
			} else {
				quantileCIApproxThreshold = 1000
			}

			b.Run(fmt.Sprintf("n=%d/approx=%v", n, approx), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					QuantileCI(n, 0.5, 0.95)
				}
			})
		}
	}
}
