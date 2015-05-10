// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import "math"

// A UDist is the discrete probability distribution of the
// Mann-Whitney U statistic for a pair of samples of sizes M and N.
//
// The details of computing this distribution can be found in Mann,
// Henry B.; Whitney, Donald R. (1947). "On a Test of Whether one of
// Two Random Variables is Stochastically Larger than the Other".
// Annals of Mathematical Statistics 18 (1): 50–60.
type UDist struct {
	M, N int
}

// p returns the p_{d.N,d.M} function defined by Mann, Whitney 1947
// for values of U from 0 up to and including the U argument.
func (d UDist) p(U int) []float64 {
	// This is a dynamic programming implementation of the
	// recursive recurrence definition given by Mann and Whitney:
	//
	//   p_{n,m}(U) = (n * p_{n-1,m}(U-m) + m * p_{n,m-1}(U)) / (n+m)
	//   p_{n,m}(U) = 0                           if U < 0
	//   p_{0,m}(U) = p{n,0}(U) = 1 / nCr(m+n, n) if U = 0
	//                          = 0               if U > 0
	//
	// (Note that there is a typo in the original paper. The first
	// recursive application of p should be for U-m, not U-M.)
	//
	// Since p{n,m} only depends on p{n-1,m} and p{n,m-1}, we only
	// need to store one "plane" of the three dimensional space at
	// a time.
	//
	// Furthermore, p_{n,m} = p_{m,n}, so we only construct values
	// for n <= m and obtain the rest through symmetry.
	//
	// We organize the computed values of p as followed:
	//
	//       n →   N
	//     m *
	//     ↓ * *
	//       * * *
	//       * * * *
	//       * * * *
	//     M * * * *
	//
	// where each * is a slice indexed by U. The code below
	// computes these left-to-right, top-to-bottom, so it only
	// stores one row of this matrix at a time. Furthermore,
	// computing an element in a given U slice only depends on the
	// same and smaller values of U, so we can overwrite the U
	// slice we're computing in place as long as we start with the
	// largest value of U. Finally, even though the recurrence
	// depends on (n,m) above the diagonal and we use symmetry to
	// mirror those across the diagonal to (m,n), the mirrored
	// indexes are always available in the current row, so this
	// mirroring does not interfere with our ability to recycle
	// state.

	N, M := d.N, d.M
	if N > M {
		N, M = M, N
	}

	memo := make([][]float64, N+1)
	for n := range memo {
		memo[n] = make([]float64, U+1)
	}

	for m := 0; m <= M; m++ {
		// Compute p_{0,m}. This is zero except for U=0.
		memo[0][0] = 1

		// Compute the remainder of this row.
		nlim := N
		if m < nlim {
			nlim = m
		}
		for n := 1; n <= nlim; n++ {
			lp := memo[n-1] // p_{n-1,m}
			var rp []float64
			if n <= m-1 {
				rp = memo[n] // p_{n,m-1}
			} else {
				rp = memo[m-1] // p{m-1,n} and m==n
			}

			// For a given n,m, U is at most n*m.
			//
			// TODO: Actually, it's at most ⌈n*m/2⌉, but
			// then we need to use more complex symmetries
			// in the inner loop below.
			ulim := n * m
			if U < ulim {
				ulim = U
			}

			out := memo[n] // p_{n,m}
			nplusm := float64(n + m)
			for U1 := ulim; U1 >= 0; U1-- {
				l := 0.0
				if U1-m >= 0 {
					l = float64(n) * lp[U1-m]
				}
				r := float64(m) * rp[U1]
				out[U1] = (l + r) / nplusm
			}
		}
	}
	return memo[N]
}

func (d UDist) PMF(U float64) float64 {
	Ui := int(math.Floor(U))
	// TODO: Use symmetry to minimize U
	if Ui < 0 || Ui > d.M*d.N {
		return 0
	}
	return d.p(Ui)[Ui]
}

func (d UDist) CDF(U float64) float64 {
	Ui := int(math.Floor(U))
	if Ui < 0 {
		return 0
	} else if Ui >= d.M*d.N {
		return 1
	}
	// The distribution is symmetric around U = m * n / 2. Sum up
	// whichever tail is smaller.
	flip := Ui >= (d.M*d.N+1)/2
	if flip {
		Ui = d.M*d.N - Ui - 1
	}
	pdfs := d.p(Ui)
	p := 0.0
	for _, pdf := range pdfs[:Ui+1] {
		p += pdf
	}
	if flip {
		p = 1 - p
	}
	return p
}

func (d UDist) Step() float64 {
	return 1
}

func (d UDist) Bounds() (float64, float64) {
	return 0, float64(d.M * d.N)
}
