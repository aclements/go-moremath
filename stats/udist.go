// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import (
	"fmt"
	"math"
)

// A UDist is the discrete probability distribution of the
// Mann-Whitney U statistic for a pair of samples of sizes N1 and N2.
//
// The details of computing this distribution with no ties can be
// found in Mann, Henry B.; Whitney, Donald R. (1947). "On a Test of
// Whether one of Two Random Variables is Stochastically Larger than
// the Other". Annals of Mathematical Statistics 18 (1): 50–60.
// Computing this distribution in the presence of ties is described in
// Klotz, J. H. (1966). "The Wilcoxon, Ties, and the Computer".
// Journal of the American Statistical Association 61 (315): 772-787
// and Cheung, Ying Kuen; Klotz, Jerome H. (1997). "The Mann Whitney
// Wilcoxon Distribution Using Linked Lists". Statistica Sinica 7:
// 805-813 (the former paper contains details that are glossed over in
// the latter paper but has mathematical typesetting issues, so it's
// easiest to get the context from the former paper and the details
// from the latter).
type UDist struct {
	N1, N2 int

	// T is the count of the number of ties at each rank in the
	// input distributions. T may be nil, in which case it is
	// assumed there are no ties (which is equivalent to an M+N
	// slice of 1s). It must be the case that Sum(T) == M+N.
	T []int
}

// hasTies returns true if d has any tied samples.
func (d UDist) hasTies() bool {
	for _, t := range d.T {
		if t > 1 {
			return true
		}
	}
	return false
}

// p returns the p_{d.N1,d.N2} function defined by Mann, Whitney 1947
// for values of U from 0 up to and including the U argument.
//
// This algorithm runs in Θ(N1*N2*U) = O(N1²N2²) time and is quite
// fast for small values of N1 and N2. However, it does not handle ties.
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

	N, M := d.N1, d.N2
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

// permCount returns the number of sample permutations under the tie
// vector d.T with U statistic <=, == or >= twoUthresh/2, which can be
// used to directly compute the PMF and CDF of the U distribution
// under ties. This uses the "graphical method" of Klotz (1966). It is
// (currently) very slow.
func (d UDist) permCount(twoUthresh int, cmp int) (count float64) {
	// TODO: This is exponential in len(T). Implementing the
	// direct enumeration algorithm from Klotz will help, but I
	// think it's still exponential. There's the "linked list"
	// method from Cheung-Klotz; I can't make sense of that paper,
	// but there's an implementation in appendix L of
	// http://www.stat.wisc.edu/~klotz/Book.pdf. There's also van
	// de Wiel, "The split-up algorithm: a fast symbolic method
	// for computing p-values of distribution-free statistics"
	// that sounds really promising, but I can't get my hands on
	// it. The "split-up algorithm" seems to be the technique used
	// by R's coin package.

	// Enumerate all u vectors such that  0 <= u_i <= t_i.
	u := make([]int, len(d.T))
	u[len(u)-1] = -1 // Get enumeration started.
	for {
		// Compute the next u vector.
		u[len(u)-1]++
		for i := len(u) - 1; i >= 0 && u[i] > d.T[i]; i-- {
			if i == 0 {
				// All u vectors have been enumerated.
				return
			}
			// Carry.
			u[i-1]++
			u[i] = 0
		}

		// Is this a legal u vector?
		if sumint(u) != d.N1 {
			// TODO: Implement optimized enumeration
			// method that does not construct illegal
			// candidates.
			continue
		}

		// TODO: Reuse partial U and prod computations.

		// Compute 2*U statistic for this u vector.
		twoU, vsum := 0, 0
		for i, u_i := range u {
			v_i := d.T[i] - u_i
			// U = U + vsum*u_i + u_i*v_i/2
			twoU += 2*vsum*u_i + u_i*v_i
			vsum += v_i
		}

		if cmp < 0 && twoU > twoUthresh {
			continue
		} else if cmp == 0 && twoU != twoUthresh {
			continue
		} else if cmp > 0 && twoU < twoUthresh {
			continue
		}

		// Compute Π choose(t_i, u_i). This is the number of
		// ways of permuting the input sample under u.
		prod := 1
		for i, u_i := range u {
			prod *= choose(d.T[i], u_i)
		}

		// Accumulate the permutations on this u path.
		count += float64(prod)

		if false {
			// Print a table in the form of Klotz's
			// "direct enumeration" example.
			//
			// Convert 2U = 2UQV' to UQt' used in Klotz
			// examples.
			UQt := float64(twoU)/2 + float64(d.N1*d.N1)/2
			fmt.Printf("%+v %f %-2d\n", u, UQt, prod)
		}
	}
}

func (d UDist) PMF(U float64) float64 {
	if U < 0 || U >= 0.5+float64(d.N1*d.N2) {
		return 0
	}

	if d.hasTies() {
		return d.permCount(int(2*U), 0) / float64(choose(d.N1+d.N2, d.N1))
	}

	// There are no ties. Use the fast algorithm. U must be integral.
	Ui := int(math.Floor(U))
	// TODO: Use symmetry to minimize U
	return d.p(Ui)[Ui]
}

func (d UDist) CDF(U float64) float64 {
	if U < 0 {
		return 0
	} else if U >= float64(d.N1*d.N2) {
		return 1
	}

	if d.hasTies() {
		return d.permCount(int(2*U), -1) / float64(choose(d.N1+d.N2, d.N1))
	}

	// There are no ties. Use the fast algorithm. U must be integral.
	Ui := int(math.Floor(U))
	// The distribution is symmetric around U = m * n / 2. Sum up
	// whichever tail is smaller.
	flip := Ui >= (d.N1*d.N2+1)/2
	if flip {
		Ui = d.N1*d.N2 - Ui - 1
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
	return 0.5
}

func (d UDist) Bounds() (float64, float64) {
	// TODO: More precise bounds when there are ties.
	return 0, float64(d.N1 * d.N2)
}
