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

func (d UDist) newMemo(U int) [][][]float64 {
	base := make([]float64, U+1)
	for i := 0; i < U+1; i++ {
		base[i] = -1
	}
	// TODO: We only use half of this table, since n <= m.
	memo := make([][][]float64, d.N+1)
	for n := 0; n < d.N+1; n++ {
		memo1 := make([][]float64, d.M+1)
		memo[n] = memo1
		for m := 0; m < d.M+1; m++ {
			memo1[m] = append([]float64(nil), base...)
		}
	}
	return memo
}

func (d UDist) pdf(memo [][][]float64, U int) float64 {
	// TODO: I'm almost certain this can be done with dynamic
	// programming, which would eliminate the huge call/stack
	// overhead and cut down the table size by a factor of max(n,
	// m).
	var rec func(n, m, U int) float64
	rec = func(n, m, U int) float64 {
		if U < 0 {
			return 0
		}
		if m > n {
			// rec(m, n, U) = rec(n, m, U), so keep m and
			// n in order in order to reuse more
			// computation.
			m, n = n, m
		}

		memop := &memo[n][m][U]
		if *memop != -1 {
			return *memop
		}

		var v float64
		if m == 0 {
			if U != 0 {
				v = 0
			} else {
				v = 1 / float64(choose(m+n, n))
			}
		} else {
			// Note that there's a typo in the Mann-Whitney paper.
			// This recurrence is written in terms of U-M. It's
			// supposed to be U-m.
			v = (float64(n)*rec(n-1, m, U-m) + float64(m)*rec(n, m-1, U)) / float64(n+m)
		}
		*memop = v
		return v
	}
	return rec(d.N, d.M, U)
}

func (d UDist) PMF(U float64) float64 {
	Ui := int(math.Floor(U))
	if Ui < 0 || Ui > d.M*d.N {
		return 0
	}
	memo := d.newMemo(Ui)
	return d.pdf(memo, Ui)
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
	memo := d.newMemo(Ui)
	p := 0.0
	for i := 0; i <= Ui; i++ {
		p += d.pdf(memo, i)
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
