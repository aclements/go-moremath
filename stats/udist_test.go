// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import (
	"fmt"
	"math"
	"testing"
)

func aeqTable(a, b [][]float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return false
		}
		for j := range a[i] {
			// "%f" precision
			if math.Abs(a[i][j]-b[i][j]) >= 0.000001 {
				return false
			}
		}
	}
	return true
}

// U distribution for N=3 up to U=5.
var udist3 = [][]float64{
	//    m=1         2         3
	{0.250000, 0.100000, 0.050000}, // U=0
	{0.500000, 0.200000, 0.100000}, // U=1
	{0.750000, 0.400000, 0.200000}, // U=2
	{1.000000, 0.600000, 0.350000}, // U=3
	{1.000000, 0.800000, 0.500000}, // U=4
	{1.000000, 0.900000, 0.650000}, // U=5
}

// U distribution for N=5 up to U=5.
var udist5 = [][]float64{
	//    m=1         2         3         4         5
	{0.166667, 0.047619, 0.017857, 0.007937, 0.003968}, // U=0
	{0.333333, 0.095238, 0.035714, 0.015873, 0.007937}, // U=1
	{0.500000, 0.190476, 0.071429, 0.031746, 0.015873}, // U=2
	{0.666667, 0.285714, 0.125000, 0.055556, 0.027778}, // U=3
	{0.833333, 0.428571, 0.196429, 0.095238, 0.047619}, // U=4
	{1.000000, 0.571429, 0.285714, 0.142857, 0.075397}, // U=5
}

func TestUDist(t *testing.T) {
	makeTable := func(n int) [][]float64 {
		out := make([][]float64, 6)
		for U := 0; U < 6; U++ {
			out[U] = make([]float64, n)
			for m := 1; m <= n; m++ {
				out[U][m-1] = UDist{N: n, M: m}.CDF(float64(U))
			}
		}
		return out
	}
	fmtTable := func(a [][]float64) string {
		out := fmt.Sprintf("%8s", "m=")
		for m := 1; m <= len(a[0]); m++ {
			out += fmt.Sprintf("%9d", m)
		}
		out += "\n"

		for U, row := range a {
			out += fmt.Sprintf("U=%-6d", U)
			for m := 1; m <= len(a[0]); m++ {
				out += fmt.Sprintf(" %f", row[m-1])
			}
			out += "\n"
		}
		return out
	}

	// Compare against tables given in Mann, Whitney (1947).
	got3 := makeTable(3)
	if !aeqTable(got3, udist3) {
		t.Errorf("For n=3, want:\n%sgot:\n%s", fmtTable(udist3), fmtTable(got3))
	}

	got5 := makeTable(5)
	if !aeqTable(got5, udist5) {
		t.Errorf("For n=5, want:\n%sgot:\n%s", fmtTable(udist5), fmtTable(got5))
	}
}

func BenchmarkUDist(b *testing.B) {
	for i := 0; i < b.N; i++ {
		// R uses the exact distribution up to N=50.
		// N*M/2=1250 is the hardest point to get the CDF for.
		UDist{N: 50, M: 50}.CDF(1250)
	}
}
