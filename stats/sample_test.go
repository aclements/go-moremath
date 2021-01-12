// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import (
	"math"
	"testing"
)

func TestSampleQuantile(t *testing.T) {
	s := Sample{Xs: []float64{15, 20, 35, 40, 50}}
	testFunc(t, "Quantile", s.Quantile, map[float64]float64{
		-1:  15,
		0:   15,
		.05: 15,
		.30: 19.666666666666666,
		.40: 27,
		.95: 50,
		1:   50,
		2:   50,
	})
}

func TestMeanCI(t *testing.T) {
	var xs []float64
	naneq := func(a, b float64) bool {
		return a == b || (math.IsNaN(a) && math.IsNaN(b))
	}
	check := func(conf, wmean, wlo, whi float64) {
		t.Helper()
		mean, lo, hi := MeanCI(xs, conf)
		if !(naneq(mean, wmean) && naneq(lo, wlo) && naneq(hi, whi)) {
			t.Errorf("for %v, want %v@[%v,%v], got %v@[%v,%v]", xs, wmean, wlo, whi, mean, lo, hi)
		}
	}

	xs = []float64{-8, 2, 3, 4, 5, 6}
	check(0, 2, 2, 2)
	check(0.95, 2, -3.351092806089359, 7.351092806089359)
	check(0.99, 2, -6.39357495385287, 10.39357495385287)
	check(1, 2, -inf, inf)

	xs = []float64{1}
	check(0, 1, 1, 1)
	check(0.95, 1, -inf, inf)
	check(1, 1, -inf, inf)

	xs = nil
	check(0, math.NaN(), math.NaN(), math.NaN())
	check(0.95, math.NaN(), math.NaN(), math.NaN())
	check(1, math.NaN(), math.NaN(), math.NaN())
}
