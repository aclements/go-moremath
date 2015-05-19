// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import "testing"

func TestMannWhitneyUTest(t *testing.T) {
	check := func(want, got *MannWhitneyUTestResult) {
		if want.N1 != got.N1 || want.N2 != got.N2 ||
			!aeq(want.U, got.U) ||
			want.AltHypothesis != got.AltHypothesis ||
			!aeq(want.P, got.P) {
			t.Errorf("want %+v, got %+v", want, got)
		}
	}

	s1 := []float64{2, 1, 3, 5}
	s2 := []float64{12, 11, 13, 15}
	s3 := []float64{0, 4, 6, 7} // Interleaved with s1, but no ties
	s4 := []float64{2, 2, 2, 2}
	s5 := []float64{1, 1, 1, 1, 1}

	// Small sample, no ties

	r, _ := MannWhitneyUTest(s1, s2, LocationDiffers)
	check(&MannWhitneyUTestResult{4, 4, 0, LocationDiffers, 0.028571428571428577}, r)

	r, _ = MannWhitneyUTest(s2, s1, LocationDiffers)
	check(&MannWhitneyUTestResult{4, 4, 16, LocationDiffers, 0.028571428571428577}, r)

	r, _ = MannWhitneyUTest(s1, s3, LocationDiffers)
	check(&MannWhitneyUTestResult{4, 4, 5, LocationDiffers, 0.485714285714285770}, r)

	// Small sample, ties
	// TODO: Check these against some other implementation.

	r, _ = MannWhitneyUTest(s1, s1, LocationDiffers)
	check(&MannWhitneyUTestResult{4, 4, 8, LocationDiffers, 1}, r)

	r, _ = MannWhitneyUTest(s1, s4, LocationDiffers)
	check(&MannWhitneyUTestResult{4, 4, 10, LocationDiffers, 0.7142857142857143}, r)

	r, _ = MannWhitneyUTest(s1, s5, LocationDiffers)
	check(&MannWhitneyUTestResult{4, 5, 17.5, LocationDiffers, 0}, r)

	r, err := MannWhitneyUTest(s4, s4, LocationDiffers)
	if err != ErrSamplesEqual {
		t.Errorf("want ErrSamplesEqual, got %+v, %+v", r, err)
	}

	// Large samples.
	l1 := make([]float64, 500)
	for i := range l1 {
		l1[i] = float64(i * 2)
	}
	l2 := make([]float64, 600)
	for i := range l2 {
		l2[i] = float64(i*2 - 41)
	}
	l3 := append([]float64{}, l2...)
	for i := 0; i < 30; i++ {
		l3[i] = l1[i]
	}
	// For comparing with R's wilcox.test:
	// l1 <- seq(0, 499)*2
	// l2 <- seq(0,599)*2-41
	// l3 <- l2; for (i in 1:30) { l3[i] = l1[i] }

	r, _ = MannWhitneyUTest(l1, l2, LocationDiffers)
	check(&MannWhitneyUTestResult{N1: 500, N2: 600, U: 135250, AltHypothesis: LocationDiffers, P: 0.0049335360814172224}, r)

	r, _ = MannWhitneyUTest(l1, l1, LocationDiffers)
	check(&MannWhitneyUTestResult{N1: 500, N2: 500, U: 125000, AltHypothesis: LocationDiffers, P: 1}, r)

	r, _ = MannWhitneyUTest(l1, l3, LocationDiffers)
	check(&MannWhitneyUTestResult{N1: 500, N2: 600, U: 134845, AltHypothesis: LocationDiffers, P: 0.0038703814239617884}, r)
}
