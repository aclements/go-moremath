// Copyright 2020 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import (
	"fmt"
	"math"
	"testing"
)

func TestBinomialDist(t *testing.T) {
	dist := BinomialDist{N: 5, P: 0.2}
	testFunc(t, fmt.Sprintf("%+v.PMF", dist), dist.PMF,
		map[float64]float64{
			-1000: 0,
			-1:    0,
			0:     0.32768,
			1:     0.4096,
			2:     0.2048,
			3:     0.0512,
			4:     0.0064,
			5:     math.Pow(dist.P, 5),
			6:     0,
			1000:  0,
		})
	testDiscreteCDF(t, fmt.Sprintf("%+v.CDF", dist), dist)

	dist = BinomialDist{N: 30, P: 0.5}
	norm := dist.NormalApprox()
	for k := 10; k <= 20; k++ {
		b := dist.PMF(float64(k))
		n := norm.CDF(float64(k)+0.5) - norm.CDF(float64(k)-0.5)

		// The normal approximation isn't actually very close,
		// even with high N and P near 0.5, so we only check
		// the center of the distribution and we're pretty
		// lax.
		err := math.Abs(b/n - 1)
		if err > 0.01 {
			t.Errorf("want %v ≅ %v at %d", b, n, k)
		}
	}
}
