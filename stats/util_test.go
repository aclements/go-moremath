// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import "math"

func aeq(expect, got float64) bool {
	return math.Abs(expect-got) < 0.00001
}
