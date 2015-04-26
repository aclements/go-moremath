// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

func aeq(expect, got float64) bool {
	if expect < 0 && got < 0 {
		expect, got = -expect, -got
	}
	return expect*0.99999999 <= got && got*0.99999999 <= expect
}
