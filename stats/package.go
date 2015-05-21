// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// stats is a grab bag of statistical routines.
package stats // import "github.com/aclements/go-moremath/stats"

import (
	"errors"
	"math"
)

var inf = math.Inf(1)
var nan = math.NaN()

// TODO: Put all errors in the same place and maybe unify them.

var (
	ErrSamplesEqual = errors.New("all samples are equal")
)
