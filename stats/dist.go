// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

// A DistCommon is a statistical distribution. DistCommon is a base
// interface provided by both continuous and discrete distributions.
type DistCommon interface {
	// CDF returns the cumulative probability Pr[X <= x].
	//
	// For continuous distributions, the CDF is the integral of
	// the PDF from -inf to x.
	//
	// For discrete distributions, the CDF is the sum of the PMF
	// at all defined points from -inf to x, inclusive. Note that
	// the CDF of a discrete distribution is defined for the whole
	// real line (unlike the PMF) but has discontinuities where
	// the PMF is non-zero.
	//
	// The CDF is a monotonically increasing function and has a
	// domain of all real numbers. If the distribution has bounded
	// support, it has a range of [0, 1]; otherwise it has a range
	// of (0, 1). Finally, CDF(-inf)==0 and CDF(inf)==1.
	CDF(x float64) float64

	// InvCDF returns the inverse of the CDF for y. That is,
	// InvCDF(CDF(x)) = x. The value of y must be in [0, 1].
	InvCDF(y float64) float64

	// Bounds returns reasonable bounds for this distribution's
	// PDF/PMF and CDF. The total weight outside of these bounds
	// should be approximately 0.
	//
	// For a discrete distribution, both bounds are integer
	// multiples of Step().
	//
	// If this distribution has finite support, it returns exact
	// bounds l, h such that CDF(l')=0 for all l' < l and
	// CDF(h')=1 for all h' >= h.
	Bounds() (float64, float64)
}

// A Dist is a continuous statistical distribution.
type Dist interface {
	DistCommon

	// PDF returns the value of the probability density function
	// of this distribution at x.
	PDF(x float64) float64
}

// A DiscreteDist is a discrete statistical distribution.
//
// Most discrete distributions are defined only at integral values of
// the random variable. However, some are defined at other intervals,
// so this interface takes a float64 value for the random variable.
// The probability mass function rounds down to the nearest defined
// point. Note that float64 values can exactly represent integer
// values between ±2**53, so this generally shouldn't be an issue for
// integer-valued distributions (likewise, for half-integer-valued
// distributions, float64 can exactly represent all values between
// ±2**52).
type DiscreteDist interface {
	DistCommon

	// PMF returns the value of the probability mass function
	// Pr[X = x'], where x' is x rounded down to the nearest
	// defined point on the distribution.
	//
	// Note for implementers: for integer-valued distributions,
	// round x using int(math.Floor(x)). Do not use int(x), since
	// that truncates toward zero (unless all x <= 0 are handled
	// the same).
	PMF(x float64) float64

	// Step returns s, where the distribution is defined for sℕ.
	Step() float64
}

// TODO: Add a Support method for finite support distributions? Or
// maybe just another return value from Bounds indicating that the
// bounds are exact?

// TODO: Plot method to return a pre-configured Plot object with
// reasonable bounds and an integral function? Have to distinguish
// PDF/CDF/InvCDF. Three methods? Argument?
//
// Doesn't have to be a method of Dist. Could be just a function that
// takes a Dist and uses Bounds.
