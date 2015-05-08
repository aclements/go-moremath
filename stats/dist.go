// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

// A Dist is a continuous statistical distribution.
type Dist interface {
	// PDF returns the value of the probability density function
	// of this distribution at x.
	PDF(x float64) float64

	// PDFEach returns PDF(xs[i]) for each i.
	PDFEach(xs []float64) []float64

	// CDF returns the value of the cumulative distribution
	// function for this distribution at x. This is the integral
	// of the PDF from 0 to x.
	CDF(x float64) float64

	// CDFEach returns CDF(xs[i]) for each i.
	CDFEach(xs []float64) []float64

	// InvCDF returns the inverse of the CDF for y. That is,
	// InvCDF(CDF(x)) = x. The value of y must be in [0, 1].
	InvCDF(y float64) float64

	// InvCDFEach returns InvCDF(ys[i]) for each i.
	InvCDFEach(ys []float64) []float64

	// Bounds returns reasonable bounds for this distribution's
	// PDF and CDF. The total weight outside of these bounds
	// should be approximately 0.
	Bounds() (float64, float64)
}

// TODO: Add a Support method for finite support distributions?

// TODO: Plot method to return a pre-configured Plot object with
// reasonable bounds and an integral function? Have to distinguish
// PDF/CDF/InvCDF. Three methods? Argument?
//
// Doesn't have to be a method of Dist. Could be just a function that
// takes a Dist and uses Bounds.
