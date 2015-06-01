// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package stats

import (
	"fmt"
	"math"
)

// TODO: Consider moving this to stats/kde.  Then I could write things
// like kde.Setup{Bandwidth: kde.Scott}.FromSample(sample)

// KDE represents options for constructing a kernel density estimate.
//
// Kernel density estimation is a method for constructing an estimate
// ƒ̂(x) of a unknown distribution ƒ(x) given a sample from that
// distribution.  Unlike many techniques, kernel density estimation is
// non-parametric: in general, it doesn't assume any particular true
// distribution (note, however, that the resulting distribution
// depends deeply on the selected bandwidth, and many bandwidth
// estimation techniques assume normal reference rules).
//
// A kernel density estimate is similar to a histogram, except that it
// is a smooth probability estimate and does not require choosing a
// bin size and discretizing the data.
//
// To construct a kernel density estimate, create an instance of KDE
// and then use the From method to provide data.
//
// The default (zero) value of KDE is a reasonable default
// configuration.
type KDE struct {
	// Kernel is the kernel to use for the KDE.
	Kernel KDEKernel

	// Bandwidth is the bandwidth to use for the KDE.
	//
	// If this is zero, the bandwidth is computed from the
	// provided data using a default bandwidth estimator
	// (currently BandwidthScott).
	Bandwidth float64

	// BoundaryMethod is the boundary correction method to use for
	// the KDE. The default value is BoundaryReflect; however, the
	// default bounds are effectively +/-inf, which is equivalent
	// to performing no boundary correction.
	BoundaryMethod KDEBoundaryMethod

	// [BoundaryMin, BoundaryMax) specify a bounded support for
	// the KDE. If both are 0 (their default values), they are
	// treated as +/-inf.
	//
	// To specify a half-bounded support, set Min to math.Inf(-1)
	// or Max to math.Inf(1).
	BoundaryMin float64
	BoundaryMax float64
}

// BandwidthSilverman is a bandwidth estimator implementing
// Silverman's Rule of Thumb. It's fast, but not very robust to
// outliers as it assumes data is approximately normal.
//
// Silverman, B. W. (1986) Density Estimation.
func BandwidthSilverman(data interface {
	StdDev() float64
	Weight() float64
}) float64 {
	return 1.06 * data.StdDev() * math.Pow(data.Weight(), -1.0/5)
}

// BandwidthScott is a bandwidth estimator implementing Scott's Rule.
// This is generally robust to outliers: it chooses the minimum
// between the sample's standard deviation and an robust estimator of
// a Gaussian distribution's standard deviation.
//
// Scott, D. W. (1992) Multivariate Density Estimation: Theory,
// Practice, and Visualization.
func BandwidthScott(data interface {
	StdDev() float64
	Weight() float64
	Percentile(float64) float64
}) float64 {
	iqr := data.Percentile(0.75) - data.Percentile(0.25)
	hScale := 1.06 * math.Pow(data.Weight(), -1.0/5)
	stdDev := data.StdDev()
	if stdDev < iqr/1.349 {
		// Use Silverman's Rule of Thumb
		return hScale * stdDev
	} else {
		// Use IQR/1.349 as a robust estimator of the standard
		// deviation of a Gaussian distribution.
		return hScale * (iqr / 1.349)
	}
}

// TODO(austin) Implement bandwidth estimator from Botev, Grotowski,
// Kroese. (2010) Kernel Density Estimation via Diffusion.

// KDEKernel represents a kernel to use for a KDE.
type KDEKernel int

//go:generate stringer -type=KDEKernel

const (
	GaussianKernel KDEKernel = iota

	// DeltaKernel is a Dirac delta function.  The PDF of such a
	// KDE is not well-defined, but the CDF will represent each
	// sample as an instantaneous increase.  This kernel ignores
	// bandwidth and never requires boundary correction.
	DeltaKernel
)

// KDEBoundaryMethod represents a boundary correction method for
// constructing a KDE with bounded support.
type KDEBoundaryMethod int

//go:generate stringer -type=KDEBoundaryMethod

const (
	// BoundaryReflect reflects the density estimate at the
	// boundaries.  For example, for a KDE with support [0, inf),
	// this is equivalent to ƒ̂ᵣ(x)=ƒ̂(x)+ƒ̂(-x) for x>=0.  This is a
	// simple and fast technique, but enforces that ƒ̂ᵣ'(0)=0, so
	// it may not be applicable to all distributions.
	BoundaryReflect KDEBoundaryMethod = iota

	// boundaryNone represents no boundary correction.
	//
	// This is used internally when the bounds are -/+inf.
	boundaryNone
)

// From returns the probability density function of the kernel density
// estimate for the sample s.
func (k KDE) From(s Sample) Dist {
	if s.Weights != nil && len(s.Xs) != len(s.Weights) {
		panic("len(xs) != len(weights)")
	}

	// Compute bandwidth
	h := k.Bandwidth
	if h == 0 {
		h = BandwidthScott(s)
	}

	// Construct kernel
	kernel := kdeKernel(nil)
	switch k.Kernel {
	default:
		panic(fmt.Sprint("unknown kernel", k))
	case GaussianKernel:
		kernel = NormalDist{0, h}
	case DeltaKernel:
		kernel = DeltaDist{0}
	}

	// Normalize boundaries
	bm := k.BoundaryMethod
	min, max := k.BoundaryMin, k.BoundaryMax
	if min == 0 && max == 0 {
		min, max = math.Inf(-1), math.Inf(1)
	}
	if math.IsInf(min, -1) && math.IsInf(max, 1) {
		bm = boundaryNone
	}

	return &kdeDist{kernel, s.Xs, s.Weights, bm, min, max}
}

// TODO: For KDEs of histograms, make histograms able to create a
// weighted Sample and simply require the caller to provide a
// good bandwidth from a StreamStats.

type kdeKernel interface {
	PDFEach(xs []float64) []float64
	CDFEach(xs []float64) []float64
}

type kdeDist struct {
	kernel      kdeKernel
	xs, weights []float64
	bm          KDEBoundaryMethod
	min, max    float64 // Support bounds
}

// normalizedXs returns x - kde.xs.  Evaluating kernels shifted by
// kde.xs all at x is equivalent to evaluating one unshifted kernel at
// x - kde.xs.
func (kde *kdeDist) normalizedXs(x float64) []float64 {
	txs := make([]float64, len(kde.xs))
	for i, xi := range kde.xs {
		txs[i] = x - xi
	}
	return txs
}

func (kde *kdeDist) PDF(x float64) float64 {
	// Apply boundary
	if x < kde.min || x >= kde.max {
		return 0
	}

	y := func(x float64) float64 {
		// Shift kernel to each of kde.xs and evaluate at x
		ys := kde.kernel.PDFEach(kde.normalizedXs(x))

		// Kernel samples are weighted according to the weights of xs
		wys := Sample{Xs: ys, Weights: kde.weights}

		return wys.Sum() / wys.Weight()
	}
	switch kde.bm {
	default:
		panic("unknown boundary correction method")
	case boundaryNone:
		return y(x)
	case BoundaryReflect:
		if math.IsInf(kde.max, 1) {
			return y(x) + y(2*kde.min-x)
		} else if math.IsInf(kde.min, -1) {
			return y(x) + y(2*kde.max-x)
		} else {
			d := 2 * (kde.max - kde.min)
			w := 2 * (x - kde.min)
			return series(func(n float64) float64 {
				// Points >= x
				return y(x+n*d) + y(x+n*d-w)
			}) + series(func(n float64) float64 {
				// Points < x
				return y(x-(n+1)*d+w) + y(x-(n+1)*d)
			})
		}
	}
}

func (cdf *kdeDist) CDF(x float64) float64 {
	// Apply boundary
	if x < cdf.min {
		return 0
	} else if x >= cdf.max {
		return 1
	}

	y := func(x float64) float64 {
		// Shift kernel integral to each of cdf.xs and evaluate at x
		ys := cdf.kernel.CDFEach(cdf.normalizedXs(x))

		// Kernel samples are weighted according to the weights of xs
		wys := Sample{Xs: ys, Weights: cdf.weights}

		return wys.Sum() / wys.Weight()
	}
	switch cdf.bm {
	default:
		panic("unknown boundary correction method")
	case boundaryNone:
		return y(x)
	case BoundaryReflect:
		if math.IsInf(cdf.max, 1) {
			return y(x) - y(2*cdf.min-x)
		} else if math.IsInf(cdf.min, -1) {
			return y(x) + (1 - y(2*cdf.max-x))
		} else {
			d := 2 * (cdf.max - cdf.min)
			w := 2 * (x - cdf.min)
			return series(func(n float64) float64 {
				// Windows >= x-w
				return y(x+n*d) - y(x+n*d-w)
			}) + series(func(n float64) float64 {
				// Windows < x-w
				return y(x-(n+1)*d) - y(x-(n+1)*d-w)
			})
		}
	}
}

func (cdf *kdeDist) Bounds() (low float64, high float64) {
	// TODO(austin) If this KDE came from a histogram, we'd better
	// not sample at a significantly higher rate than the
	// histogram.  Maybe we want to just return the bounds of the
	// histogram?

	// TODO(austin) It would be nice if this could be instructed
	// to include all original data points, even if they are in
	// the tail.  Probably that should just be up to the caller to
	// pass an axis derived from the bounds of the original data.

	// Use the lowest and highest samples as starting points
	lowX, highX := Sample{Xs: cdf.xs, Weights: cdf.weights}.Bounds()
	if lowX == highX {
		lowX -= 1
		highX += 1
	}

	// Find the end points that contain 99% of the CDF's weight.
	// Since bisect requires that the root be bracketed, start by
	// expanding our range if necessary.  TODO(austin) This can
	// definitely be done faster.
	const (
		lowY      = 0.005
		highY     = 0.995
		tolerance = 0.001
	)
	for cdf.CDF(lowX) > lowY {
		lowX -= highX - lowX
	}
	for cdf.CDF(highX) < highY {
		highX += highX - lowX
	}
	// Explicitly accept discontinuities, since we may be using a
	// discontiguous kernel.
	low, _ = bisect(func(x float64) float64 { return cdf.CDF(x) - lowY }, lowX, highX, tolerance)
	high, _ = bisect(func(x float64) float64 { return cdf.CDF(x) - highY }, lowX, highX, tolerance)

	// Expand width by 20% to give some margins
	width := high - low
	low, high = low-0.1*width, high+0.1*width

	// Limit to bounds
	low, high = math.Max(low, cdf.min), math.Min(high, cdf.max)

	return
}
