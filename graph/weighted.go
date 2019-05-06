// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graph

// Weighted represented a weighted directed graph.
type Weighted interface {
	Graph

	// OutWeight returns the weight of the e'th edge out from node
	// i. e must be in the range [0, len(Out(i))).
	OutWeight(i, e int) float64
}

// WeightedUnit wraps a graph as a weighted graph where all edges have
// weight 1.
type WeightedUnit struct {
	Graph
}

// OutWeight returns 1.
func (w WeightedUnit) OutWeight(i, e int) float64 {
	return 1
}
