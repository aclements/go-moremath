// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graphalg

import "github.com/aclements/go-moremath/graph"

// SimplifyMulti simplifies a multigraph to an simple graph.
//
// If g is a weighted graph, each edge in the result receives the sum
// of the weights of the combined edges in g. If g is not weighted,
// each edge in g is assumed to have a weight of 1.
func SimplifyMulti(g graph.Graph) graph.Weighted {
	gw, ok := g.(graph.Weighted)
	if !ok {
		gw = graph.WeightedUnit{g}
	}

	indexes := make([]int, gw.NumNodes()+1)
	var edges []int
	var weights []float64
	edgeMap := make(map[int]int)
	for n := range indexes[1:] {
		for k := range edgeMap {
			delete(edgeMap, k)
		}
		for i, o := range gw.Out(n) {
			if idx, ok := edgeMap[o]; ok {
				// Already have an edge.
				weights[idx] += gw.OutWeight(n, i)
			} else {
				edgeMap[o] = len(edges)
				edges = append(edges, o)
				weights = append(weights, gw.OutWeight(n, i))
			}
		}
		indexes[n+1] = len(edges)
	}
	return &simplified{indexes, edges, weights}
}

type simplified struct {
	indexes []int
	edges   []int
	weights []float64
}

func (g *simplified) NumNodes() int {
	return len(g.indexes) - 1
}

func (g *simplified) Out(n int) []int {
	return g.edges[g.indexes[n]:g.indexes[n+1]]
}

func (g *simplified) OutWeight(n, e int) float64 {
	return g.weights[g.indexes[n]:g.indexes[n+1]][e]
}
