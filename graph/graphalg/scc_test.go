// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graphalg

import (
	"reflect"
	"sort"
	"testing"

	"github.com/aclements/go-moremath/graph"
	"github.com/aclements/go-moremath/graph/graphout"
)

type sccTest struct {
	g          graph.Graph
	components [][]int // Component -> Sub-node IDs (sorted)
	edges      graph.Graph
}

// SCC example from CLRS.
var clrsSCC = sccTest{
	graph.IntGraph{
		0: {1},
		1: {2, 4, 5},
		2: {3, 6},
		3: {2, 7},
		4: {0, 5},
		5: {6},
		6: {5, 7},
		7: {7},
	},
	[][]int{
		0: {7},
		1: {5, 6},
		2: {2, 3},
		3: {0, 1, 4},
	},
	graph.IntGraph{
		0: {},
		1: {0},
		2: {0, 1},
		3: {1, 2},
	},
}

// SCC example from Sedgewick, Algorithms in C, Part 5, Third Edition, p. 199.
var sedgewickSCC = sccTest{
	graph.IntGraph{
		0:  {2},
		1:  {0},
		2:  {3, 4},
		3:  {2, 4},
		4:  {5, 6},
		5:  {0, 3},
		6:  {0, 7},
		7:  {8},
		8:  {7},
		9:  {6, 8, 12},
		10: {9},
		11: {4, 9},
		12: {10, 11},
	},
	[][]int{
		3: {9, 10, 11, 12},
		2: {1},
		1: {0, 2, 3, 4, 5, 6},
		0: {7, 8},
	},
	graph.IntGraph{
		3: {0, 1},
		2: {1},
		1: {0},
		0: {},
	},
}

// SCC example from
// https://algs4.cs.princeton.edu/lectures/42DirectedGraphs-2x2.pdf
//
// This is very similar to the Sedgewick graph, but not the same.
// (Maybe this is from the fourth edition?)
var sedgewick2SCC = sccTest{
	graph.IntGraph{
		0:  {1, 5},
		1:  {},
		2:  {0, 3},
		3:  {2, 5},
		4:  {2, 3},
		5:  {4},
		6:  {0, 4, 8, 9},
		7:  {6, 9},
		8:  {6},
		9:  {10, 11},
		10: {12},
		11: {4, 12},
		12: {9},
	},
	[][]int{
		0: {1},
		1: {0, 2, 3, 4, 5},
		2: {9, 10, 11, 12},
		3: {6, 8},
		4: {7},
	},
	graph.IntGraph{
		0: {},
		1: {0},
		2: {1},
		3: {1, 2},
		4: {2, 3},
	},
}

func TestSCC(t *testing.T) {
	t.Run("clrs", func(t *testing.T) { testSCC(t, clrsSCC) })
	t.Run("sedgewick", func(t *testing.T) { testSCC(t, sedgewickSCC) })
	t.Run("sedgewick2", func(t *testing.T) { testSCC(t, sedgewick2SCC) })
}

func testSCC(t *testing.T, test sccTest) {
	scc := SCC(test.g, SCCEdges)

	// Check components.
	var components [][]int
	for i := 0; i < scc.NumNodes(); i++ {
		comp := append([]int{}, scc.SubNodes(i)...)
		sort.Ints(comp)
		components = append(components, comp)
	}
	if !reflect.DeepEqual(test.components, components) {
		t.Errorf("want components:\n%v\ngot:\n%v", test.components, components)
	}

	// Check the edges.
	if !graph.Equal(test.edges, scc) {
		t.Errorf("want edges:\n%s\ngot:\n%s\n",
			graphout.Dot{}.Sprint(test.edges),
			graphout.Dot{}.Sprint(scc),
		)
	}
}
