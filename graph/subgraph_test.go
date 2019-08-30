// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graph

import (
	"bytes"
	"fmt"
	"testing"
)

var graph1 = IntGraph{
	0: {0, 1},
	1: {0, 2},
	2: {3, 4},
	3: {4},
	4: {},
}

func TestSubgraphKeep(t *testing.T) {
	var want = IntGraph{
		0: {},     // Was node 4
		1: {2},    // Was node 1
		2: {1, 2}, // Was node 0
	}
	g2 := SubgraphKeep(graph1, []int{4, 1, 0}, []Edge{{0, 1}, {0, 0}, {1, 0}})
	if !Equal(want, g2) {
		t.Fatalf("want:\n%sgot:\n%s", pgraph(want), pgraph(g2))
	}

	nodeMap := g2.NodeMap(func(node int) interface{} { return node })
	newToOldNode := []int{4, 1, 0}
	for newNode, oldNode := range newToOldNode {
		if got := nodeMap(newNode).(int); got != oldNode {
			t.Errorf("nodeMap(%d) = %d, want %d", newNode, got, oldNode)
		}
	}

	edgeMap := g2.EdgeMap(func(node, edge int) interface{} { return Edge{node, edge} })
	newToOldEdge := map[Edge]Edge{
		Edge{1, 0}: Edge{1, 0},
		Edge{2, 0}: Edge{0, 1},
		Edge{2, 1}: Edge{0, 0},
	}
	for newEdge, oldEdge := range newToOldEdge {
		if got := edgeMap(newEdge.Node, newEdge.Edge); got != oldEdge {
			t.Errorf("edgeMap(%d, %d) = %v, want %v", newEdge.Node, newEdge.Edge, got, oldEdge)
		}
	}
}

func TestSubgraphRemove(t *testing.T) {
	// Test automatic edge removal.
	var want = IntGraph{
		0: {0}, // Was node 0
		1: {2}, // Was node 3
		2: {},  // Was node 4
	}
	g2 := SubgraphRemove(graph1, []int{1, 2}, nil)
	if !Equal(want, g2) {
		t.Fatalf("want:\n%sgot:\n%s", pgraph(want), pgraph(g2))
	}

	// Test edge removal.
	want = IntGraph{
		0: {1},
		1: {0},
		2: {3, 4},
		3: {4},
		4: {},
	}
	g2 = SubgraphRemove(graph1, nil, []Edge{{0, 0}, {1, 1}})
	if !Equal(want, g2) {
		t.Fatalf("want:\n%sgot:\n%s", pgraph(want), pgraph(g2))
	}
}

func pgraph(g Graph) string {
	var buf bytes.Buffer
	for nid := 0; nid < g.NumNodes(); nid++ {
		fmt.Fprintf(&buf, "%d ->", nid)
		for _, n2 := range g.Out(nid) {
			fmt.Fprintf(&buf, " %d", n2)
		}
		fmt.Fprintf(&buf, "\n")
	}
	return buf.String()
}
