// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graph

import "sort"

// Equal returns true if g1 and g2 have identical nodes and edges,
// including the IDs of all nodes.
func Equal(g1, g2 Graph) bool {
	n := g1.NumNodes()
	if n != g2.NumNodes() {
		return false
	}
	var temp []int
	for i := 0; i < n; i++ {
		e1 := g1.Out(i)
		e2 := g2.Out(i)
		if len(e1) != len(e2) {
			return false
		}
		// Quick check to see if they're identical without
		// sorting.
		eq := true
		for ei, x := range e1 {
			if e2[ei] != x {
				eq = false
				break
			}
		}
		if eq {
			continue
		}
		// Sort the adjacency list and check equality again.
		temp = append(append(temp[:0], e1...), e2...)
		e1, e2 = temp[:len(e1)], temp[len(e1):]
		sort.Ints(e1)
		sort.Ints(e2)
		for ei, x := range e1 {
			if e2[ei] != x {
				return false
			}
		}
	}

	return true
}
