// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graphalg

import (
	"github.com/aclements/go-moremath/graph"
)

// Euler visits a graph using an Euler tour.
//
// For a tree, the Euler tour is well-defined and unique (given an
// ordering of the children of a node). For a general graph, this uses
// the tree formed by the pre-order traversal of the graph.
type Euler struct {
	// Enter is called when a node a first visited. It may be nil.
	Enter func(n int)

	// Exit is called when all of the children of n have been
	// visited. It may be nil.
	//
	// Calls to Enter and Exit are always paired in nested order.
	Exit func(n int)
}

// Visit performs a Euler tour over g starting at root and invokes the
// callbacks on e.
func (e Euler) Visit(g graph.Graph, root int) {
	visited := NewNodeMarks()
	var visit func(n int)
	visit = func(n int) {
		if e.Enter != nil {
			e.Enter(n)
		}
		visited.Mark(n)
		for _, succ := range g.Out(n) {
			if !visited.Test(succ) {
				visit(succ)
			}
		}
		if e.Exit != nil {
			e.Exit(n)
		}
	}
	visit(root)
}
