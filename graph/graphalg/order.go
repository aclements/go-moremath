// Copyright 2018 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graphalg

import (
	"github.com/aclements/go-moremath/graph"
)

// PreOrder returns the nodes of g visited in pre-order.
func PreOrder(g graph.Graph, root int) []int {
	visited := NewNodeMarks()
	out := []int{}
	var visit func(n int)
	visit = func(n int) {
		out = append(out, n)
		visited.Mark(n)
		for _, succ := range g.Out(n) {
			if !visited.Test(succ) {
				visit(succ)
			}
		}
	}
	visit(root)

	return out
}

// PostOrder returns the nodes of g visited in post-order.
func PostOrder(g graph.Graph, root int) []int {
	visited := NewNodeMarks()
	out := []int{}
	var visit func(n int)
	visit = func(n int) {
		visited.Mark(n)
		for _, succ := range g.Out(n) {
			if !visited.Test(succ) {
				visit(succ)
			}
		}
		out = append(out, n)
	}
	visit(root)

	return out
}

// Reverse reverses xs in place and returns the slice. This is useful
// in conjunction with PreOrder and PostOrder to compute reverse
// post-order and reverse pre-order.
func Reverse(xs []int) []int {
	for i, j := 0, len(xs)-1; i < j; i, j = i+1, j-1 {
		xs[i], xs[j] = xs[j], xs[i]
	}
	return xs
}
