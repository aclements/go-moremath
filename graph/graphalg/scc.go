// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graphalg

import (
	"sort"

	"github.com/aclements/go-moremath/graph"
)

// SCCFlags is a set of optional analyses to perform when constructing
// strongly-connected components.
type SCCFlags int

const (
	// SCCSubnodeComponent instructs SCC to record a mapping from
	// subnode to component ID containing that subnode.
	SCCSubnodeComponent SCCFlags = 1 << iota

	// SCCEdges instructs SCC to record edges between components.
	// Otherwise, the resulting SCC graph will have a node for
	// each strongly-connected component, but no edges.
	SCCEdges
)

// SCC computes the strongly-connected components of graph g.
//
// This implements Tarjan's strongly connected components algorithm
// [1]. It runs in O(V + E) time and O(V) space.
//
// [1] Tarjan, R. E. (1972), "Depth-first search and linear graph
// algorithms", SIAM Journal on Computing, 1 (2): 146–160.
func SCC(g graph.Graph, flags SCCFlags) *SCCGraph {
	var sccs SCCGraph

	if flags&SCCEdges != 0 {
		// Edge construction requires sub-graph ID ->
		// component ID mapping.
		flags |= SCCSubnodeComponent
	}

	// This is based on the presentation of Tarjan's algorithm in
	// Sedgewick, Algorithms in C, Part 5, Third Edition, p. 202.
	// This is a fair bit simpler than Tarjan's original
	// presentation. We further simplify it by combining "pre" and
	// "low" into just "low", since pre is only ever used as a
	// visited mark. We instead start node indexing at 1 and use
	// low[x] == 0 to indicate node x has not been visited.
	numNodes := g.NumNodes()
	// For low, 0 means "not visited", ^uint(0) means "processed".
	low := make([]uint, numNodes)
	stack := []int{}
	index := uint(1)

	// We construct out-edges of each component by maintaining a
	// parallel stack of seen out-edges. As we pop a component off
	// the primary stack, we match it with this sack by recording
	// the length of the stack when an edge was pushed.
	type outEdge struct {
		cid      int
		stackLen int
	}
	var out []outEdge

	var connect func(nid int)
	connect = func(nid int) {
		// Set the depth of n to the next unused index.
		low[nid] = index
		min := index
		index++
		stackPos := len(stack)
		stack = append(stack, nid)

		// Process successors of n.
		for _, oid := range g.Out(nid) {
			if low[oid] == 0 {
				// Successor has not yet been visited.
				connect(oid)
			}
			if low[oid] < min {
				min = low[oid]
			}

			if flags&SCCEdges != 0 && low[oid] == ^uint(0) {
				// Successor is in another component.
				// Record the out-edge from this
				// component.
				cid := sccs.subnodeComponent[oid]
				out = append(out, outEdge{cid, stackPos})
			}
		}

		if min < low[nid] {
			// Node n is not the root of an SCC.
			low[nid] = min
			return
		}

		// Node n is a root of an SCC. Pop the stack to
		// construct the component.
		cid := len(sccs.subnodeIndexes)
		var i int
		for i = len(stack) - 1; i >= 0; i-- {
			oid := stack[i]
			// Set low such that it can never be less than
			// min. This also indicates we've connected
			// oid to a component.
			low[oid] = ^uint(0)
			if flags&SCCSubnodeComponent != 0 {
				sccs.subnodeComponent[oid] = cid
			}
			if oid == nid {
				break
			}
		}
		sccs.subnodeIndexes = append(sccs.subnodeIndexes, len(sccs.subnodes))
		sccs.subnodes = append(sccs.subnodes, stack[i:]...)
		stack = stack[:i]

		// Collect out-edges of this SCC.
		if flags&SCCEdges != 0 {
			outStart := len(sccs.out)
			sccs.outIndexes = append(sccs.outIndexes, outStart)
			// Pop the out-edge stack until it
			// aligns with the node stack.
			for i = len(out) - 1; i >= 0; i-- {
				if out[i].stackLen < len(stack) {
					break
				}
				sccs.out = append(sccs.out, out[i].cid)
			}
			i++
			out = out[:i]
			// Dedup component IDs.
			sort.Ints(sccs.out[outStart:])
			i = outStart
			for j := outStart; j < len(sccs.out); j++ {
				if i == outStart || sccs.out[i-1] != sccs.out[j] {
					sccs.out[i] = sccs.out[j]
					i++
				}
			}
			sccs.out = sccs.out[:i]
		}
	}

	sccs.subnodes = make([]int, 0, numNodes)
	if flags&SCCSubnodeComponent != 0 {
		sccs.subnodeComponent = make([]int, numNodes)
	}

	for nid := range low {
		// If node n is not yet visited, then connect it.
		if low[nid] == 0 {
			connect(nid)
		}
	}
	sccs.subnodeIndexes = append(sccs.subnodeIndexes, len(sccs.subnodes))
	if flags&SCCEdges != 0 {
		sccs.outIndexes = append(sccs.outIndexes, len(sccs.out))
	}

	return &sccs
}

// SCCGraph is a set of strongly-connected components of another
// graph.
//
// Each strongly-connected component is a node in this graph.
// The components are numbered in reverse topological sort order.
//
// If the graph was constructed with flag SCCEdges, then it also has
// edges between the components that follow the edges in the
// underlying graph.
type SCCGraph struct {
	subnodes       []int // Concatenated list of sub-graph nodes in each component
	subnodeIndexes []int // Component ID -> subnodes base index

	subnodeComponent []int // Sub-node ID -> component ID

	out        []int // Concatenated list of out-edges of each component
	outIndexes []int // Component ID -> out base index
}

// Subnodes returns the IDs of the nodes in the underlying graph that
// comprise component cid.
func (g *SCCGraph) Subnodes(cid int) []int {
	return g.subnodes[g.subnodeIndexes[cid]:g.subnodeIndexes[cid+1]]
}

// SubnodeComponent returns the component ID (a node ID in g) of
// sub-graph node subID (a node ID in the underlying graph).
//
// Graph g must have been constructed with flag SCCSubnodeComponent.
func (g *SCCGraph) SubnodeComponent(subID int) (componentID int) {
	if g.subnodeComponent == nil {
		panic("SCCGraph constructed without SCCSubnodeComponent flag")
	}
	return g.subnodeComponent[subID]
}

// NumNodes returns the number of strongly-connected components in g.
func (g *SCCGraph) NumNodes() int {
	return len(g.subnodeIndexes) - 1
}

// Out returns the IDs of the components for which there are any edges
// in the underlying graph from component cid.
//
// Graph g must have been constructed with flag SCCEdges. Otherwise
// this returns nil.
func (g *SCCGraph) Out(cid int) []int {
	if g.out == nil {
		return nil
	}
	return g.out[g.outIndexes[cid]:g.outIndexes[cid+1]]
}
