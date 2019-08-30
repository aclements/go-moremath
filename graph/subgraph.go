// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graph

// A Subgraph is a Graph that consists of a subset of the nodes and
// vertices from another, underlying Graph.
type Subgraph interface {
	Graph

	// Underlying returns the underlying graph that this is a
	// subgraph of.
	Underlying() Graph

	// NodeMap transduces a node property map on the underlying
	// graph into a node property map on this graph.
	NodeMap(underlyingMap func(node int) interface{}) func(node int) interface{}

	// EdgeMap transduces an edge property map on the underlying
	// graph into an edge property map on this graph.
	EdgeMap(underlyingMap func(node, edge int) interface{}) func(node, edge int) interface{}
}

// SubgraphKeep returns a subgraph of g that keeps the given nodes and
// edges. Subgraph node i corresponds to nodes[i] in g.
func SubgraphKeep(g Graph, nodes []int, edges []Edge) Subgraph {
	// Create old-to-new node mapping.
	gNodes := g.NumNodes()
	oldToNew := make(map[int]int, len(nodes))
	for newNode, oldNode := range nodes {
		if oldNode < 0 || oldNode >= gNodes {
			panic("node not in underlying graph")
		}
		if _, ok := oldToNew[oldNode]; ok {
			panic("duplicate node")
		}
		oldToNew[oldNode] = newNode
	}

	// Construct new nodes.
	newNodes := make([]listSubgraphNode, len(nodes))
	for i, oldNode := range nodes {
		newNodes[i].oldNode = oldNode
	}

	// Map old edge indexes to new node IDs.
	for _, oldEdge := range edges {
		newNode := &newNodes[oldToNew[oldEdge.Node]]
		oldTo := g.Out(oldEdge.Node)[oldEdge.Edge]
		newTo := oldToNew[oldTo]

		newNode.out = append(newNode.out, newTo)
		newNode.oldEdges = append(newNode.oldEdges, oldEdge.Edge)
	}

	return &listSubgraph{g, newNodes}
}

// SubgraphRemove returns a subgraph of g that removes the given nodes
// and edges from g, as well as all edges incident to those nodes.
func SubgraphRemove(g Graph, nodes []int, edges []Edge) Subgraph {
	// Collect the set of nodes and edges to remove.
	rmNodes := make(map[int]struct{}, len(nodes))
	for _, node := range nodes {
		rmNodes[node] = struct{}{}
	}
	rmEdges := make(map[Edge]struct{}, len(edges))
	for _, edge := range edges {
		rmEdges[edge] = struct{}{}
	}

	// Create new-to-old and old-to-new node mappings.
	newNodes := make([]listSubgraphNode, 0, g.NumNodes()-len(rmNodes))
	oldToNew := make(map[int]int, cap(newNodes))
	for oldNode := 0; oldNode < g.NumNodes(); oldNode++ {
		if _, ok := rmNodes[oldNode]; ok {
			continue
		}
		newNode := len(newNodes)
		newNodes = append(newNodes, listSubgraphNode{oldNode: oldNode})
		oldToNew[oldNode] = newNode
	}

	// Create edge mappings.
	for i := range newNodes {
		newNode := &newNodes[i]
		oldNode := newNode.oldNode
		oldOut := g.Out(oldNode)
		for j, oldNode2 := range oldOut {
			if _, ok := rmNodes[oldNode2]; ok {
				// Target node removed.
				continue
			}
			if _, ok := rmEdges[Edge{oldNode, j}]; ok {
				// Edge removed.
				continue
			}
			newNode.out = append(newNode.out, oldToNew[oldNode2])
			newNode.oldEdges = append(newNode.oldEdges, j)
		}
	}

	return &listSubgraph{g, newNodes}
}

type listSubgraph struct {
	underlying Graph
	nodes      []listSubgraphNode
}

type listSubgraphNode struct {
	out      []int // Adjacency list
	oldNode  int   // Node ID in underlying graph
	oldEdges []int // New edge index -> old edge index
}

func (s *listSubgraph) NumNodes() int {
	return len(s.nodes)
}

func (s *listSubgraph) Out(node int) []int {
	return s.nodes[node].out
}

func (s *listSubgraph) Underlying() Graph {
	return s.underlying
}

func (s *listSubgraph) NodeMap(underlyingMap func(node int) interface{}) func(node int) interface{} {
	return func(node int) interface{} {
		return underlyingMap(s.nodes[node].oldNode)
	}
}

func (s *listSubgraph) EdgeMap(underlyingMap func(node, edge int) interface{}) func(node, edge int) interface{} {
	return func(node, edge int) interface{} {
		newNode := &s.nodes[node]
		return underlyingMap(newNode.oldNode, newNode.oldEdges[edge])
	}
}
