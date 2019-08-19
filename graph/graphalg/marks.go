// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graphalg

// NodeMarks is a structure for marking nodes in a graph.
type NodeMarks struct {
	marks []uint32
}

// Test returns whether node i is marked.
func (m NodeMarks) Test(i int) bool {
	if i < 0 || i/32 >= len(m.marks) {
		return false
	}
	return m.marks[i/32]&(1<<uint(i%32)) != 0
}

// Mark marks node i.
func (m *NodeMarks) Mark(i int) {
	if i/32 >= len(m.marks) {
		m.grow(i)
	}
	m.marks[i/32] |= 1 << uint(i%32)
}

// Unmark clears the mark on node i.
func (m *NodeMarks) Unmark(i int) {
	if i/32 >= len(m.marks) {
		return
	}
	m.marks[i/32] &^= 1 << uint(i%32)
}

func (m *NodeMarks) grow(i int) {
	n := i/32 + 1
	// Round n up to a power of two.
	k := 1
	for k > n {
		k <<= 1
	}
	marks := make([]uint32, k)
	copy(marks, m.marks)
	m.marks = marks
}

// NewNodeMarks returns a node mark set with no marks set.
func NewNodeMarks() *NodeMarks {
	// This is small enough to get inlined, allowing the initial
	// marks slice to get stack-allocated.
	return &NodeMarks{make([]uint32, 1024/32)}
}
