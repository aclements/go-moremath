// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graphalg

import (
	"reflect"
	"testing"
)

func TestNodeMarksNext(t *testing.T) {
	tests := [][]int{
		{0},
		{1},
		{0, 4},
		{},       // No marks
		{0, 100}, // Big gap
	}

	for _, test := range tests {
		m := NewNodeMarks()
		for _, id := range test {
			m.Mark(id)
		}
		got := []int{}
		for i := m.Next(-1); i >= 0; i = m.Next(i) {
			got = append(got, i)
		}
		if !reflect.DeepEqual(test, got) {
			t.Errorf("want %v, got %v", test, got)
		}
	}
}
