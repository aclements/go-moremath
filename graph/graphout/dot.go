// Copyright 2018 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package graphout

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/aclements/go-moremath/graph"
)

// Dot contains options for generating a Graphviz Dot graph from a
// Graph.
type Dot struct {
	// Name is the name given to the graph. Usually this can be
	// left blank.
	Name string

	// Label returns the string to use as a label for the given
	// node. If nil, nodes are labeled with their node numbers.
	Label func(node int) string

	// NodeAttrs, if non-nil, returns a set of attributes for a
	// node. If this includes a "label" attribute, it overrides
	// the label returned by Label.
	NodeAttrs func(node int) []DotAttr

	// EdgeAttrs, if non-nil, returns a set of attributes for an
	// edge.
	EdgeAttrs func(node, edge int) []DotAttr
}

// DotAttr is an attribute for a Dot node or edge.
type DotAttr struct {
	Name string
	// Val is the value of this attribute. It may be a string
	// (which will be escaped), bool, int, uint, float64 or
	// DotLiteral.
	Val interface{}
}

// DotLiteral is a string literal that should be passed to dot
// unescaped.
type DotLiteral string

func defaultLabel(node int) string {
	return fmt.Sprintf("%d", node)
}

// Print writes the Dot form of g to os.Stdout.
func (d Dot) Print(g graph.Graph) error {
	return d.Fprint(g, os.Stdout)
}

// Sprint returns the Dot form of g as a string.
func (d Dot) Sprint(g graph.Graph) string {
	var buf strings.Builder
	d.Fprint(g, &buf)
	return buf.String()
}

// Fprint writes the Dot form of g to w.
func (d Dot) Fprint(g graph.Graph, w io.Writer) error {
	label := d.Label
	if label == nil {
		label = defaultLabel
	}

	_, err := fmt.Fprintf(w, "digraph %s {\n", DotString(d.Name))
	if err != nil {
		return err
	}

	for i := 0; i < g.NumNodes(); i++ {
		// Define node.
		var attrList []DotAttr
		var haveLabel bool
		if d.NodeAttrs != nil {
			attrList = d.NodeAttrs(i)
			for _, attr := range attrList {
				if attr.Name == "label" {
					haveLabel = true
					break
				}
			}
		}
		if !haveLabel {
			attrList = attrList[:len(attrList):len(attrList)]
			attrList = append(attrList, DotAttr{"label", label(i)})
		}
		_, err = fmt.Fprintf(w, "n%d%s;\n", i, formatAttrs(attrList))
		if err != nil {
			return err
		}

		// Connect node.
		for j, out := range g.Out(i) {
			var attrs string
			if d.EdgeAttrs != nil {
				attrs = formatAttrs(d.EdgeAttrs(i, j))
			}
			_, err = fmt.Fprintf(w, "n%d -> n%d%s;\n", i, out, attrs)
			if err != nil {
				return err
			}
		}
	}

	_, err = fmt.Fprintf(w, "}\n")
	return err
}

// DotString returns s as a quoted dot string.
//
// Users of the Dot type don't need to call this, since it will
// automatically quote strings. However, this is useful for building
// custom dot output.
func DotString(s string) string {
	buf := []byte{'"'}
	for i := 0; i < len(s); i++ {
		switch s[i] {
		case '\n':
			buf = append(buf, '\\', 'n')
		case '\\', '"', '{', '}', '<', '>', '|':
			// TODO: Option to allow formatting
			// characters? Maybe private use code points
			// to encode formatting characters? Or
			// something more usefully structured?
			buf = append(buf, '\\', s[i])
		default:
			buf = append(buf, s[i])
		}
	}
	buf = append(buf, '"')
	return string(buf)
}

// formatAttrs formats attrs as a dot attribute set, including the
// surrounding brackets. If attrs is empty, it returns an empty
// string.
func formatAttrs(attrs []DotAttr) string {
	if len(attrs) == 0 {
		return ""
	}
	var buf strings.Builder
	buf.WriteString(" [")
	for i, attr := range attrs {
		if i > 0 {
			buf.WriteString(",")
		}
		buf.WriteString(attr.Name)
		buf.WriteString("=")
		switch val := attr.Val.(type) {
		case string:
			buf.WriteString(DotString(val))
		case int, uint, float64:
			fmt.Fprintf(&buf, "%v", val)
		case DotLiteral:
			buf.WriteString(string(val))
		default:
			panic(fmt.Sprintf("dot attribute %s had unknown type %T", attr.Name, attr.Val))
		}
	}
	buf.WriteString("]")
	return buf.String()
}
