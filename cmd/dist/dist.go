// dist reads newline-separated numbers from stdin and describes their
// distribution.
package main

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"

	"github.com/aclements/go-moremath/stats"
)

func main() {
	s := readInput(os.Stdin)
	s.Sort()

	fmt.Printf("N %d  sum %.6g  mean %.6g", len(s.Xs), s.Sum(), s.Mean())
	gmean := s.GeoMean()
	if !math.IsNaN(gmean) {
		fmt.Printf("  gmean %.6g", gmean)
	}
	fmt.Printf("  std dev %.6g  variance %.6g\n", s.StdDev(), s.Variance())
	fmt.Println()

	// Quartiles and tails.
	labels := map[int]string{0: "min", 50: "median", 100: "max"}
	for _, p := range []int{0, 1, 5, 25, 50, 75, 95, 99, 100} {
		label, ok := labels[p]
		if !ok {
			label = fmt.Sprintf("%d%%ile", p)
		}
		fmt.Printf("%8s %.6g\n", label, s.Percentile(float64(p)/100))
	}
	fmt.Println()

	// Kernel density estimate.
	kde := &stats.KDE{Sample: s}
	FprintPDF(os.Stdout, kde)
}

func readInput(r io.Reader) (sample stats.Sample) {
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		l := scanner.Text()
		value, err := strconv.ParseFloat(l, 64)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}

		sample.Xs = append(sample.Xs, value)
	}
	if err := scanner.Err(); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	return
}
