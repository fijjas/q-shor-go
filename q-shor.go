package main

import (
	"flag"
	"fmt"
	"github.com/itsubaki/q"
	"github.com/itsubaki/q/math/number"
	"github.com/itsubaki/q/math/rand"
	"strings"
)

func main() {
	var N, t, a, seed int
	flag.IntVar(&N, "N", 15, "positive integer")
	flag.IntVar(&t, "t", 3, "precision bits")
	flag.IntVar(&a, "a", -1, "coprime number of N")
	flag.IntVar(&seed, "seed", -1, "PRNG seed for measurements")
	flag.Parse()

	if N < 2 {
		fmt.Printf("N=%d. N must be greater than 1.\n", N)
		return
	}

	if number.IsPrime(N) {
		fmt.Printf("N=%d is prime.\n", N)
		return
	}

	if number.IsEven(N) {
		fmt.Printf("N=%d is even. p=%d, q=%d.\n", N, 2, N/2)
		return
	}

	if a, b, ok := number.BaseExp(N); ok {
		fmt.Printf("N=%d. N is exponentiation. %d^%d.\n", N, a, b)
		return
	}

	if a < 0 {
		a = rand.Coprime(N)
	}

	if a < 2 || a > N-1 {
		fmt.Printf("N=%d, a=%d. a must be 1 < a < N.\n", N, a)
		return
	}

	if number.GCD(N, a) != 1 {
		fmt.Printf("N=%d, a=%d. a is not coprime. a is non-trivial factor.\n", N, a)
		return
	}

	fmt.Printf("N=%d, a=%d, t=%d, seed=%d.\n\n", N, a, t, seed)

	qs := q.New()
	if seed > 0 {
		//qs.Rand = rand.Math(seed)
		qs.Rand = rand.Math
	}

	r0 := qs.ZeroWith(t)
	r1 := qs.ZeroLog2(N)

	qs.X(r1[len(r1)-1])
	dbg("initial state", qs, r0, r1)

	qs.H(r0...)
	dbg("create superposition", qs, r0, r1)

	// qs.CModExp2(a, N, r0, r1)
	// dbg("apply controlled-U", qs, r0, r1)
	for i := 0; i < len(r0); i++ {
		qs.ControlledModExp2(a, i, N, r0[i], r1)
		dbg(fmt.Sprintf("apply controlled-U[%d]", i), qs, r0, r1)
	}

	qs.InvQFT(r0...)
	dbg("apply inverse QFT", qs, r0, r1)

	qs.Measure(r1...)
	dbg("measure reg1", qs, r0, r1)

	for _, state := range qs.State(r0) {
		i, m := state.Value()
		s, r, d, ok := number.FindOrder(a, N, fmt.Sprintf("0.%s", m))
		if !ok || number.IsOdd(r) {
			fmt.Printf("  i=%3d: N=%d, a=%d. s/r=%2d/%2d ([0.%v]~%.4f).\n", i, N, a, s, r, m, d)
			continue
		}

		p0 := number.GCD(number.Pow(a, r/2)-1, N)
		p1 := number.GCD(number.Pow(a, r/2)+1, N)
		if number.IsTrivial(N, p0, p1) {
			fmt.Printf("  i=%3d: N=%d, a=%d. s/r=%2d/%2d ([0.%v]~%.4f). p=%v, q=%v.\n", i, N, a, s, r, m, d, p0, p1)
			continue
		}

		fmt.Printf("* i=%3d: N=%d, a=%d. s/r=%2d/%2d ([0.%v]~%.4f). p=%v, q=%v.\n", i, N, a, s, r, m, d, p0, p1)
	}
}

func dbg(desc string, qsim *q.Q, reg ...any) {
	fmt.Println(desc)

	max := number.Max(qsim.Probability())
	for _, s := range qsim.State(reg...) {
		p := strings.Repeat("*", int(s.Probability/max*32))
		fmt.Printf("%s: %s\n", s, p)
	}

	fmt.Println()
}
