package matrix_test

import (
	"math/rand"
	"testing"
	"time"

	matrix "github.com/andres8bit/numerical_methods/src/matrix"
)

func TestConstructor(t *testing.T) {
	s := rand.NewSource(time.Now().UnixNano())
	r := rand.New(s)

	col := int(r.Int63n(100))
	row := int(r.Int63n(100))

	m := matrix.NewMatrix(row, col)

	if m.ColCount() != col {
		t.Errorf("col count does not match")
	}
	if m.RowCount() != row {
		t.Errorf("row count does not match")
	}
}

func TestGetSet(t *testing.T) {
	s := rand.NewSource(time.Now().UnixNano())
	r := rand.New(s)
	col := int(r.Int63n(100))
	row := int(r.Int63n(100))
	x := int(r.Int63n(int64(row)))
	y := int(r.Int63n(int64(col)))
	val := 1 + 2i

	m := matrix.NewMatrix(row, col)
	m.Set(x, y, val)
	test, _ := m.Get(x, y)

	if test != val {
		t.Errorf("val was not set")

	}

}

func TestGetSetRow(t *testing.T) {
	m := matrix.NewMatrix(10, 10)
	list := [...]complex128{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

	s := rand.NewSource(time.Now().UnixNano())
	r := rand.New(s)
	row := int(r.Int63n(9))

	m.SetRow(list[:], row)
	test, _ := m.GetRow(row)

	for i := 0; i < 10; i++ {
		if test[i] != list[i] {
			t.Errorf("row was not set")

		}
	}

}

func TestGetSetCol(t *testing.T) {
	m := matrix.NewMatrix(10, 10)
	list := [...]complex128{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

	s := rand.NewSource(time.Now().UnixNano())
	r := rand.New(s)
	col := int(r.Int63n(9))

	m.SetCol(list[:], col)
	test, _ := m.GetCol(col)

	for i := 0; i < 10; i++ {
		if test[i] != list[i] {
			t.Errorf("col was not set")

		}
	}
}

func TestCopy(t *testing.T) {
	m := matrix.NewMatrix(10, 10)
	list := [...]complex128{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

	s := rand.NewSource(time.Now().UnixNano())
	r := rand.New(s)
	col := int(r.Int63n(9))
	row := int(r.Int63n(9))

	m.SetCol(list[:], col)
	m.SetRow(list[:], row)
	test := matrix.NewMatrix(9, 9)
	test.Copy(m)

	for i := 0; i < 10; i++ {
		for j := 0; j < 10; j++ {
			val, _ := test.Get(i, j)
			truth, _ := m.Get(i, j)
			if val != truth {
				t.Errorf("copy failed")
			}
		}
	}

}

func TestNaiveElimination(t *testing.T) {
	r1 := [...]complex128{6, -2, 2, 4}
	r2 := [...]complex128{12, -8, 6, 10}
	r3 := [...]complex128{3, -13, 9, 3}
	r4 := [...]complex128{-6, 4, 1, -18}
	b := [...]complex128{16, 26, -19, -34}
	test := [...]complex128{3, 1, -2, 1}

	a := matrix.NewMatrix(4, 4)
	xmat := matrix.NewMatrix(0, 4)
	a.SetRow(r1[:], 0)
	a.SetRow(r2[:], 1)
	a.SetRow(r3[:], 2)
	a.SetRow(r4[:], 3)
	bmat := matrix.NewMatrix(0, 4)
	bmat.SetRow(b[:], 0)
	matrix.Naive(a, xmat, bmat)
	for i := 0; i < 4; i++ {
		truth, _ := xmat.Get(0, i)
		if test[i] != truth {
			t.Errorf("incorrect val")

		}
	}
}

// func TestScaledPartialPivoting(t *testing.T) {
// 	r1 := [...]complex128{3, -13, 9, 3}
// 	r2 := [...]complex128{-6, 4, 1, -18}
// 	r3 := [...]complex128{6, -2, 2, 4}
// 	r4 := [...]complex128{12, -8, 6, 10}
// 	b := [...]complex128{-19, -34, 16, 26}
// 	test := [...]complex128{3, 1, -2, 1}
// 	a := matrix.NewMatrix(4, 4)
// 	xmat := matrix.NewMatrix(0, 4)
// 	a.SetRow(r1[:], 0)
// 	a.SetRow(r2[:], 1)
// 	a.SetRow(r3[:], 2)
// 	a.SetRow(r4[:], 3)
// 	bmat := matrix.NewMatrix(0, 4)
// 	bmat.SetRow(b[:], 0)

// 	matrix.ScaledPartialPivoting(a, xmat, bmat)
// 	for i := 0; i < 4; i++ {
// 		truth, _ := xmat.Get(0, i)
// 		if test[i] != rtruth {
// 			t.Errorf("incorrect val")

// 		}
// 	}

//}
