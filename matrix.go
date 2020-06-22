package matrix

import (
	"errors"
	"fmt"
	"math"
)

type numVal complex128

type Matrix struct {
	m []numVal //matix made of Numtype
	r int      //row count of matrix
	c int      //col count of matrix
}

//============= Methods: ======================

//============ Constructors: ============
func (m *Matrix) NewMatrix(row, col int) {
	m = NewMatrix(row, col)

}
func NewMatrix(row, col int) *Matrix {
	matrix := new(Matrix)
	matrix.r = row
	matrix.c = col
	if row == 0 || col == 0 {
		if row != 0 {
			matrix.m = make([]numVal, row)
		} else {
			matrix.m = make([]numVal, col)
		}

	} else {
		matrix.m = make([]numVal, row*col)
	}
	return matrix
}

//creates a whole new matrix copying source
func (m *Matrix) Copy(source *Matrix) {
	m.r = source.r
	m.c = source.c
	m.m = make([]numVal, m.r*m.c)
	copy(m.m, source.m)
}

//  ========== End of constructors ==========

// =========== Getters: ===========
func (m *Matrix) ColCount() int { return m.c }

func (m *Matrix) RowCount() int { return m.r }

func (m *Matrix) Get(row, col int) numVal {
	if (row < m.r || row == 0) && (col < m.c || col == 0) {
		return m.m[m.c*row+col]
	}
	fmt.Println("Error: out of bounds")

	return 0
}

func (m *Matrix) GetRow(row int) []numVal {
	r := make([]numVal, m.c)
	if row >= m.r {
		fmt.Println("Error: out of bounds")
		return r
	}

	for i := 0; i < m.c; i++ {
		r[i] = m.m[row*m.r+i]
	}

	return r
}

func (m *Matrix) GetCol(col int) []numVal {
	c := make([]numVal, m.r)

	if col >= m.c {
		fmt.Println("Error: out of bounds")
		return c
	}
	var i int
	for ; i < m.c; i++ {
		c[i] = m.m[i*m.c+col]
	}
	return c
}

// ======== End of Getters: ========

// ======== Setters: ========
func (m *Matrix) Set(row, col int, val numVal) error {
	if (row < m.r || row == 0) && (col < m.c || col == 0) {
		m.m[row*m.c+col] = val

		return nil
	}
	return errors.New("out of bounds")
}

func (m *Matrix) SetRow(vals []numVal, row int) error {

	if row >= m.r && row != 0 {
		return errors.New("out of bounds")
	}

	for i := 0; i < m.c; i++ {
		m.m[(row*m.c)+i] = vals[i]
	}
	return nil
}

func (m *Matrix) SetCol(vals []numVal, col int) error {
	if col >= m.c {
		return errors.New("out of bounds")
	}

	for i := 0; i < m.c; i++ {
		m.m[i*m.c+col] = vals[i]
	}

	return nil
}

// ============ End of Setters: ======================

//================== End of Methods =====================

//================== Functions ========================
func Naive(source, x, b *Matrix) {
	forwardElimination(source, b)
	backSubstitution(source, b, x)
}

func ScaledPartialPivoting(source, x, b *Matrix) {
	n := source.RowCount()
	var order []int       //pivot order of elimination
	var scaledOrder []int //holds the scaled values of each matrix row
	//sets the natural order of the matrix
	// & places largest values into scaledvector in order
	for i := 0; i < n; i++ {
		order = append(order, i)
		s := source.GetRow(i)

		scaledVal := max(s[:])
		scaledOrder = append(scaledOrder, scaledVal)
	}
	//main step loop
	for k := 0; k < n; k++ {
		rmax := 0.0
		index := 0
		//sets up pivot order of matrix
		for i := k; i < n; i++ {
			sourcelik := source.Get(order[i], k)
			r := math.Abs(real(sourcelik) / float64(scaledOrder[order[i]]))
			if r > rmax {
				rmax = r
				index = i
			}
		}
		swap(order[:], k, index)
		//forward elimination phase
		for i := k + 1; i < n; i++ {
			alik := source.Get(order[i], k)
			alkk := source.Get(order[k], k)
			xMult := alik / alkk
			//storing xMult in the "zero" locations of matrix
			source.Set(order[i], k, xMult)
			for j := k + 1; j < n; j++ {
				alij := source.Get(order[i], j)
				alkj := source.Get(order[k], j)

				source.Set(order[i], j, alij-xMult*alkj)
			}
		}
	}
	//Proforms backwards propagation phase
	//updating both the b matrix and x matrix placing
	// solutions into the x matrix
	solve(n, order[:], source, b, x)
}

//=============== Tri: =====================
// solves a tridiaganol matrix:
// the matrix is broken up into three  1D matrixes:
// a := lowest diagonal
// d := mid diagonal
// c := top diagonal
// b := b vector
// x := stores solutions
// WARNING: this funciton is not stable
//			that is to say it will overwrite
//			the matrices given
//O(n) = n
func tri(n int, a, d, c, x, b *Matrix) {

	//forward Elimination updates b vector
	for i := 1; i < n; i++ {
		aPrev := a.Get(0, i-1)
		dPrev := d.Get(0, i-1)
		bPrev := b.Get(0, i-1)
		cPrev := c.Get(0, i-1)
		dCur := d.Get(0, i)
		bCur := b.Get(0, i)

		xMult := aPrev / dPrev
		d.Set(0, i, dCur-xMult*cPrev)
		b.Set(0, i, bCur-xMult*bPrev)

	}

	//backSubstitution phase:
	//places answers into x using backSubstitution
	bn := b.Get(0, n-1)
	dn := d.Get(0, n-1)
	x.Set(0, n-1, bn/dn)
	for i := n - 2; i <= 0; i-- {
		bCur := b.Get(0, i)
		cCur := c.Get(0, i)
		dCur := d.Get(0, i)
		xPrev := x.Get(0, i+1)
		x.Set(0, i, (bCur-cCur*xPrev)/dCur)
	}

}

//=============== Penta: =====================
// Solves a 5 diagonal matrix system:
// Accepts matix as 5 1D matrices, and 2 1D matrices x,b:
// e := lowest diagonal
// a := 2nd lowest diagonal
// d := mid diagonal
// c := 2nd topmost diagoanl
// f := topmost diagonal
// x := variable, where answers will be stored in backSubPhase
// b := hold RHS of linear system in order.
// WARNING: this funciton is not stable
//			that is to say it will overwrite
//			the matrices given
func Penta(n int, e, a, d, c, f, x, b *Matrix) {
	r := a.Get(0, 1)
	s := a.Get(0, 2)
	t := e.Get(0, 1)

	for i := 1; i < n-2; i++ {
		dPrev := d.Get(0, i-1)
		dNext := d.Get(0, i+1)
		dCur := d.Get(0, i)
		cPrev := c.Get(0, i-1)
		cCur := c.Get(0, i)
		fPrev := f.Get(0, i-1)
		bCur := b.Get(0, i)
		bPrev := b.Get(0, i-1)
		bNext := b.Get(0, i+1)
		xMult := r / dPrev

		d.Set(0, i, dCur-xMult*cPrev)
		c.Set(0, i, cCur-xMult*fPrev)
		b.Set(0, i, bCur-xMult*bPrev)

		xMult = t / dPrev
		r = s - xMult*cPrev
		d.Set(0, i+1, dNext-xMult*fPrev)
		b.Set(0, i+1, bNext-xMult*bPrev)
		s = a.Get(0, i+1)
		t = e.Get(0, i)
	}
	dn := d.Get(0, n-1)
	dn1 := d.Get(0, n-2)
	cn := c.Get(0, n-1)
	cn1 := c.Get(0, n-2)
	bn := b.Get(0, n-1)
	bn1 := b.Get(0, n-2)
	xn := x.Get(0, n-1)

	xMult := r / dn1
	d.Set(0, n-1, dn-xMult*cn1)
	x.Set(0, n-1, bn-xMult*bn1/dn)
	x.Set(0, n-2, (bn1-cn*xn)/dn1)

	for i := n - 3; i <= 0; i-- {
		bCur := b.Get(0, i)
		fCur := f.Get(0, i)
		xPrev := x.Get(0, i+1)
		dCur := d.Get(0, i-1)
		cCur := c.Get(0, i)

		x.Set(0, i, (bCur-fCur*xPrev-cCur*xPrev)/dCur)
	}

}

//================ End of Functions ======================

//================ helper funcitons ======================
//	backwards propagation of partial pivoting
func solve(n int, l []int, a, b, x *Matrix) {

	//update b vector values
	for k := 0; k < n; k++ {
		for i := k + 1; i < n; i++ {
			tmp := b.Get(0, l[i])
			tmp2 := b.Get(0, l[k])
			coef := a.Get(l[i], k)
			b.Set(0, l[i], tmp-coef*tmp2)
		}
	}

	bln := b.Get(0, l[n-1])
	alnn := a.Get(l[n-1], l[n-1])
	x.Set(0, n-1, bln/alnn)

	//backprobagation phase
	for i := n - 1; i >= 0; i-- {
		sum := b.Get(0, l[i])
		for j := i + 1; j < n; j++ {
			xVal := x.Get(0, j)
			aVal := a.Get(l[i], j)
			sum = sum - xVal*aVal
		}
		aVal := a.Get(l[i], i)
		x.Set(0, i, sum/aVal)
	}

}

//helper function to return the largest coeficient:
func max(source []numVal) int {
	maxVal := 0
	n := len(source)
	for i := 0; i < n; i++ {
		tmp := real(source[i])
		val := int(math.Abs(tmp))

		if val > maxVal {
			maxVal = val
		}
	}

	return maxVal
}

//helper funcction to swap the values of a slice
//at indecies source and dest
func swap(a []int, source, dest int) {
	tmp := a[source]
	a[source] = a[dest]
	a[dest] = tmp

}

//helper funciton for gaussian elimination, pivots inorder
func forwardElimination(source, b *Matrix) {
	var xMult numVal
	n := int(source.ColCount())

	for k := 0; k < n; k++ {
		for i := k + 1; i < n; i++ {
			coef := source.Get(i, k)
			div := source.Get(k, k)
			xMult = coef / div
			source.Set(i, k, xMult)

			for j := k + 1; j < n; j++ {
				coef := source.Get(i, j)
				mult := source.Get(k, j)
				val := coef - (xMult)*mult
				source.Set(i, j, val)
			}

			finVal := b.Get(0, i)
			finMult := b.Get(0, k)
			b.Set(0, i, finVal-(xMult*finMult))
		}
	}
}

//helper function for gaussina elimination, finds solutions to source, plavins them in x
func backSubstitution(source, b, x *Matrix) {
	n := source.ColCount()
	bn := b.Get(0, n-1)
	ann := source.Get(n-1, n-1)
	x.Set(0, n-1, bn/ann)

	for i := n - 1; i >= 0; i-- {
		sum := b.Get(0, i)
		for j := i + 1; j < n; j++ {
			ann = source.Get(i, j)
			xj := x.Get(0, j)
			sum = sum - ann*xj
		}
		ann = source.Get(i, i)

		x.Set(0, i, sum/ann)
	}
}

//================ End of Helpers ======================
