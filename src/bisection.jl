"""
- x_lo: lower bound for solution, should be finite
- x_hi: upper bound for solution, should be finite
- pred(x): predicate function
- mid(x,y): a midpoint function which returns a point between x <= y
- dist(x,y): a a distance functions between x <= y
- tol: Stop when dist(x_lo,x_hi) < tol
- Watch out for numerical overflow. eventually mid(x,y) will no longer be
  strictly between x and y

# Example: find the smallest x in [0,1] where x^2 >= 0.5
eucdist(x,y) = y - x
midpoint(x,y) = (x + y) / 2
pred(x) = (x^2 >= 0.5)
x = DirectSampling.bisection(0, 1, pred, midpoint, eucdist, tol = 1e-20)
"""

function bisection(x_lo, x_hi, pred::Function, mid::Function, dist::Function; tol = 1e-10)
	x = mid(x_lo, x_hi)

	while (dist(x_lo, x_hi) > tol && x > x_lo && x < x_hi)
		ind = pred(x)
		x_lo = ind*x_lo + (1-ind)*x
		x_hi = ind*x + (1-ind)*x_hi
		x = mid(x_lo, x_hi)
	end

	return x
end



