template < class Matrix, class Vector, class Preconditioner, class Real >
int
CG(const Matrix &A, Vector &x, const Vector &b,
const Preconditioner &M, int &max_iter, Real &tol)
{
	Real resid;
	Vector p, z, q;
	Vector alpha(1), beta(1), rho(1), rho_1(1);

	Real normb = norm(b);
	Vector r = b - A*x;

	if (normb == 0.0)
		normb = 1;

	if ((resid = norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		z = M.solve(r);
		rho(0) = dot(r, z);

		if (i == 1)
			p = z;
		else {
			beta(0) = rho(0) / rho_1(0);
			p = z + beta(0) * p;
		}

		q = A*p;
		alpha(0) = rho(0) / dot(p, q);

		x += alpha(0) * p;
		r -= alpha(0) * q;

		if ((resid = norm(r) / normb) <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		rho_1(0) = rho(0);
	}

	tol = resid;
	return 1;
}
