/*571
 https://people.math.sc.edu/Burkardt/cpp_src/r8lib/r8lib.html
 List of Routines:
 GAMMA_VALUES returns some values of the Gamma function.
 GAMMA_LOG_VALUES returns some values of the Log Gamma function.
 I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
 I4_MAX returns the maximum of two I4's.
 I4_MIN returns the minimum of two I4's.
 I4_MODP returns the nonnegative remainder of I4 division.
 I4_POWER returns the value of I^J.
 I4_SIGN returns the sign of an I4.
 I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
 I4_WRAP forces an I4 to lie between given limits by wrapping.
 I4INT_TO_R8INT maps an I4 interval to an R8 interval.
 I4VEC_COPY copies an I4VEC.
 I4VEC_INDICATOR0_NEW sets an I4VEC to the indicator vector (0,1,2,...).
 I4VEC_INDICATOR1_NEW sets an I4VEC to the indicator vector (1,2,3,...).
 I4VEC_PERMUTE permutes an I4VEC in place.
 I4VEC_PRINT prints an I4VEC.
 I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
 I4VEC_ZEROS zeroes an I4VEC.
 I4VEC_ZEROS_NEW creates and zeroes an I4VEC.
 LEGENDRE_ZEROS returns the zeros of the Legendre polynomial of degree N.
 PERM0_CHECK checks a permutation of ( 0, ..., N-1 ).
 PERM0_UNIFORM_NEW selects a random permutation of 0,...,N-1.
 PERM1_CHECK checks a permutation of (1, ..., N ).
 PERM1_UNIFORM_NEW selects a random permutation of 1,...,N.
 R8_ABS returns the absolute value of an R8.
 R8_ACOS computes the arc cosine function, with argument truncation.
 R8_ACOSH returns the inverse hyperbolic cosine of a number.
 R8_ADD adds two R8's.
 R8_AGM computes the arithmetic-geometric mean of A and B.
 R8_AINT truncates an R8 argument to an integer.
 R8_ASIN computes the arc sine function, with argument truncation.
 R8_ASINH returns the inverse hyperbolic sine of a number.
 R8_ATAN computes the inverse tangent of the ratio Y / X.
 R8_ATANH returns the inverse hyperbolic tangent of a number.
 R8_BIG returns a "big" R8.
 R8_CAS returns the "casine" of an R8.
 R8_CEILING rounds an R8 up to the nearest integral R8.
 R8_CHOOSE computes the combinatorial coefficient C(N,K).
 R8_CHOP chops an R8 to a given number of binary places.
 R8_COSD returns the cosine of an angle given in degrees.
 R8_COT returns the cotangent of an angle.
 R8_COTD returns the cotangent of an angle given in degrees.
 R8_CSC returns the cosecant of X.
 R8_CSCD returns the cosecant of an angle given in degrees.
 R8_CUBE_ROOT returns the cube root of an R8.
 R8_DEGREES converts an angle from radian to degree measure.
 R8_DIFF computes (X-Y) to a specified accuracy.
 R8_DIGIT returns a particular decimal digit of an R8.
 R8_DIVIDE_I4 returns an I4 fraction as an R8.
 R8_E returns the value of the base of the natural logarithm system.
 R8_EPSILON returns the R8 roundoff unit.
 R8_EPSILON_COMPUTE computes the R8 roundoff unit.
 R8_EXP computes the exponential function, avoiding overflow and underflow.
 R8_FACTORIAL computes the factorial of N.
 R8_FACTORIAL_VALUES returns values of the real factorial function.
 R8_FACTORIAL2 computes the double factorial function.
 R8_FACTORIAL2_VALUES returns values of the double factorial function.
 R8_FALL computes the falling factorial function [X]_N.
 R8_FALL_VALUES returns some values of the falling factorial function.
 R8_FLOOR rounds an R8 down to the nearest integral R8.
 R8_FRACTION uses real arithmetic on an integer ratio.
 R8_FRACTIONAL returns the fractional part of an R8.
 R8_GAMMA evaluates Gamma(X) for an R8.
 R8_GAMMA_LOG evaluates the logarithm of the gamma function.
 R8_HUGE returns a "huge" R8.
 R8_HYPOT returns the value of sqrt ( X^2 + Y^2 ).
 R8_IN_01 is TRUE if an R8 is in the range [0,1].
 R8_INSIGNIFICANT determines if an R8 is insignificant.
 R8_IS_INT determines if an R8 represents an integer value.
 R8_LOG_10 returns the logarithm base 10 of the absolute value of an R8.
 R8_LOG_2 returns the logarithm base 2 of the absolute value of an R8.
 R8_LOG_B returns the logarithm base B of an R8.
 R8_MANT computes the "mantissa" or "fraction part" of an R8.
 R8_MAX returns the maximum of two R8's.
 R8_MIN returns the minimum of two R8's.
 R8_MOD returns the remainder of R8 division.
 R8_MODP returns the nonnegative remainder of R8 division.
 R8_MOP returns the I-th power of -1 as an R8 value.
 R8_NINT returns the nearest integer to an R8.
 R8_NORMAL_01 samples the standard normal probability distribution.
 R8_NORMAL_AB returns a scaled pseudonormal R8.
 R8_PI returns the value of PI as an R8.
 R8_PI_SQRT returns the square root of PI as an R8.
 R8_POWER computes an integer power of an R8.
 R8_POWER_FAST computes the P-th power of R, for real R and integer P.
 R8_PRINT prints an R8.
 R8_PYTHAG computes sqrt ( A*A + B*B ), avoiding overflow and underflow.
 R8_RADIANS converts an angle from degree to radian measure.
 R8_REVERSE_BYTES reverses the bytes in an R8.
 R8_RISE computes the rising factorial function [X]^N.
 R8_RISE_VALUES returns some values of the rising factorial function.
 R8_ROUND rounds an R8 to the nearest integral value.
 R8_ROUND_I4 rounds an R8, returning an I4.
 R8_ROUND2 rounds an R8 in base 2.
 R8_ROUNDB rounds an R8 in a given base.
 R8_ROUNDX rounds an R8 in base 10.
 R8_SECD returns the secant of an angle given in degrees.
 R8_SECH evaluates the hyperbolic secant, while avoiding COSH overflow.
 R8_SIGN returns the sign of an R8.
 R8_SIGN3 returns the three-way sign of an R8.
 R8_SIGN_CHAR returns a character indicating the sign of an R8.
 R8_SIGN_MATCH is TRUE if two R8's are of the same sign.
 R8_SIGN_MATCH_STRICT is TRUE if two R8's are of the same strict sign.
 R8_SIGN_OPPOSITE is TRUE if two R8's are not of the same sign.
 R8_SIGN_OPPOSITE_STRICT is TRUE if two R8's are strictly of opposite sign.
 R8_SIGN2 returns the first argument with the sign of the second.
 R8_SINCOS_SUM simplifies a*sin(cx)+b*cos(cx).
 R8_SIND returns the sine of an angle given in degrees.
 R8_SQRT_I4 returns the square root of an I4 as an R8.
 R8_SUM returns the sum of two R8's.
 R8_SWAP switches two R8's.
 R8_SWAP3 swaps three R8's.
 R8_TAND returns the tangent of an angle given in degrees.
 R8_TINY returns a "tiny" R8.
 R8_TO_DHMS converts an R8 day value into days, hours, minutes, seconds.
 R8_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
 R8_TO_R8_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
 R8_UNIFORM_01 returns a unit pseudorandom R8.
 R8_UNIFORM_AB returns a scaled pseudorandom R8.
 R8_UNSWAP3 unswaps three R8's.
 R8_WALSH_1D evaluates the Walsh function of a real scalar argument.
 R8_WRAP forces an R8 to lie between given limits by wrapping.
 R82_DIST_L2 returns the L2 distance between a pair of R82's.
 R82_PRINT prints an R82.
 R82_UNIFORM_AB returns a random R82 value in a given range.
 R82COL_PRINT_PART prints "part" of an R82COL.
 R82POLY2_PRINT prints a second order polynomial in two variables.
 R82POLY2_TYPE analyzes a second order polynomial in two variables.
 R82POLY2_TYPE_PRINT prints the meaning of the output from R82POLY2_TYPE.
 R82ROW_MAX returns the maximum value in an R82ROW.
 R82ROW_MIN returns the minimum value in an R82ROW.
 R82ROW_ORDER_TYPE finds if an R82ROW is (non)strictly ascending/descending.
 R82ROW_PART_QUICK_A reorders an R82ROW as part of a quick sort.
 R82ROW_PERMUTE permutes an R82ROW in place.
 R82ROW_PRINT prints an R82ROW.
 R82ROW_PRINT_PART prints "part" of an R82ROW.
 R82ROW_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82ROW.
 R82ROW_SORT_QUICK_A ascending sorts an R82ROW using quick sort.
 R83_NORM returns the Euclidean norm of an R83.
 R83COL_PRINT_PART prints "part" of an R83COL.
 R83ROW_MAX returns the maximum value in an R83ROW.
 R83ROW_MIN returns the minimum value in an R83ROW.
 R83ROW_PART_QUICK_A reorders an R83ROW as part of a quick sort.
 R83ROW_PRINT_PART prints "part" of an R83ROW.
 R83ROW_SORT_QUICK_A ascending sorts an R83ROW using quick sort.
 R8BLOCK_DELETE frees memory associated with an R8BLOCK.
 R8BLOCK_EXPAND_LINEAR linearly interpolates new data into a 3D block.
 R8BLOCK_NEW allocates a new R8BLOCK.
 R8BLOCK_PRINT prints an R8BLOCK block (a 3D matrix).
 R8BLOCK_ZEROS_NEW returns a new zeroed R8BLOCK.
 R8CMAT_DELETE frees memory associated with an R8CMAT.
 R8CMAT_NEW allocates a new R8CMAT.
 R8CMAT_PRINT prints an R8CMAT.
 R8CMAT_PRINT_SOME prints some of an R8CMAT.
 R8CMAT_TO_R8MAT_NEW copies data from an R8CMAT to an R8MAT.
 R8CMAT_ZEROS_NEW allocates and zeros a new R8CMAT.
 R8INT_TO_R8INT maps one R8 interval to another.
 R8INT_TO_I4INT maps an R8 interval to an integer interval.
 R8MAT_ADD computes C = alpha * A + beta * B for R8MAT's.
 R8MAT_ADD_NEW computes C = alpha * A + beta * B for R8MAT's.
 R8MAT_AMAX returns the maximum absolute value entry of an R8MAT.
 R8MAT_BORDER_ADD adds a "border" to an R8MAT.
 R8MAT_BORDER_CUT cuts the "border" of an R8MAT.
 R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric R8MAT.
 R8MAT_CHOLESKY_FACTOR_UPPER: upper Cholesky factor of a symmetric R8MAT.
 R8MAT_CHOLESKY_INVERSE computes the inverse of a symmetric matrix.
 R8MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
 R8MAT_CHOLESKY_SOLVE_UPPER solves Cholesky factored linear system A * x = b.
 R8MAT_COPY copies one R8MAT to another.
 R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
 R8MAT_COVARIANCE computes the sample covariance of a set of vector data.
 R8MAT_DET computes the determinant of an R8MAT.
 R8MAT_DET_2D computes the determinant of a 2 by 2 R8MAT.
 R8MAT_DET_3D computes the determinant of a 3 by 3 R8MAT.
 R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
 R8MAT_DET_5D computes the determinant of a 5 by 5 R8MAT.
 R8MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R8MAT.
 R8MAT_DIAG_ADD_VECTOR adds a vector to the diagonal of an R8MAT.
 R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
 R8MAT_DIAG_GET_VECTOR_NEW gets the value of the diagonal of an R8MAT.
 R8MAT_DIAG_SET_SCALAR sets the diagonal of an R8MAT to a scalar value.
 R8MAT_DIAG_SET_VECTOR sets the diagonal of an R8MAT to a vector.
 R8MAT_DIAGONAL_NEW returns a diagonal matrix.
 R8MAT_DIFF_FROBENIUS returns the Frobenius norm of the difference of R8MAT's.
 R8MAT_EXPAND_LINEAR linearly interpolates new data into an R8MAT.
 R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.
 R8MAT_FLIP_COLS_NEW makes a new copy of an R8MAT with reversed column order.
 R8MAT_FLIP_ROWS_NEW makes a new copy of an R8MAT with reversed row order.
 R8MAT_FS factors and solves a system with one right hand side.
 R8MAT_FS_NEW factors and solves a system with one right hand side.
 R8MAT_FSS factors and solves a system with multiple right hand sides.
 R8MAT_FSS_NEW factors and solves a system with multiple right hand sides.
 R8MAT_GIVENS_POST computes the Givens postmultiplier rotation matrix.
 R8MAT_GIVENS_PRE computes the Givens premultiplier rotation matrix.
 R8MAT_HESS approximates a Hessian matrix via finite differences.
 R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
 R8MAT_HOUSE_AXH_NEW computes A*H where H is a compact Householder matrix.
 R8MAT_HOUSE_FORM constructs a Householder matrix from its compact form.
 R8MAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.
 R8MAT_HOUSE_POST computes a Householder post-multiplier matrix.
 R8MAT_HOUSE_PRE computes a Householder pre-multiplier matrix.
 R8MAT_IDENTITY sets the square matrix A to the identity.
 R8MAT_IDENTITY_NEW returns an identity matrix.
 R8MAT_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].
 R8MAT_INDICATOR_NEW sets up an "indicator" R8MAT.
 R8MAT_INSIGNIFICANT determines if an R8MAT is insignificant.
 R8MAT_INVERSE_2D inverts a 2 by 2 matrix using Cramer's rule.
 R8MAT_INVERSE_3D inverts a 3 by 3 matrix using Cramer's rule.
 R8MAT_INVERSE_4D inverts a 4 by 4 matrix using Cramer's rule.
 R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
 R8MAT_IS_SYMMETRIC checks an R8MAT for symmetry.
 R8MAT_JAC estimates a dense jacobian matrix of the function FX.
 R8MAT_KRONECKER computes the Kronecker product of two R8MAT's.
 R8MAT_L_INVERSE inverts a lower triangular R8MAT.
 R8MAT_L_PRINT prints a lower triangular R8MAT.
 R8MAT_L_SOLVE solves a lower triangular linear system.
 R8MAT_L1_INVERSE inverts a unit lower triangular R8MAT.
 R8MAT_LT_SOLVE solves a transposed lower triangular linear system.
 R8MAT_LU computes the LU factorization of a rectangular R8MAT.
 R8MAT_MAX returns the maximum entry of an R8MAT.
 R8MAT_MAX_INDEX returns the location of the maximum entry of an R8MAT.
 R8MAT_MAXCOL_MINROW gets the maximum column minimum row of an M by N matrix.
 R8MAT_MAXROW_MINCOL gets the maximum row minimum column of an M by N matrix.
 R8MAT_MEAN returns the mean of an R8MAT.
 R8MAT_MIN returns the minimum entry of an R8MAT.
 R8MAT_MIN_INDEX returns the location of the minimum entry of an R8MAT.
 R8MAT_MINCOL_MAXROW gets the minimum column maximum row of an M by N matrix.
 R8MAT_MINROW_MAXCOL gets the minimum row maximum column of an M by N matrix.
 R8MAT_MINVM computes inverse(A) * B for R8MAT's.
 R8MAT_MINVM_NEW returns inverse(A) * B for R8MAT's.
 R8MAT_MM multiplies two matrices.
 R8MAT_MM_NEW multiplies two matrices.
 R8MAT_MMT_NEW computes C = A * B'.
 R8MAT_MTM_NEW computes C = A' * B.
 R8MAT_MTV multiplies a transposed matrix times a vector.
 R8MAT_MTV_NEW multiplies a transposed matrix times a vector.
 R8MAT_MV multiplies a matrix times a vector.
 R8MAT_MV_NEW multiplies a matrix times a vector.
 R8MAT_NINT rounds the entries of an R8MAT.
 R8MAT_NONZEROS returns the number of nonzeros in an R8MAT.
 R8MAT_NORM_EIS returns the EISPACK norm of an R8MAT.
 R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
 R8MAT_NORM_FRO_AFFINE returns the Frobenius norm of an R8MAT difference.
 R8MAT_NORM_L1 returns the matrix L1 norm of an R8MAT.
 R8MAT_NORM_L2 returns the matrix L2 norm of an R8MAT.
 R8MAT_NORM_LI returns the matrix L-oo norm of an R8MAT.
 R8MAT_NORMAL_01_NEW returns a unit pseudonormal R8MAT.
 R8MAT_NULLSPACE computes the nullspace of a matrix.
 R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.
 R8MAT_ORTH_UNIFORM_NEW returns a random orthogonal matrix.
 R8MAT_PLOT "plots" an R8MAT.
 R8MAT_PLOT_SYMBOL returns a symbol for entries of an R8MAT.
 R8MAT_POLY_CHAR computes the characteristic polynomial of an R8MAT.
 R8MAT_POWER computes a nonnegative power of an R8MAT.
 R8MAT_POWER_METHOD applies the power method to a matrix.
 R8MAT_PRINT prints an R8MAT.
 R8MAT_PRINT_SOME prints some of an R8MAT.
 R8MAT_REF computes the row echelon form of a matrix.
 R8MAT_RMS returns the RMS norm of an R8MAT.
 R8MAT_ROW_COPY copies a vector into a row of an R8MAT.
 R8MAT_RREF computes the reduced row echelon form of a matrix.
 R8MAT_SCALE multiplies an R8MAT by a scalar.
 R8MAT_SIGNIFICANT determines if an R8MAT is significant compared to another.
 R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
 R8MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
 R8MAT_SOLVE_3D solves a 3 by 3 linear system using Cramer's rule.
 R8MAT_SOLVE2 computes the solution of an N by N linear system.
 R8MAT_SUB_NEW computes C = A - B.
 R8MAT_SUM returns the sum of an R8MAT.
 R8MAT_SYMM_EIGEN returns a symmetric matrix with given eigensystem.
 R8MAT_SYMM_JACOBI applies Jacobi eigenvalue iteration to a symmetric matrix.
 R8MAT_TO_R8CMAT_NEW copies data from an R8MAT to an R8CMAT.
 R8MAT_TO_R8PLU factors a general matrix.
 R8MAT_TO_R8RMAT copies data from an R8MAT to an R8RMAT.
 R8MAT_TRACE computes the trace of an R8MAT.
 R8MAT_TRANSPOSE_IN_PLACE transposes a square R8MAT in place.
 R8MAT_TRANSPOSE_NEW returns the transpose of an R8MAT.
 R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
 R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
 R8MAT_U_INVERSE inverts an upper triangular R8MAT.
 R8MAT_U_SOLVE solves an upper triangular linear system.
 R8MAT_U1_INVERSE inverts a unit upper triangular R8MAT.
 R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
 R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
 R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
 R8MAT_UNIFORM_AB_NEW returns a new scaled pseudorandom R8MAT.
 R8MAT_UNIFORM_ABVEC returns a scaled pseudorandom R8MAT.
 R8MAT_UNIFORM_ABVEC_NEW returns a new scaled pseudorandom R8MAT.
 R8MAT_UT_SOLVE solves a transposed upper triangular linear system.
 R8MAT_VAND2 returns the N by N row Vandermonde matrix A.
 R8MAT_VTMV multiplies computes the scalar x' * A * y.
 R8MAT_ZEROS zeroes an R8MAT.
 R8MAT_ZEROS_NEW returns a new zeroed R8MAT.
 R8PLU_DET computes the determinant of a real PLU matrix.
 R8PLU_INVERSE computes the inverse of a real PLU matrix.
 R8PLU_MUL computes A * x using the PLU factors of A.
 R8PLU_SOL solves a linear system A*x=b from the PLU factors.
 R8PLU_TO_R8MAT recovers the matrix A that was factored by R8MAT_TO_R8PLU.
 R8POLY_DEGREE returns the degree of a polynomial.
 R8POLY_DERIV returns the derivative of a polynomial.
 R8POLY_LAGRANGE_0 evaluates the Lagrange factor at a point.
 R8POLY_LAGRANGE_1 evaluates the first derivative of the Lagrange factor.
 R8POLY_LAGRANGE_2 evaluates the second derivative of the Lagrange factor.
 R8POLY_LAGRANGE_COEF returns the coefficients of a Lagrange polynomial.
 R8POLY_LAGRANGE_FACTOR evaluates the polynomial Lagrange factor at a point.
 R8POLY_LAGRANGE_VAL evaluates the IPOL-th Lagrange polynomial.
 R8POLY_ORDER returns the order of a polynomial.
 R8POLY_PRINT prints out a polynomial.
 R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
 R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
 R8POLY_VALUES_HORNER evaluates a polynomial using Horner's method.
 R8POLY_VALUE_2D evaluates a polynomial in 2 variables, X and Y.
 R8POLY2_EX finds the extremal point of a parabola determined by three points.
 R8POLY2_EX2 finds the extremal point of a parabola determined by three points.
 R8POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
 R8POLY2_VAL evaluates a parabola defined by three data values.
 R8POLY2_VAL2 evaluates a parabolic function through 3 points in a table.
 R8PP_DELETE frees the memory set aside by R8PP_NEW.
 R8PP_NEW allocates a new R8PP.
 R8R8_COMPARE compares two R8R8's.
 R8R8_PRINT prints an R8R8.
 R8R8R8_COMPARE compares two R8R8R8's.
 R8R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8R8 value in an indexed sorted list.
 R8R8R8VEC_INDEX_SEARCH searches for an R8R8R8 value in an indexed sorted list.
 R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8 value in an indexed sorted list.
 R8R8VEC_INDEX_SEARCH searches for an R8R8 value in an indexed sorted list.
 R8RMAT_COPY_NEW makes a new copy of an R8RMAT .
 R8RMAT_DELETE frees memory associated with an R8RMAT.
 R8RMAT_FS_NEW factors and solves an R8RMAT system with one right hand side.
 R8RMAT_NEW allocates a new R8RMAT.
 R8RMAT_PRINT prints an R8RMAT.
 R8RMAT_PRINT_SOME prints some of an R8RMAT.
 R8RMAT_TO_R8MAT copies data from an R8RMAT to an R8MAT.
 R8RMAT_ZEROS allocates and zeroes a new R8RMAT.
 R8ROW_COMPARE compares two rows in an R8ROW.
 R8ROW_MAX returns the row maximums of an R8ROW.
 R8ROW_MEAN returns the row means of an R8ROW.
 R8ROW_MIN returns the row minimums of an R8ROW.
 R8ROW_PRINT prints an R8ROW.
 R8ROW_PRINT_SOME prints some of an R8ROW.
 R8ROW_REVERSE reverses the order of the rows of an R8MAT.
 R8ROW_RUNNING_AVERAGE computes the running averages of an R8ROW.
 R8ROW_RUNNING_SUM computes the running sum of an R8ROW.
 R8ROW_SORT_HEAP_A ascending heapsorts an R8ROW.
 R8ROW_SUM returns the sums of the rows of an R8ROW.
 R8ROW_SWAP swaps two rows of an R8ROW.
 R8ROW_TO_R8VEC converts an R8ROW into an R8VEC.
 R8ROW_UNIFORM_AB_NEW returns a new scaled pseudorandom R8ROW.
 R8ROW_VARIANCE returns the variances of an R8ROW.
 R8SLMAT_PRINT prints a strict lower triangular R8MAT.
 R8VEC_01_TO_AB shifts and rescales data to lie within given bounds.
 R8VEC_AB_TO_01 shifts and rescales data to lie within [0,1].
 R8VEC_AB_TO_CD shifts and rescales data to lie within a given pair of bounds.
 R8VEC_ADD adds one R8VEC to another.
 R8VEC_ALL_NONPOSITIVE: ( all ( A <= 0 ) ) for R8VEC's.
 R8VEC_AMAX returns the maximum absolute value in an R8VEC.
 R8VEC_AMAX_INDEX returns the index of the maximum absolute value in an R8VEC.
 R8VEC_AMIN returns the minimum absolute value in an R8VEC.
 R8VEC_AMIN_INDEX returns the index of the minimum absolute value in an R8VEC.
 R8VEC_ANY_NEGATIVE: ( any ( A < 0 ) ) for R8VEC's.
 R8VEC_ANY_NONZERO: ( any A nonzero ) for R8VEC's.
 R8VEC_ANY_NORMAL returns some normal vector to V1.
 R8VEC_ASCENDS determines if an R8VEC is (weakly) ascending.
 R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
 R8VEC_ASUM sums the absolute values of the entries of an R8VEC.
 R8VEC_BIN computes bins based on a given R8VEC.
 R8VEC_BRACKET searches a sorted array for successive brackets of a value.
 R8VEC_BRACKET2 searches a sorted array for successive brackets of a value.
 R8VEC_BRACKET3 finds the interval containing or nearest a given value.
 R8VEC_BRACKET4 finds the interval containing or nearest a given value.
 R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
 R8VEC_BRACKET6 brackets data between successive entries of a sorted R8VEC.
 R8VEC_CHEBYSPACE_NEW creates a vector of Chebyshev spaced values in [A,B].
 R8VEC_CHEBY1SPACE_NEW creates Type 1 Chebyshev spaced values in [A,B].
 R8VEC_CHEBY2SPACE_NEW creates Type 2 Chebyshev spaced values in [A,B].
 R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC.
 R8VEC_COMPARE compares two R8VEC's.
 R8VEC_CONCATENATE concatenates two R8VEC's.
 R8VEC_CONCATENATE_NEW concatenates two R8VEC's.
 R8VEC_CONVOLUTION returns the convolution of two R8VEC's.
 R8VEC_CONVOLUTION_CIRC returns the discrete circular convolution of two R8VEC's.
 R8VEC_COPY copies an R8VEC.
 R8VEC_COPY_NEW copies an R8VEC to a new R8VEC.
 R8VEC_CORRELATION returns the correlation of two R8VEC's.
 R8VEC_COVAR computes the covariance of two vectors.
 R8VEC_CROSS_PRODUCT_2D finds the cross product of a pair of R8VEC's in 2D.
 R8VEC_CROSS_PRODUCT_AFFINE_2D finds the affine cross product in 2D.
 R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.
 R8VEC_CROSS_PRODUCT_AFFINE_3D computes the affine cross product in 3D.
 R8VEC_CUM_NEW computes the cumulutive sums of an R8VEC.
 R8VEC_CUM0_NEW computes the cumulutive sums of an R8VEC.
 R8VEC_DIF computes coefficients for estimating the N-th derivative.
 R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
 R8VEC_DIFF_NORM_L1 returns the L1 norm of the difference of R8VEC's.
 R8VEC_DIFF_NORM_L2 returns the L2 norm of the difference of R8VEC's.
 R8VEC_DIFF_NORM_LI returns the L-oo norm of the difference of R8VEC's.
 R8VEC_DIFF_NORM_SQUARED: square of the L2 norm of the difference of R8VEC's.
 R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
 R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
 R8VEC_DISTANCE returns the Euclidean distance between two R8VEC's.
 R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
 R8VEC_DIVIDE divides an R8VEC by a nonzero scalar.
 R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
 R8VEC_DOT_PRODUCT_AFFINE computes the affine dot product.
 R8VEC_ENTROPY computes the entropy of an R8VEC.
 R8VEC_EQ is true if every pair of entries in two R8VEC's is equal.
 R8VEC_EVEN returns an R8VEC of values evenly spaced between ALO and AHI.
 R8VEC_EVEN_NEW returns an R8VEC of values evenly spaced between ALO and AHI.
 R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
 R8VEC_EVEN2 linearly interpolates new numbers into an R8VECa.
 R8VEC_EVEN2_SELECT returns the I-th of N evenly spaced midpoint values.
 R8VEC_EVEN3 evenly interpolates new data into an R8VEC.
 R8VEC_EXPAND_LINEAR linearly interpolates new data into an R8VEC.
 R8VEC_EXPAND_LINEAR2 linearly interpolates new data into an R8VEC.
 R8VEC_FIRST_INDEX indexes the first occurrence of values in an R8VEC.
 R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.
 R8VEC_FRACTION returns the fraction parts of an R8VEC.
 R8VEC_GT == ( A1 > A2 ) for two R8VEC's.
 R8VEC_HEAP_A reorders an R8VEC into a ascending heap.
 R8VEC_HEAP_D reorders an R8VEC into a descending heap.
 R8VEC_HISTOGRAM histograms an R8VEC.
 R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
 R8VEC_I4VEC_DOT_PRODUCT computes the dot product of an R8VEC and an I4VEC.
 R8VEC_IN_01 is TRUE if the entries of an R8VEC are in the range [0,1].
 R8VEC_IN_AB is TRUE if the entries of an R8VEC are in the range [A,B].
 R8VEC_INDEX_DELETE_ALL deletes all occurrences of a value from an indexed sorted list.
 R8VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted list.
 R8VEC_INDEX_DELETE_ONE deletes one copy of a value from an indexed sorted list.
 R8VEC_INDEX_INSERT inserts a value in an indexed sorted list.
 R8VEC_INDEX_INSERT_UNIQUE inserts a unique value in an indexed sorted list.
 R8VEC_INDEX_ORDER sorts an R8VEC using an index vector.
 R8VEC_INDEX_SEARCH searches for a value in an indexed sorted R8VEC.
 R8VEC_INDEX_SORT_UNIQUE creates a sort index for an R8VEC.
 R8VEC_INDEX_SORTED_RANGE: search index sorted vector for elements in a range.
 R8VEC_INDEXED_HEAP_D creates a descending heap from an indexed R8VEC.
 R8VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed R8VEC.
 R8VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed R8VEC.
 R8VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed R8VEC.
 R8VEC_INDICATOR0 sets an R8VEC to the indicator vector (0,1,2,...)
 R8VEC_INDICATOR0_NEW sets an R8VEC to the indicator vector {0,1,2,...}.
 R8VEC_INDICATOR1 sets an R8VEC to the indicator vector (1,2,3,...)
 R8VEC_INDICATOR1_NEW sets an R8VEC to the indicator vector {1,2,3,...}.
 R8VEC_INSERT inserts a value into an R8VEC.
 R8VEC_INSIGNIFICANT determines if an R8VEC is insignificant.
 R8VEC_IS_INT is TRUE if an R8VEC is integral.
 R8VEC_IS_NONNEGATIVE is true if all entries in an R8VEC are nonnegative.
 R8VEC_IS_ZERO is true if the entries in an R8VEC are all zero.
 R8VEC_LEGENDRE_NEW creates a vector of Chebyshev spaced values.
 R8VEC_LINSPACE creates a vector of linearly spaced values.
 R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
 R8VEC_LINSPACE2_NEW creates a vector of linearly spaced values.
 R8VEC_LT == ( A1 < A2 ) for two R8VEC's.
 R8VEC_MASK_PRINT prints a masked R8VEC.
 R8VEC_MAX returns the value of the maximum element in an R8VEC.
 R8VEC_MAX_ABS_INDEX returns the index of the maximum absolute value in an R8VEC.
 R8VEC_MAX_INDEX returns the index of the maximum value in an R8VEC.
 R8VEC_MEAN returns the mean of an R8VEC.
 R8VEC_MEAN_GEOMETRIC returns the geometric mean of an R8VEC.
 R8VEC_MEDIAN returns the median of an unsorted R8VEC.
 R8VEC_MESH_2D creates a 2D mesh from X and Y vectors.
 R8VEC_MIDSPACE_NEW creates a vector of linearly spaced values.
 R8VEC_MIN returns the value of the minimum element in an R8VEC.
 R8VEC_MIN_INDEX returns the index of the minimum value in an R8VEC.
 R8VEC_MIN_POS returns the minimum positive value of an R8VEC.
 R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
 R8VEC_NEGATIVE_STRICT: all entries of R8VEC are strictly negative.
 R8VEC_NINT rounds the entries of an R8VEC.
 R8VEC_NINT_NEW rounds the entries of an R8VEC.
 R8VEC_NORM returns the L2 norm of an R8VEC.
 R8VEC_NORM_AFFINE returns the affine L2 norm of an R8VEC.
 R8VEC_NORM_L0 returns the l0 "norm" of an R8VEC.
 R8VEC_NORM_L1 returns the L1 norm of an R8VEC.
 R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
 R8VEC_NORM_LI returns the L-oo norm of an R8VEC.
 R8VEC_NORM_LP returns the LP norm of an R8VEC.
 R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
 R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
 R8VEC_NORMAL_AB_NEW returns a scaled pseudonormal R8VEC.
 R8VEC_NORMALIZE normalizes an R8VEC.
 R8VEC_NORMALIZE_L1 normalizes an R8VEC to have unit sum.
 R8VEC_NORMSQ returns the squared L2 norm of an R8VEC.
 R8VEC_NORMSQ_AFFINE returns the squared affine L2 norm of an R8VEC.
 R8VEC_ONES_NEW creates a vector of 1's.
 R8VEC_ORDER_TYPE determines if an R8VEC is (non)strictly ascending/descending.
 R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.
 R8VEC_PERMUTE applies a 0-based permutation to an R8VEC.
 R8VEC_PERMUTE_CYCLIC performs a cyclic permutation of an R8VEC.
 R8VEC_PERMUTE_UNIFORM randomly permutes an R8VEC.
 R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
 R8VEC_POSITIVE_STRICT: all entries of R8VEC are strictly positive.
 R8VEC_PRINT prints an R8VEC.
 R8VEC_PRINT_16 prints an R8VEC to 16 decimal places.
 R8VEC_PRINT_PART prints "part" of an R8VEC.
 R8VEC_PRINT_SOME prints "some" of an R8VEC.
 R8VEC_PRODUCT returns the product of the entries of an R8VEC.
 R8VEC_RANGE finds the range of Y's within a restricted X range.
 R8VEC_RANGE_2 updates a range to include a new R8VEC
 R8VEC_REVERSE reverses the elements of an R8VEC.
 R8VEC_RMS returns the RMS norm of an R8VEC.
 R8VEC_ROTATE "rotates" the entries of an R8VEC in place.
 R8VEC_RUNNING_AVERAGE computes the running average of an R8VEC.
 R8VEC_RUNNING_SIGN3 computes the running threeway sign of an R8VEC.
 R8VEC_RUNNING_SUM computes the running sum of an R8VEC.
 R8VEC_SCALAR_TRIPLE_PRODUCT computes the scalar triple product.
 R8VEC_SCALE multiplies an R8VEC by a scale factor.
 R8VEC_SEARCH_BINARY_A searches an ascending sorted R8VEC.
 R8VEC_SHIFT performs a shift on an R8VEC.
 R8VEC_SHIFT_CIRCULAR performs a circular shift on an R8VEC.
 R8VEC_SORT_BUBBLE_A ascending sorts an R8VEC using bubble sort.
 R8VEC_SORT_BUBBLE_D descending sorts an R8VEC using bubble sort.
 R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.
 R8VEC_SORT_HEAP_D descending sorts an R8VEC using heap sort.
 R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC
 R8VEC_SORT_HEAP_INDEX_A_NEW: indexed heap ascending sort of an R8VEC
 R8VEC_SORT_HEAP_INDEX_D_NEW: indexed heap descending sort of an R8VEC.
 R8VEC_SORT_HEAP_INDEX_D_NEW: indexed heap descending sort of an R8VEC.
 R8VEC_SORT_HEAP_MASK_A: indexed heap ascending sort of a masked R8VEC.
 R8VEC_SORT_INSERT_A ascending sorts an R8VEC using an insertion sort.
 R8VEC_SORT_INSERT_INDEX_A ascending index sorts an R8VEC using insertion.
 R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.
 R8VEC_SORT_SHELL_A ascending sorts an R8VEC using Shell's sort.
 R8VEC_SORTED_MERGE_A merges two ascending sorted R8VEC's.
 R8VEC_SORTED_NEAREST returns the nearest element in a sorted R8VEC.
 R8VEC_SORTED_RANGE searches a sorted vector for elements in a range.
 R8VEC_SORTED_SPLIT "splits" a sorted R8VEC, given a splitting value.
 R8VEC_SORTED_UNDEX returns unique sorted indexes for a sorted R8VEC.
 R8VEC_SORTED_UNIQUE finds the unique elements in a sorted R8VEC.
 R8VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8VEC.
 R8VEC_SORTED_UNIQUE_HIST histograms unique elements of a sorted R8VEC.
 R8VEC_SPLIT "splits" an unsorted R8VEC based on a splitting value.
 R8VEC_STD returns the standard deviation of an R8VEC.
 R8VEC_STEP evaluates a unit step function.
 R8VEC_STUTTER makes a "stuttering" copy of an R8VEC.
 R8VEC_STUTTER_NEW makes a "stuttering" copy of an R8VEC.
 R8VEC_SUM returns the sum of an R8VEC.
 R8VEC_SWAP swaps the entries of two R8VEC's.
 R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
 R8VEC_UNDEX returns unique sorted indexes for an R8VEC.
 R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
 R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
 R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
 R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.
 R8VEC_UNIFORM_ABVEC returns a scaled pseudorandom R8VEC.
 R8VEC_UNIFORM_ABVEC_NEW returns a scaled pseudorandom R8VEC.
 R8VEC_UNIFORM_UNIT_NEW generates a random unit vector.
 R8VEC_UNIQUE_COUNT counts the unique elements in an unsorted R8VEC.
 R8VEC_UNIQUE_INDEX indexes the unique occurrence of values in an R8VEC.
 R8VEC_VARIANCE returns the variance of an R8VEC.
 R8VEC_VECTOR_TRIPLE_PRODUCT computes the vector triple product.
 R8VEC_WRITE writes an R8VEC to a file.
 R8VEC_ZEROS zeroes an R8VEC.
 R8VEC_ZEROS_NEW creates and zeroes an R8VEC.
 R8VEC2_COMPARE compares two elements of an R8VEC2.
 R8VEC2_PRINT prints an R8VEC2.
 R8VEC2_PRINT_SOME prints "some" of an R8VEC2.
 R8VEC2_SORT_A ascending sorts an R8VEC2.
 R8VEC2_SORT_D descending sorts an R8VEC2.
 R8VEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC2.
 R8VEC2_SORTED_UNIQUE keeps the unique elements in an R8VEC2.
 R8VEC2_SORTED_UNIQUE_INDEX indexes unique elements in a sorted R8VEC2.
 R8VEC2_SUM_MAX_INDEX returns the index of the maximum sum of two R8VEC's.
 R8VEC3_PRINT prints a triple of real vectors.
 ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
 S_LEN_TRIM returns the length of a string to the last nonblank.
 SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
 TIMESTAMP prints the current YMDHMS date as a time stamp.
 */
using namespace::std;
void int_mat_print ( int m, int n, int **a, string title ); // Adapted by Gilberto Rouxinol
void int_print_some ( int m, int n, int **a, int ilo, int jlo, int ihi, int jhi, string title); // Adapted by Gilberto Rouxinol
void gamma_values ( int &n_data, double &x, double &fx );
void gamma_log_values ( int &n_data, double &x, double &fx );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
int i4_sign ( int i );
int i4_uniform_ab ( int a, int b, int &seed );
int i4_wrap ( int ival, int ilo, int ihi );
double i4int_to_r8int ( int imin, int imax, int i, double rmin, double rmax );
void i4mat_print ( int m, int n, int a[], string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void i4vec_copy ( int n, int a1[], int a2[] );
int *i4vec_indicator0_new ( int n );
int *i4vec_indicator1_new ( int n );
void i4vec_permute ( int n, int p[], int a[] );
void i4vec_print ( int n, int a[], string title );
void i4vec_transpose_print ( int n, int a[], string title );
void i4vec_zeros ( int n, int a[] );
int *i4vec_zeros_new ( int n );
void jacobi_eigenvalue ( int n, double a[], int it_max, double v[],
  double d[], int &it_num, int &rot_num );
double r8mat_is_eigen_right ( int n, int k, double a[], double x[],
  double lambda[] );
double *legendre_zeros ( int order );
bool perm0_check ( int n, int p[] );
int *perm0_uniform_new ( int n, int &seed );
bool perm1_check ( int n, int p[] );
int *perm1_uniform_new ( int n, int &seed );
double r8_abs ( double x );
double r8_acos ( double c );
double r8_acosh ( double c );
double r8_add ( double y, double x );
double r8_aint ( double x );
double r8_asin ( double s );
double r8_asinh ( double x );
double r8_atan ( double y, double x );
double r8_atanh ( double x );
double r8_big ( );
double r8_cas ( double x );
double r8_ceiling ( double x );
double r8_choose ( int n, int k );
double r8_chop ( int place, double x );
double r8_cosd ( double x );
double r8_cot ( double angle );
double r8_cotd ( double x );
double r8_csc ( double theta );
double r8_cscd ( double x );
double r8_cube_root ( double x );
double r8_degrees ( double radians );
double r8_diff ( double x, double y, int n );
int r8_digit ( double x, int idigit );
double r8_divide_i4 ( int i, int j );
double r8_e ( );
double r8_epsilon ( );
double r8_epsilon_compute ( );
double r8_exp ( double x );
double r8_factorial ( int n );
double r8_factorial_stirling ( int n );
void r8_factorial_values ( int &n_data, int &n, double &fn );
double r8_factorial2 ( int n );
void r8_factorial2_values ( int &n_data, int &n, double &f );
double r8_fall ( double x, int n );
void r8_fall_values ( int &n_data, double &x, int &n, double &f );
double r8_floor ( double x );
double r8_fraction ( int i, int j );
double r8_fractional ( double x );
double r8_gamma ( double x );
double r8_gamma_log ( double x );
double r8_huge ( );
double r8_hypot ( double a, double b );
bool r8_is_in_01 ( double a );
bool r8_is_inf ( double r );
bool r8_is_insignificant ( double r, double s );
bool r8_is_integer ( double r );
bool r8_is_nan ( double r );
double r8_log_10 ( double x );
double r8_log_2 ( double x );
double r8_log_b ( double x, double b );
void r8_mant ( double x, int &s, double &r, int &l );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_mod ( double x, double y );
double r8_modp ( double x, double y );
double r8_mop ( int i );
int r8_nint ( double x );
double r8_normal_01 ( int &seed );
double r8_normal_ab ( double a, double b, int &seed );
double r8_nth_root ( double x, int n );
double r8_pi ( );
double r8_pi_sqrt ( );
double r8_power ( double r, int p );
double r8_power_fast ( double r, int p, int &mults );
void r8_print ( double r, string title );
double r8_radians ( double degrees );
double r8_reverse_bytes ( double x );
double r8_rise ( double x, int n );
void r8_rise_values ( int &n_data, double &x, int &n, double &f );
double r8_round ( double x );
int r8_round_i4 ( double x );
double r8_round2 ( int nplace, double x );
double r8_roundb ( int base, int nplace, double x );
double r8_roundx ( int nplace, double x );
double r8_secd ( double x );
double r8_sech ( double x );
double r8_sigmoid ( double x );
double r8_sign ( double x );
double r8_sign3 ( double x );
char r8_sign_char ( double x );
bool r8_sign_match ( bool r1, bool r2 );
bool r8_sign_match_strict ( bool r1, bool r2 );
bool r8_sign_opposite ( double r1, double r2 );
bool r8_sign_opposite_strict ( double r1, double r2 );
double r8_sign2 ( double x, double y );
void r8_sincos_sum ( double a, double b, double &d, double &e, double &f );
double r8_sind ( double x );
double r8_softplus ( double x );
double r8_sqrt_i4 ( int i );
double r8_sum ( double x, double y );
void r8_swap ( double &x, double &y );
void r8_swap3 ( double &x, double &y, double &z );
double r8_tand ( double x );
double r8_tiny ( );
void r8_to_dhms ( double r, int &d, int &h, int &m, int &s );
int r8_to_i4 ( double xmin, double xmax, double x, int ixmin, int ixmax );
double r8_to_r8_discrete ( double r, double rmin, double rmax, int nr );
double r8_uniform_ab ( double b, double c, int &seed );
double r8_uniform_01 ( int &seed );
void r8_unswap3 ( double &x, double &y, double &z );
double r8_walsh_1d ( double x, int digit );
double r8_wrap ( double r, double rlo, double rhi );
double r82_dist_l2 ( double a1[2], double a2[2] );
void r82_print ( double a[2], string title );
void r82_uniform_ab ( double b, double c, int &seed, double r[] );
void r82col_print_part ( int n, double a[], int max_print, string title );
double *r82row_max ( int n, double a[] );
double *r82row_min ( int n, double a[] );
int r82row_order_type ( int n, double a[] );
void r82row_part_quick_a ( int n, double a[], int &l, int &r );
void r82row_permute ( int n, int p[], double a[] );
void r82row_print ( int n, double a[], string title );
void r82row_print_part ( int n, double a[], int max_print, string title );
int *r82row_sort_heap_index_a ( int n, double a[] );
void r82row_sort_quick_a ( int n, double a[] );
double r83_norm ( double x, double y, double z );
void r83col_print_part ( int n, double a[], int max_print, string title );
double *r83row_max ( int n, double a[] );
double *r83row_min ( int n, double a[] );
void r83row_part_quick_a ( int n, double a[], int &l, int &r );
void r83row_print_part ( int n, double a[], int max_print, string title );
void r83row_sort_quick_a ( int n, double a[] );
void r8block_delete ( int l, int m, int n, double ***a );
double *r8block_expand_linear ( int l, int m, int n, double x[], int lfat, int mfat,
  int nfat );
double ***r8block_new ( int l, int m, int n );
void r8block_print ( int l, int m, int n, double a[], string title );
double *r8block_zeros_new ( int l, int m, int n );
void r8cmat_delete ( int m, int n, double **a );
double **r8cmat_new ( int m, int n );
void r8cmat_print ( int m, int n, double **a, string title );
void r8cmat_print_some ( int m, int n, double **a, int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8cmat_to_r8mat_new ( int m, int n, double **a );
double **r8cmat_zeros_new ( int m, int n );
double r8int_to_r8int ( double rmin, double rmax, double r, double r2min,
  double r2max );
int r8int_to_i4int ( double rmin, double rmax, double r, int imin, int imax );
void r8mat_add ( int m, int n, double alpha, double a[], double beta, 
  double b[], double c[] );
double *r8mat_add_new ( int m, int n, double alpha, double a[], double beta, 
  double b[] );
double r8mat_amax ( int m, int n, double a[] );
double *r8mat_border_add ( int m, int n, double table[] );
double *r8mat_border_cut ( int m, int n, double table[] );
double *r8mat_cholesky_factor ( int n, double a[], int &flag );
double *r8mat_cholesky_factor_upper ( int n, double a[], int &flag );
void r8mat_cholesky_inverse ( int n, double a[] );
double *r8mat_cholesky_solve ( int n, double a[], double b[] );
double *r8mat_cholesky_solve_upper ( int n, double r[], double b[] );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
double *r8mat_copy_new ( int m, int n, double a1[] );
double *r8mat_covariance ( int m, int n, double x[] );
double r8mat_det ( int n, double a[] );
double r8mat_det_2d ( double a[] );
double r8mat_det_3d ( double a[] );
double r8mat_det_4d ( double a[] );
double r8mat_det_5d ( double a[] );
void r8mat_diag_add_scalar ( int n, double a[], double s );
void r8mat_diag_add_vector ( int n, double a[], double v[] );
void r8mat_diag_get_vector ( int n, double a[], double v[] );
double *r8mat_diag_get_vector_new ( int n, double a[] );
void r8mat_diag_set_scalar ( int n, double a[], double s );
void r8mat_diag_set_vector ( int n, double a[], double v[] );
double *r8mat_diagonal_new ( int n, double diag[] );
double r8mat_diff_frobenius ( int m, int n, double a[], double b[] );
double *r8mat_expand_linear ( int m, int n, double x[], int mfat, int nfat );
double *r8mat_expand_linear2 ( int m, int n, double a[], int m2, int n2 );
double *r8mat_flip_cols_new ( int m, int n, double a[] );
double *r8mat_flip_rows_new ( int m, int n, double a[] );
void r8mat_fs ( int n, double a[], double x[] );
double *r8mat_fs_new ( int n, double a[], double b[] );
void r8mat_fss ( int n, double a[], int nb, double x[] );
double *r8mat_fss_new ( int n, double a[], int nb, double b[] );
double *r8mat_givens_post ( int n, double a[], int row, int col );
double *r8mat_givens_pre ( int n, double a[], int row, int col );
double *r8mat_hess ( double (*fx )( int n, double x[] ), int n, double x[] );
void r8mat_house_axh ( int n, double a[], double v[] );
double *r8mat_house_axh_new ( int n, double a[], double v[] );
double *r8mat_house_form ( int n, double v[] );
double *r8mat_house_hxa ( int n, double a[], double v[] );
double *r8mat_house_post ( int n, double a[], int row, int col );
double *r8mat_house_pre ( int n, double a[], int row, int col );
void r8mat_identity ( int n, double a[] );
double *r8mat_identity_new ( int n );
double *r8mat_indicator_new ( int m, int n );
double *r8mat_inverse_2d ( double a[] );
double *r8mat_inverse_3d ( double a[] );
double *r8mat_inverse_4d ( double a[] );
bool r8mat_is_binary ( int m, int n, double a[] );
double r8mat_is_identity ( int n, double a[] );
bool r8mat_is_in_01 ( int m, int n, double a[] );
bool r8mat_is_insignificant ( int m, int n, double r[], double s[] );
bool r8mat_is_integer ( int m, int n, double a[] );
bool r8mat_is_significant ( int m, int n, double r[], double s[] );
double r8mat_is_symmetric ( int m, int n, double a[] );
double *r8mat_jac ( int m, int n, double eps,
  double *(*fx) ( int m, int n, double x[] ), double x[] );
double *r8mat_kronecker ( int m1, int n1, double a[], int m2, int n2,
  double b[] );
double *r8mat_l_inverse ( int n, double a[] );
void r8mat_l_print ( int m, int n, double a[], string title );
double *r8mat_l_solve ( int n, double a[], double b[] );
double *r8mat_l1_inverse ( int n, double a[] );
double *r8mat_lt_solve ( int n, double a[], double b[] );
void r8mat_lu ( int m, int n, double a[], double l[], double p[], double u[] );
double r8mat_max ( int m, int n, double a[] );
double *r8mat_max_columns ( int m, int n, double a[] );
double *r8mat_max_rows ( int m, int n, double a[] );
void r8mat_max_index ( int m, int n, double a[], int &i_max, int &j_max );
double r8mat_maxcol_minrow ( int m, int n, double a[] );
double r8mat_maxrow_mincol ( int m, int n, double a[] );
double r8mat_mean ( int m, int n, double a[] );
double *r8mat_mean_columns ( int m, int n, double a[] );
double *r8mat_mean_rows ( int m, int n, double a[] );
double r8mat_min ( int m, int n, double a[] );
double *r8mat_min_columns ( int m, int n, double a[] );
double *r8mat_min_rows ( int m, int n, double a[] );
void r8mat_min_index ( int m, int n, double a[], int &i_min, int &j_min );
double r8mat_mincol_maxrow ( int m, int n, double a[] );
double r8mat_minrow_maxcol ( int m, int n, double a[] );
void r8mat_minvm ( int n1, int n2, double a[], double b[], double c[] );
double *r8mat_minvm_new ( int n1, int n2, double a[], double b[] );
void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] );
void r8mat_mtv ( int m, int n, double a[], double x[], double atx[] );
double *r8mat_mtv_new ( int m, int n, double a[], double x[] );
void r8mat_mv ( int m, int n, double a[], double x[], double ax[] );
double *r8mat_mv_new ( int m, int n, double a[], double x[] );
void r8mat_nint ( int m, int n, double a[] );
int r8mat_nonzeros ( int m, int n, double a[] );
double r8mat_norm_eis ( int m, int n, double a[] );
double r8mat_norm_fro ( int m, int n, double a[] );
double r8mat_norm_fro_affine ( int m, int n, double a1[], double a2[] );
double r8mat_norm_l1 ( int m, int n, double a[] );
double r8mat_norm_l2 ( int m, int n, double a[] );
double r8mat_norm_li ( int m, int n, double a[] );
double *r8mat_normal_01_new ( int m, int n, int &seed );
double *r8mat_nullspace ( int m, int n, double a[], int nullspace_size );
int r8mat_nullspace_size ( int m, int n, double a[] );
double *r8mat_orth_uniform_new ( int n, int &seed );
void r8mat_plot ( int m, int n, double a[], string title );
char r8mat_plot_symbol ( double r );
double *r8mat_poly_char ( int n, double a[] );
double *r8mat_power ( int n, double a[], int npow );
void r8mat_power_method ( int n, double a[], double *r, double v[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8mat_product_elementwise ( int m, int n, double a[], double b[] );
double r8mat_ref ( int m, int n, double a[] );
double r8mat_rms ( int m, int n, double a[] );
void r8mat_row_copy ( int m, int n, int i, double v[], double a[] );
double r8mat_rref ( int m, int n, double a[] );
void r8mat_scale ( int m, int n, double s, double a[] );
double *r8mat_scale_01 ( int m, int n, double x[] );
double *r8mat_scale_ab ( int m, int n, double x[], double a, double b );
int r8mat_solve ( int n, int nrhs, double a[] );
double *r8mat_solve_2d ( double a[], double b[], double *det );
double *r8mat_solve_3d ( double a[], double b[], double *det );
double *r8mat_solve2 ( int n, double a[], double b[], int &ierror );
double *r8mat_standardize ( int m, int n, double x[] );
double *r8mat_std_columns ( int m, int n, double a[] );
double *r8mat_std_rows ( int m, int n, double a[] );
double r8mat_sum ( int m, int n, double a[] );
double *r8mat_sum_columns ( int m, int n, double a[] );
double *r8mat_sum_rows ( int m, int n, double a[] );
double *r8mat_symm_eigen ( int n, double x[], double q[] );
void r8mat_symm_jacobi ( int n, double a[] );
double **r8mat_to_r8cmat_new (  int m, int n, double a[] );
int r8mat_to_r8plu ( int n, double a[], int pivot[], double lu[] );
double **r8mat_to_r8rmat ( int m, int n, double a[] );
double r8mat_trace ( int n, double a[] );
void r8mat_transpose_in_place ( int n, double a[] );
double *r8mat_transpose_new ( int m, int n, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
double *r8mat_u_inverse ( int n, double a[] );
double *r8mat_u_solve ( int n, double a[], double b[] );
double *r8mat_u1_inverse ( int n, double a[] );
void r8mat_uniform_01 ( int m, int n, int &seed, double r[] );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
double *r8mat_uniform_ab_new ( int m, int n, double a, double b, int &seed );
void r8mat_uniform_ab ( int m, int n, double a[], double b[], int &seed, double r[] );
double *r8mat_uniform_ab_new ( int m, int n, double a[], double b[], int &seed );
void r8mat_uniform_abvec ( int m, int n, double a[], double b[], int &seed, 
  double r[] );
double *r8mat_uniform_abvec_new ( int m, int n, double a[], double b[], int &seed );
double *r8row_uniform_new ( int m, int n, double a[], double b[], int &seed );
double *r8mat_ut_solve ( int n, double a[], double b[] );
double *r8mat_vand2 ( int n, double x[] );
double *r8mat_variance_columns ( int m, int n, double a[] );
double *r8mat_variance_rows ( int m, int n, double a[] );
double r8mat_vtmv ( int m, int n, double x[], double a[], double y[] );
void r8mat_zeros ( int m, int n, double a[] );
double *r8mat_zeros_new ( int m, int n );
double r8plu_det ( int n, int pivot[], double lu[] );
void r8plu_inverse ( int n, int pivot[], double lu[], double a_inverse[] );
void r8plu_mul ( int n, int pivot[], double lu[], double x[], double b[] );
void r8plu_sol ( int n, int pivot[], double lu[], double b[], double x[] );
void r8plu_to_r8mat ( int n, int pivot[], double lu[], double a[] );
void r8pp_delete ( int m, int n, double **a );
double **r8pp_new ( int m, int n );
int r8r8_compare ( double x1, double y1, double x2, double y2 );
void r8r8_print ( double a1, double a2, string title );
int r8r8r8_compare ( double x1, double y1, double z1, double x2, double y2,
  double z2 );
void r8r8r8vec_index_insert_unique ( int maxn, int &n, double x[], double y[],
  double z[], int indx[], double xval, double yval, double zval, int &ival,
  int &ierror );
void r8r8r8vec_index_search ( int n, double x[], double y[], double z[],
  int indx[], double xval, double yval, double zval, int &less, int &equal,
  int &more );
void r8r8vec_index_insert_unique ( int maxn, int &n, double x[], double y[],
  int indx[], double xval, double yval, int &ival, int &ierror );
void r8r8vec_index_search ( int n, double x[], double y[], int indx[],
  double xval, double yval, int &less, int &equal, int &more );
double **r8rmat_copy_new ( int m, int n, double **a );
void r8rmat_delete ( int m, int n, double **a );
double *r8rmat_fs_new ( int n, double **a, double b[] );
double **r8rmat_new ( int m, int n );
void r8rmat_print ( int m, int n, double **a, string title );
void r8rmat_print_some ( int m, int n, double **a, int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8rmat_to_r8mat ( int m, int n, double **a );
double **r8rmat_zeros ( int m, int n );
double *r8rows_to_r8mat ( int m, int n, double r8rows[] );
void r8slmat_print ( int m, int n, double a[], string title );
void r8vec_01_to_ab ( int n, double a[], double amax, double amin );
void r8vec_add ( int n, double a1[], double a2[] );
double r8vec_amax ( int n, double r8vec[] );
int r8vec_amax_index ( int n, double a[] );
double r8vec_amin ( int n, double r8vec[] );
int r8vec_amin_index ( int n, double a[] );
double *r8vec_any_normal ( int dim_num, double v1[] );
void r8vec_append ( int *n, double **a, double value );
double *r8vec_append_new ( int n, double a[], double value );
double r8vec_asum ( int n, double x[] );
void r8vec_bin ( int n, double x[], int bin_num, double bin_min, double bin_max,
  int bin[], double bin_limit[] );
void r8vec_binary_next ( int n, double x[] );
void r8vec_bracket ( int n, double x[], double xval, int &left,
  int &right );
void r8vec_bracket2 ( int n, double x[], double xval, int start, int &left,
  int &right );
void r8vec_bracket3 ( int n, double t[], double tval, int &left );
void r8vec_bracket4 ( int nt, double t[], int ns, double s[], int left[] );
int r8vec_bracket5 ( int nd, double xd[], double xi );
int *r8vec_bracket6 ( int nd, double xd[], int ni, double xi[] );
double *r8vec_cheby_extreme_new ( int n, double a, double b );
double *r8vec_cheby_zero_new ( int n, double a, double b );
double *r8vec_cheby1space_new ( int n, double a, double b );
double *r8vec_cheby2space_new ( int n, double a, double b );
int r8vec_compare ( int n, double a[], double b[] );
void r8vec_concatenate ( int n1, double a[], int n2, double b[], double c[] );
double *r8vec_concatenate_new ( int n1, double a[], int n2, double b[] );
double *r8vec_convolution ( int m, double x[], int n, double y[] );
double *r8vec_convolution_circ ( int n, double x[], double y[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double *r8vec_copy_new ( int n, double a1[] );
double r8vec_correlation ( int n, double x[], double y[] );
double r8vec_covariance ( int n, double x[], double y[] );
double r8vec_cross_product_2d ( double v1[2], double v2[2] );
double r8vec_cross_product_affine_2d ( double v0[2], double v1[2], double v2[2] );
double *r8vec_cross_product_3d ( double v1[3], double v2[3] );
double *r8vec_cross_product_affine_3d ( double v0[3], double v1[3], double v2[3] );
double *r8vec_cum_new ( int n, double a[] );
double *r8vec_cum0_new ( int n, double a[] );
double *r8vec_dif ( int n, double h );
double r8vec_diff_norm ( int n, double a[], double b[] );
double r8vec_diff_norm_l1 ( int n, double a[], double b[] );
double r8vec_diff_norm_l2 ( int n, double a[], double b[] );
double r8vec_diff_norm_li ( int n, double a[], double b[] );
double r8vec_diff_norm_squared ( int n, double a[], double b[] );
void r8vec_direct_product ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double x[] );
void r8vec_direct_product2 ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double w[] );
double r8vec_distance ( int dim_num, double v1[], double v2[] );
void r8vec_divide ( int n, double a[], double s );
double r8vec_dot_product ( int n, double v1[], double v2[] );
double r8vec_dot_product_affine ( int n, double v0[], double v1[], double v2[] );
double r8vec_entropy ( int n, double x[] );
bool r8vec_eq ( int n, double a1[], double a2[] );
void r8vec_even ( int n, double alo, double ahi, double a[] );
double *r8vec_even_new ( int n, double alo, double ahi );
double r8vec_even_select ( int n, double xlo, double xhi, int ival );
void r8vec_even2 ( int maxval, int nfill[], int nold, double xold[],
  int &nval, double xval[] );
double r8vec_even2_select ( int n, double xlo, double xhi, int ival );
void r8vec_even3 ( int nold, int nval, double xold[], double xval[] );
double *r8vec_expand_linear ( int n, double x[], int fat );
double *r8vec_expand_linear2 ( int n, double x[], int before, int fat, 
  int after );
void r8vec_fill ( int n, double value, double x[] );
double *r8vec_fill_new ( int n, double value );
int *r8vec_first_index ( int n, double a[], double tol );
double r8vec_frac ( int n, double a[], int k );
double *r8vec_fraction ( int n, double x[] );
bool r8vec_gt ( int n, double a1[], double a2[] );
void r8vec_heap_a ( int n, double a[] );
void r8vec_heap_d ( int n, double a[] );
int *r8vec_histogram ( int n, double a[], double a_lo, double a_hi, 
  int histo_num );
double *r8vec_house_column ( int n, double a[], int k );
double r8vec_i4vec_dot_product ( int n, double r8vec[], int i4vec[] );
double *r8vec_identity_row_new ( int n, int i );
void r8vec_index_delete_all ( int n, double x[], int indx[], double xval,
  int &n2, double x2[], int indx2[] );
void r8vec_index_delete_dupes ( int n, double x[], int indx[],
  int &n2, double x2[], int indx2[] );
void r8vec_index_delete_one ( int n, double x[], int indx[], double xval,
  int &n2, double x2[], int indx2[] );
void r8vec_index_insert ( int &n, double x[], int indx[], double xval );
void r8vec_index_insert_unique ( int &n, double x[], int indx[], double xval );
void r8vec_index_order ( int n, double x[], int indx[] );
void r8vec_index_search ( int n, double x[], int indx[], double xval, int &less,
  int &equal, int &more );
void r8vec_index_sort_unique ( int n, double x[], int &n2, double x2[],
  int indx2[] );
void r8vec_index_sorted_range ( int n, double r[], int indx[], double r_lo,
  double r_hi, int &i_lo, int &i_hi );
void r8vec_indexed_heap_d ( int n, double a[], int indx[] );
int r8vec_indexed_heap_d_extract ( int &n, double a[], int indx[] );
void r8vec_indexed_heap_d_insert ( int &n, double a[], int indx[],
  int indx_insert );
int r8vec_indexed_heap_d_max ( int n, double a[], int indx[] );
void r8vec_indicator0 ( int n, double a[] );
double *r8vec_indicator0_new ( int n );
void r8vec_indicator1 ( int n, double a[] );
double *r8vec_indicator1_new ( int n );
void r8vec_insert ( int n, double a[], int pos, double value );
bool r8vec_is_ascending ( int n, double x[] );
bool r8vec_is_ascending_strictly ( int n, double x[] );
bool r8vec_is_binary ( int n, double a[] );
bool r8vec_is_distinct ( int n, double x[] );
bool r8vec_is_in_01 ( int n, double x[] );
bool r8vec_is_in_ab ( int n, double x[], double a, double b );
bool r8vec_is_insignificant ( int n, double r[], double s[] );
bool r8vec_is_integer ( int n, double a[] );
bool r8vec_is_negative ( int n, double a[] );
bool r8vec_is_negative_any ( int n, double a[] );
bool r8vec_is_nonnegative ( int n, double x[] );
bool r8vec_is_nonpositive ( int n, double a[] );
bool r8vec_is_nonzero_any ( int n, double a[] );
bool r8vec_is_one ( int n, double x[] );
bool r8vec_is_positive ( int n, double a[] );
bool r8vec_is_zero ( int n, double x[] );
double *r8vec_legendre_new ( int n, double a_first, double a_last );
void r8vec_linspace ( int n, double a_lo, double a_hi, double x[] );
double *r8vec_linspace_new ( int n, double a_lo, double a_hi );
double *r8vec_linspace2_new ( int n, double a_lo, double a_hi );
bool r8vec_lt ( int n, double a1[], double a2[] );
void r8vec_mask_print ( int n, double a[], int mask_num, int mask[],
  string title );
double r8vec_max ( int n, double r8vec[] );
int r8vec_max_abs_index ( int n, double a[] );
int r8vec_max_index ( int n, double a[] );
double r8vec_mean ( int n, double x[] );
double r8vec_mean_geometric ( int n, double x[] );
double *r8vec_mean_running ( int n, double v[] );
double r8vec_median ( int n, double a[] );
void r8vec_mesh_2d ( int nx, int ny, double xvec[], double yvec[], 
  double xmat[], double ymat[] );
double *r8vec_midspace_new ( int n, double a_lo, double a_hi );
double r8vec_min ( int n, double r8vec[] );
int r8vec_min_index ( int n, double a[] );
double r8vec_min_pos ( int n, double a[] );
bool r8vec_mirror_next ( int n, double a[] );
void r8vec_mirror_ab_next ( int m, double a[], double b[], double x[], bool &done );
void r8vec_mm_to_01 ( int n, double a[] );
double *r8vec_mm_to_cd ( int n, double a[], double bmin, double bmax );
void r8vec_nint ( int n, double a[] );
double *r8vec_nint_new ( int n, double a[] );
double r8vec_norm ( int n, double v[] );
double r8vec_norm_affine ( int n, double v0[], double v1[] );
double r8vec_norm_l0 ( int n, double a[] );
double r8vec_norm_l1 ( int n, double a[] );
double r8vec_norm_l2 ( int n, double a[] );
double r8vec_norm_li ( int n, double a[] );
double r8vec_norm_lp ( int n, double a[], double p );
double r8vec_norm_rms ( int n, double a[] );
void r8vec_normal_01 ( int n, int &seed, double x[] );
double *r8vec_normal_01_new ( int n, int &seed );
double *r8vec_normal_ab_new ( int n, double b, double c, int &seed );
void r8vec_normalize ( int n, double a[] );
void r8vec_normalize_l1 ( int n, double a[] );
double r8vec_normsq ( int n, double a[] );
double r8vec_normsq_affine ( int n, double v0[], double v1[] );
double *r8vec_ones_new ( int n );
int r8vec_order_type ( int n, double x[] );
void r8vec_part_quick_a ( int n, double a[], int &l, int &r );
void r8vec_permute ( int n, int p[], double a[] );
void r8vec_permute_cyclic ( int n, int k, double a[] );
void r8vec_permute_uniform ( int n, double a[], int &seed );
void r8vec_polarize ( int n, double a[], double p[], double a_normal[],
  double a_parallel[] );
void r8vec_print ( int n, double a[], string title );
void r8vec_print_16 ( int n, double a[], string title );
void r8vec_print_part ( int n, double a[], int i_lo, int i_hi, string title );
void r8vec_print_some ( int n, double a[], int max_print, string title );
double r8vec_product ( int n, double a[] );
void r8vec_range ( int n, double x[], double xmin, double xmax, double y[],
  double *ymin, double *ymax );
void r8vec_range_2 ( int n, double a[], double *amin, double *amax );
void r8vec_reverse ( int n, double a[] );
void r8vec_rotate ( int n, double a[], int m );
double r8vec_scalar_triple_product ( double v1[3], double v2[3], double v3[3] );
int r8vec_search_binary_a ( int n, double a[], double aval );
void r8vec_scale ( double s, int n, double a[] );
double *r8vec_scale_01 ( int n, double x[] );
double *r8vec_scale_ab ( int n, double x[], double a, double b );
void r8vec_shift ( int shift, int n, double x[] );
void r8vec_shift_circular ( int shift, int n, double x[] );
double *r8vec_sign3_running ( int n, double v[] );
double *r8vec_softmax ( int n, double x[] );
void r8vec_sort_bubble_a ( int n, double a[] );
void r8vec_sort_bubble_d ( int n, double a[] );
void r8vec_sort_heap_a ( int n, double a[] );
void r8vec_sort_heap_d ( int n, double a[] );
void r8vec_sort_heap_index_a ( int n, double a[], int indx[] );
int *r8vec_sort_heap_index_a_new ( int n, double a[] );
void r8vec_sort_heap_index_d ( int n, double a[], int indx[] );
int *r8vec_sort_heap_index_d_new ( int n, double a[] );
int *r8vec_sort_heap_mask_a ( int n, double a[], int mask_num, int mask[] );
void r8vec_sort_insert_a ( int n, double a[] );
int *r8vec_sort_insert_index_a ( int n, double a[] );
void r8vec_sort_quick_a ( int n, double a[] );
void r8vec_sort_shell_a ( int n, double a[] );
double *r8vec_sorted_merge_a ( int na, double a[], int nb, double b[], int &nc );
int r8vec_sorted_nearest ( int n, double a[], double value );
void r8vec_sorted_range ( int n, double r[], double r_lo, double r_hi,
  int &i_lo, int &i_hi );
void r8vec_sorted_split ( int n, double a[], double split, int &i_lt, int &i_gt );
void r8vec_sorted_undex ( int x_num, double x_val[], int x_unique_num,
  double tol, int undx[], int xdnu[] );
double *r8vec_sorted_unique ( int n, double a[], double tol, int &unique_num );
int r8vec_sorted_unique_count ( int n, double a[], double tol );
void r8vec_sorted_unique_hist ( int n, double a[], double tol, int maxuniq,
  int &unique_num, double auniq[], int acount[] );
int r8vec_split ( int n, double a[], double split );
double *r8vec_standardize ( int n, double x[] );
double r8vec_std ( int n, double a[] );
double r8vec_std_sample ( int n, double a[] );
void r8vec_std_update ( int nm1, double mean_nm1, double std_nm1, double xn,
  int &n, double &mean_n, double &std_n );
void r8vec_step ( double x0, int n, double x[], double fx[] );
void r8vec_stutter ( int n, double a[], int m, double am[] );
double *r8vec_stutter_new ( int n, double a[], int m );
double r8vec_sum ( int n, double a[] );
double *r8vec_sum_running ( int n, double v[] );
void r8vec_swap ( int n, double a1[], double a2[] );
void r8vec_transpose_print ( int n, double a[], string title );
void r8vec_undex ( int x_num, double x_val[], int x_unique_num, double tol,
  int undx[], int xdnu[] );
void r8vec_uniform_01 ( int n, int &seed, double r[] );
double *r8vec_uniform_01_new ( int n, int &seed );
void r8vec_uniform_ab ( int n, double a, double b, int &seed, double x[] );
double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed );
void r8vec_uniform_abvec ( int n, double a[], double b[], int &seed, double x[] );
double *r8vec_uniform_abvec_new ( int n, double a[], double b[], int &seed );
double *r8vec_uniform_unit_new ( int dim_num, int &seed );
int r8vec_unique_count ( int n, double a[], double tol );
int *r8vec_unique_index ( int n, double a[], double tol );
double r8vec_variance ( int n, double x[] );
double r8vec_variance_circular ( int n, double x[] );
double r8vec_variance_sample ( int n, double x[] );
void r8vec_variance_update ( int nm1, double mean_nm1, double variance_nm1, 
  double xn, int &n, double &mean_n, double &variance_n );
double *r8vec_vector_triple_product ( double v1[3], double v2[3], double v3[3] );
void r8vec_write ( int n, double r[], string output_file );
void r8vec_zeros ( int n, double a1[] );
double *r8vec_zeros_new ( int n );
int r8vec2_compare ( int n, double a1[], double a2[], int i, int j );
void r8vec2_print ( int n, double a1[], double a2[], string title );
void r8vec2_print_some ( int n, double x1[], double x2[], int max_print,
  string title );
void r8vec2_sort_a ( int n, double a1[], double a2[] );
void r8vec2_sort_d ( int n, double a1[], double a2[] );
int *r8vec2_sort_heap_index_a ( int n, double x[], double y[] );
void r8vec2_sorted_unique ( int n, double a1[], double a2[], int &unique_num );
void r8vec2_sorted_unique_index ( int n, double a1[], double a2[],
  int &unique_num, int indx[] );
int r8vec2_sum_max_index ( int n, double a[], double b[] );
void r8vec3_print ( int n, double a1[], double a2[], double a3[], string title );
int s_len_trim ( string s );
void sort_heap_external ( int n, int &indx, int &i, int &j, int isgn );
void timestamp ( );
