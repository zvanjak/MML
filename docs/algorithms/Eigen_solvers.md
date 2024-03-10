# Eigen solvers

MML implements two classes for solving eigenvalue problems: `SymmMatEigenSolver` for real symmetric matrices and `EigenSolver` for real unsymmetric matrices. 

## Real symmetric matrix eigen solver

Class `SymmMatEigenSolver` computes all eigenvalues and eigenvectors of a real symmetric matrix by reduction to tridiagonal form followed by QL iteration.
Calculation happens in constructor, which takes the matrix as input. If the second argument is set to `true`, the eigenvectors are also computed.

## Example usage
~~~C++
MatrixSym<Real>     symm_mat{ 5,  {1.4, 2.1, 2.1, 7.4, 9.6,
                                        1.5, 1.1, 0.7, 5.0,
                                             9.6, 5.4, 8.8,
                                                  0.4, 8.0,
                                                       7.7 } };
Matrix<Real>        matcopy = symm_mat.GetAsMatrix();

SymmMatEigenSolver eigen_solver(symm_mat, true);

int width = 15;
int prec = 10;

std::cout << "Initial matrix:\n";  matcopy.Print(std::cout, 10, 3);  std::cout << std::endl;
std::cout << "Eigenvalues         : "; eigen_solver.getEigenvalues().Print(std::cout, 10, 3); std::cout << "\n\n";

for (int i = 0; i < symm_mat.Dim(); i++)
{
		double        eigen_val = eigen_solver.getEigenvalues()[i];
		Vector<Real>  eigen_vec = eigen_solver.getEigenvector(i);

		std::cout << "Eigen value  " << i + 1 << " (real)    = " << eigen_val << std::endl;
		std::cout << "Eigen vector " << i + 1 << " (real)    = "; eigen_vec.Print(std::cout, width, prec); std::cout << std::endl;

		std::cout << "Mat * eig_vec            = "; (matcopy * eigen_vec).Print(std::cout, width, prec); std::cout << std::endl;
		std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";
}

/* OUTPUT
Initial matrix:
Rows: 5 Cols: 5
[        1.4,        2.1,        7.4,        1.1,        5.4,  ]
[        2.1,        2.1,        9.6,        0.7,        8.8,  ]
[        7.4,        9.6,        1.5,          5,        0.4,  ]
[        1.1,        0.7,          5,        9.6,          8,  ]
[        5.4,        8.8,        0.4,          8,        7.7,  ]

Eigenvalues         : [      24.6,       8.01,       3.54,      -0.57,      -13.3]

Eigen value  1 (real)    = 24.6
Eigen vector 1 (real)    = [   0.3207830408,    0.4359308891,    0.3974911379,    0.4784849759,    0.5657874369]
Mat * eig_vec            = [    7.887571178,     10.71888311,     9.773707598,     11.76522392,     13.91185977]
eig_val * eig_vec        = [    7.887571178,     10.71888311,     9.773707598,     11.76522392,     13.91185977]

Eigen value  2 (real)    = 8.008425134
Eigen vector 2 (real)    = [    0.323313998,    0.4160284348,    0.4267540611,   -0.7051973277,   -0.2072826582]
Mat * eig_vec            = [    2.589235948,     3.331732574,     3.417627949,    -5.647520004,     -1.66000765]
eig_val * eig_vec        = [    2.589235948,     3.331732574,     3.417627949,    -5.647520004,     -1.66000765]

Eigen value  3 (real)    = 3.541665812
Eigen vector 3 (real)    = [   0.1398659018,   -0.2065443391,    0.5233038154,    0.4572026692,   -0.6744596367]
Mat * eig_vec            = [   0.4953582827,   -0.7315110245,     1.853367233,     1.619259063,    -2.388710637]
eig_val * eig_vec        = [   0.4953582827,   -0.7315110245,     1.853367233,     1.619259063,    -2.388710637]

Eigen value  4 (real)    = -0.5700940778
Eigen vector 4 (real)    = [   0.8045643778,   -0.5452049336,   -0.1559364946,  -0.09146643767,    0.1508159234]
Mat * eig_vec            = [   -0.458677387,    0.3108181039,    0.0888984721,   0.05214447444,  -0.08597926477]
eig_val * eig_vec        = [   -0.458677387,    0.3108181039,    0.0888984721,   0.05214447444,  -0.08597926477]

Eigen value  5 (real)    = -13.2684887
Eigen vector 5 (real)    = [   0.3545137387,    0.5449546933,   -0.6014305388,    0.2373783201,   -0.3990955168]
Mat * eig_vec            = [   -4.703861536,    -7.230725189,     7.980074307,    -3.149651558,     5.295394355]
eig_val * eig_vec        = [   -4.703861536,    -7.230725189,     7.980074307,    -3.149651558,     5.295394355]*
*/
~~~

## Real unsymmetric matrix eigen solver

Class `EigenSolver` computes all eigenvalues and eigenvectors of a real nonsymmetric matrix by reduction to Hessenberg form followed by QR iteration.

## Example usage
~~~C++
Matrix<Real>     origMat{ 5, 5, { 3.2,-4.1, 2.7, 3.4, 4.6,
                                  2.1, 3.5, 1.6,-0.7, 5.0,
                                  3.8,-1.3,-6.6,-5.4, 3.8,
                                  4.6, 8.2,-8.4, 0.4, 8.0,
                                  2.6, 2.9, 0.1, 9.6,-7.7 } };
Matrix<Real>  matcopy(origMat);

EigenSolver eigen_solver(matcopy, true, false);

int width = 15;
int prec = 10;
std::cout << "Initial matrix:\n";  matcopy.Print(std::cout, 10, 3);  std::cout << std::endl;
	
std::cout << "Num real eigenvalues    : " << eigen_solver.getNumReal() << std::endl;
std::cout << "Num complex eigenvalues : " << eigen_solver.getNumComplex() << "\n\n";
	
std::cout << "Real eigenvalues    : "; eigen_solver.getRealEigenvalues().Print(std::cout, width, prec); std::cout << std::endl;
std::cout << "Complex eigenvalues : "; eigen_solver.getComplexEigenvalues().Print(std::cout, width, prec); std::cout << "\n\n";

Matrix<Complex> mat_cmplx = MatrixUtils::CmplxMatFromRealMat(matcopy);

for (int i = 0; i < origMat.RowNum(); i++)
{
	if (eigen_solver.isRealEigenvalue(i))
	{
		double        eigen_val = eigen_solver.getEigenvalues()[i].real();
		Vector<Real>  eigen_vec = eigen_solver.getRealPartEigenvector(i);

		std::cout << "Eigen value  " << i + 1 << " (real)    = " << eigen_val << std::endl;
		std::cout << "Eigen vector " << i + 1 << " (real)    = "; eigen_vec.Print(std::cout, width, prec); std::cout << std::endl;

		std::cout << "Mat * eig_vec            = "; (matcopy * eigen_vec).Print(std::cout, width, prec); std::cout << std::endl;
		std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";
	}
	else
	{
		Complex          eigen_val = eigen_solver.getEigenvalues()[i];
		Vector<Complex>  eigen_vec = eigen_solver.getEigenvector(i);

		std::cout << "Eigen value  " << i + 1 << " (complex) = " << eigen_val << std::endl;
		std::cout << "Eigen vector " << i + 1 << " (complex) = "; eigen_solver.getEigenvector(i).Print(std::cout, width, prec); std::cout << std::endl;

		std::cout << "Mat * eig_vec            = "; (mat_cmplx * eigen_vec).Print(std::cout, width, prec); std::cout << std::endl;
		std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";
	}
}

/* OUTPUT
Initial matrix:
Rows: 5 Cols: 5
[        3.2,       -4.1,        2.7,        3.4,        4.6,  ]
[        2.1,        3.5,        1.6,       -0.7,          5,  ]
[        3.8,       -1.3,       -6.6,       -5.4,        3.8,  ]
[        4.6,        8.2,       -8.4,        0.4,          8,  ]
[        2.6,        2.9,        0.1,        9.6,       -7.7,  ]

Num real eigenvalues    : 3
Num complex eigenvalues : 2

Real eigenvalues    : [    11.76342368,    -7.588677062,    -16.23751082]
Complex eigenvalues : [(2.431382106,3.808290416), (2.431382106,-3.808290416)]

Eigen value  1 (real)    = 11.76342368
Eigen vector 1 (real)    = [   0.3339842192,    0.2643676911,  -0.05813239922,    0.6554857426,    0.4070134471]
Mat * eig_vec            = [    3.928797872,     3.109869156,   -0.6838360413,     7.710756503,      4.78787162]
eig_val * eig_vec        = [    3.928797872,     3.109869156,   -0.6838360413,     7.710756503,      4.78787162]

Eigen value  2 (complex) = (2.431382106,3.808290416)
Eigen vector 2 (complex) = [(-0.2890819407,-0.9295540505), (-0.5477470685,0.2573057905), (-0.4005918445,-0.4458267692), (0.44828404,0.3937093584), (0.2334509995,0.1160082809)]
Mat * eig_vec            = [(2.837143124,-3.361009069), (-2.311677597,-1.460371216), (0.7238459701,-2.609545311), (-0.4094097833,2.664453702), (0.1258153585,1.171109662)]
eig_val * eig_vec        = [(2.837143124,-3.361009069), (-2.311677597,-1.460371216), (0.7238459701,-2.609545311), (-0.4094097833,2.664453702), (0.1258153585,1.171109662)]

Eigen value  3 (complex) = (2.431382106,-3.808290416)
Eigen vector 3 (complex) = [(-0.2890819407,0.9295540505), (-0.5477470685,-0.2573057905), (-0.4005918445,0.4458267692), (0.44828404,-0.3937093584), (0.2334509995,-0.1160082809)]
Mat * eig_vec            = [(2.837143124,3.361009069), (-2.311677597,1.460371216), (0.7238459701,2.609545311), (-0.4094097833,-2.664453702), (0.1258153585,-1.171109662)]
eig_val * eig_vec        = [(2.837143124,3.361009069), (-2.311677597,1.460371216), (0.7238459701,2.609545311), (-0.4094097833,-2.664453702), (0.1258153585,-1.171109662)]

Eigen value  4 (real)    = -7.588677062
Eigen vector 4 (real)    = [   0.7226386328,    0.3241902749,   -0.4495106792,   -0.2996399104,    -0.920582649]
Mat * eig_vec            = [   -5.483871217,    -2.460175303,      3.41119138,     2.273870515,     6.986004432]
eig_val * eig_vec        = [   -5.483871217,    -2.460175303,      3.41119138,     2.273870515,     6.986004432]

Eigen value  5 (real)    = -16.23751082
Eigen vector 5 (real)    = [ -0.01376371301,   -0.1611250262,   -0.6690299644,   -0.6241081545,    0.7685365962]
Mat * eig_vec            = [    0.223488439,     2.616269357,     10.86338129,     10.13396291,     -12.4791213]
eig_val * eig_vec        = [    0.223488439,     2.616269357,     10.86338129,     10.13396291,     -12.4791213]*
*/
~~~

