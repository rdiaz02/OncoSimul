/*
	All functions needed for Linear Algrebra
	This is especially used for Least Square Estimates
*/


#ifndef __LINEAR_ALGEBRA__
#define __LINEAR_ALGEBRA__

struct matrix {
	int r; 
	int c;
	double **val;
};

/*
	Memory functions
*/

struct matrix MemMat( int m, int n );
void IdentMat(struct matrix *M);
void FreeMat( struct matrix *M );

/*
	Output a matrix
*/
void PrintMat( struct matrix *M, char *name );

/*
	Basic matrix operators
	when opt_free is set to 1, it also free the entry matrix
*/

// :: WHAT ABOUT MEMORY HERE?:: WOULD IT BE BETTER TO MAKE ALL THESE FUNCTIONS VOID AND PASS THE MATRIX AS REFERENCE?:: //

struct matrix MatrixProduct(struct matrix *M_left, struct matrix *M_right, char opt_free);
void MatrixProductRef(struct matrix *M_left, struct matrix *M_right, struct matrix *Prod, char opt_free);

struct matrix MatrixSum(struct matrix *M1, struct matrix *M2, char opt_free);
void MatrixSumRef(struct matrix *M1, struct matrix *M2, struct matrix *Sum, char opt_free);

struct matrix MatrixDiff(struct matrix *M1, struct matrix *M2, char opt_free);
void MatrixDiffRef(struct matrix *M1, struct matrix *M2, struct matrix *Diff, char opt_free);


struct matrix MatrixTranspose( struct matrix *M, char opt_free );
struct matrix Invert_MatrixDiag(struct matrix *D, char opt_free );
struct matrix MatrixHadarmardProduct(struct matrix *M1, struct matrix *M2, char opt_free);
void MatrixHadarmardProductRef(struct matrix *M1, struct matrix *M2, struct matrix *HadarmardProduct, char opt_free);
struct matrix scaleMatrix(struct matrix *M1, float scale, char opt_free);
void scaleMatrixRef(struct matrix *M1, 	struct matrix *M2, float scale, char opt_free);
void inverseMatrix(struct matrix *M1, struct matrix *inverseM);

/*
	From a vector, create a Diag Matrix
*/
struct matrix CreateMatrixDiag( struct matrix *Vector, char opt_free );

/*
	From a Matrix, create a Diag Matrix
 */
struct matrix CreateMatrixDiagMatrix( struct matrix *M1, char opt_free );

/*
	Return a vector [n,1] from a matrix X [m,n] and a vector Y[m,1]
	m i the number of observation
	n is the number of parameter to estimate
*/
struct matrix LeastSquare( struct matrix *X, struct matrix *Y );

#endif
