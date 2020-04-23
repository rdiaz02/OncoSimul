/*
	All functions needed for Linear Algrebra
	This is especially used for Least Square Estimates
*/

/*
	All functions needed for Linear Algrebra
	This is especially used for Least Square Estimates
	
	
	date: Nov 2103
	author: gachaz
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LinearAlgebra.h"


/*
	Memory functions
*/

struct matrix MemMat( int m, int n ){

	int i;
	struct matrix M;
	
	M.r=m;
	M.c=n;
	
	M.val = (double **) malloc( (size_t) M.r * sizeof(double *) );
	if( ! M.val ) fprintf(stderr, "MemMat: cannot allocate matrix, bye\n"), exit(1);
	
	for(i=0; i< M.r ; i++){
		M.val[i] = (double *) malloc( (size_t) M.c * sizeof(double) );
		if(!M.val[i]) fprintf(stderr, "MemMat: cannot allocate M[%d], bye\n",i), exit(1);
	}

	return M;
}

void IdentMat(struct matrix *M)
{
	
	int i, j;
	
	for (i = 0; i < M->r; i++)
	{
		for (j = 0; j < M->c; j++)
		{
			if (i==j)
			{
				M->val[i][i] = 1.;
			}
			else
			{
				M->val[i][j] = 0.;
			}
		}
	}
}


void FreeMat( struct matrix *M ){

	int i;
	
	if(M->val == NULL)
		return;
	
	for(i=0;i<M->r; i++)
		free( M->val[i] );

	free(M->val);
	
	M->val =NULL;
}


/*
	Output a matrix
*/
void PrintMat( struct matrix *M, char *name ){

	int i,j;

	printf("%s [%d,%d]\n", name, M->r, M->c  );

	for(i=0;i<M->r; i++, printf("\n") )
		for(j=0;j<M->c;j++)
			printf("%8.3f", M->val[i][j]);

	printf("//\n");
}



struct matrix MatrixTranspose( struct matrix *M, char opt_free ){

	int i,j;
	struct matrix T = MemMat( M->c, M->r );	
	
	for(i=0 ; i<M->r ; i++)
		for(j=0 ; j<M->c ; j++)
			T.val[j][i] = M->val[i][j];
	
	
	if(opt_free)
		FreeMat(M);
	
	return T;
	
}

/*
	From a vector, create a Diag Matrix
*/
struct matrix CreateMatrixDiag( struct matrix *V, char opt_free )
{
	int i,j;
	
	struct matrix D = MemMat( V->r,V->r );
	
	for(i=0; i<D.r; i++ )
	{
		for(j=0; j<D.c; j++ )
		{
			D.val[i][j]=0;
			
			if( i == j )
			{
				D.val[i][j] = V->val[i][0];
			}
		}
		
	}
	
	if(opt_free)
	{
		FreeMat(V);
	}
	
	return D;
}


struct matrix CreateMatrixDiagMatrix( struct matrix *M1, char opt_free )
{
	int i,j;
	
	struct matrix D = MemMat( M1->r,M1->c );
	
	for(i=0; i<D.r; i++ )
	{
		for(j=0; j<D.c; j++ )
		{
			D.val[i][j] = 0;
			if( i == j )
			{
				D.val[i][j] = M1->val[i][i];
			}
			
		}
	}
	
	if(opt_free)
	{
		FreeMat(M1);
	}
	
		return (D);
}



struct matrix Invert_MatrixDiag(struct matrix *D, char opt_free )
{
	int i,j;
	struct matrix M;
	
	M = MemMat( D->r, D->c );
		
	for(i=0; i<M.r; i++ )
		for(j=0; j<M.c; j++ )
		{
			M.val[i][j]=0;
			
			if( i == j )
				M.val[i][j] = 1.0/D->val[i][j];
			
		}
	
	if(opt_free)
		FreeMat(D);
	

	return M;
		
}


struct matrix MatrixProduct(struct matrix *M_left, struct matrix *M_right, char opt_free)
{

	int i,j,k;
	struct matrix P;
	
	if( M_left->c != M_right->r )
		fprintf(stderr, "MatrixProduct: The matrices cannot be multiplied together. Incompatible dimensions, bye\n"),exit(1);
	
	
	P=MemMat( M_left->r, M_right->c );
	
	for(i=0; i<P.r; i++ )
	{
		for(j=0; j<P.c; j++ )
		{
			P.val[i][j]=0;
			
			for(k=0; k<M_left->c; k++ )
			{
				P.val[i][j] += M_left->val[i][k] * M_right->val[k][j];
			}
		}
	}
	
	if(opt_free)
	{
		FreeMat(M_right);
		FreeMat(M_left);
	}
	
	return P;
			
}


void MatrixProductRef(struct matrix *M_left, struct matrix *M_right, struct matrix * Prod, char opt_free)
{
	
  int i,j,k;
	
  if( M_left->c != M_right->r )
    fprintf(stderr, "MatrixProduct: The matrices cannot be multiplied together. Incompatible dimensions, bye\n"),exit(1);
	
		
  for(i=0; i<Prod->r; i++ )
    {
      for(j=0; j<Prod->c; j++ )
	{
	  Prod->val[i][j]=0;
		
	  for(k=0; k<M_left->c; k++ )
	    {
	      Prod->val[i][j] += M_left->val[i][k] * M_right->val[k][j];
	    }
	}
    }
	
  if(opt_free)
    {
      FreeMat(M_right);
      FreeMat(M_left);
    }
	
	
}

struct matrix MatrixSum(struct matrix *M1, struct matrix *M2, char opt_free)
{

	int i,j;
	struct matrix S;
	
	
	if( M1->r != M2->r || M1->c != M2->c )
		fprintf(stderr, "MatrixSum: The matrices cannot be summed together. Incompatible dimensions, bye\n"),exit(1);
	

	S = MemMat( M1->r, M1->c );
		
	for(i=0; i<S.r; i++ )
	{
		for(j=0; j<S.c; j++ )
		{
			S.val[i][j] = M1->val[i][j] + M2->val[i][j];
		}
	}
	
	if(opt_free)
	{
		FreeMat(M1);
		FreeMat(M2);
	}
	
	return S;
}

void MatrixSumRef(struct matrix *M1, struct matrix *M2, struct matrix *Sum, char opt_free)
{
	
  int i,j;
	
  if( M1->r != M2->r || M1->c != M2->c )
    fprintf(stderr, "MatrixSum: The matrices cannot be summed together. Incompatible dimensions, bye\n"),exit(1);
		
	
		
  for(i=0; i<Sum->r; i++ )
    {
      for(j=0; j<Sum->c; j++ )
	{
	  Sum->val[i][j] = M1->val[i][j] + M2->val[i][j];
	}
    }
				
  if(opt_free)
    {
      FreeMat(M1);
      FreeMat(M2);
    }
}




struct matrix MatrixDiff(struct matrix *M1, struct matrix *M2, char opt_free)
{
	
  int i,j;
  struct matrix S;
	
	
  if( M1->r != M2->r || M1->c != M2->c )
    fprintf(stderr, "MatrixSum: The matrices cannot be substracted. Incompatible dimensions, bye\n"),exit(1);
		
		
  S = MemMat( M1->r, M1->c );
		
  for(i=0; i<S.r; i++ )
    {
      for(j=0; j<S.c; j++ )
	{
	  S.val[i][j] = M1->val[i][j] - M2->val[i][j];
	}
    }
				
  if(opt_free)
    {
      FreeMat(M1);
      FreeMat(M2);
    }
	
  return S;
}

void MatrixDiffRef(struct matrix *M1, struct matrix *M2, struct matrix * Diff, char opt_free)
{
	
	int i,j;

	if( M1->r != M2->r || M1->c != M2->c )
		fprintf(stderr, "MatrixSum: The matrices cannot be substracted. Incompatible dimensions, bye\n"),exit(1);
	
	
	for(i=0; i<Diff->r; i++ )
	{
		for(j=0; j<Diff->c; j++ )
		{
			Diff->val[i][j] = M1->val[i][j] - M2->val[i][j];
		}
	}
	
	
	if(opt_free)
	{
		FreeMat(M1);
		FreeMat(M2);
	}
	
}


struct matrix MatrixHadarmardProduct(struct matrix *M1, struct matrix *M2, char opt_free)
{
	struct matrix HadarmardProduct;
	int i, j;

	if( M1->r != M2->r || M1->c != M2->c )
	{
		fprintf(stderr, "MatrixHadarmardProduct: The matrices cannot be multiplied together. Incompatible dimensions, bye\n"),exit(1);
	}

	HadarmardProduct = MemMat(M1->r, M1->c);
	
	for(i = 0; i < M1->r; i++)
	{
		for(j = 0; j < M1->c; j++)
		{
			HadarmardProduct.val[i][j] = M1->val[i][j] * M2->val[i][j];
		}
	}
	
	if(opt_free)
	{
		FreeMat(M1);
		FreeMat(M2);
	}
	
	return (HadarmardProduct);
}


void MatrixHadarmardProductRef(struct matrix *M1, struct matrix *M2, struct matrix * HadarmardProduct, char opt_free)
{
	int i, j;
	
	if( M1->r != M2->r || M1->c != M2->c )
	{
		fprintf(stderr, "MatrixHadarmardProduct: The matrices cannot be multiplied together. Incompatible dimensions, bye\n"),exit(1);
	}
	
	for(i = 0; i < M1->r; i++)
	{
		for(j = 0; j < M1->c; j++)
		{
			HadarmardProduct->val[i][j] = M1->val[i][j] * M2->val[i][j];
		}
	}
	
	if(opt_free)
	{
		FreeMat(M1);
		FreeMat(M2);
	}
	
}


struct matrix scaleMatrix(struct matrix *M1, float scale, char opt_free)
{
	struct matrix M2;
	int i, j;
	
	M2 = MemMat(M1->r, M1->c);
	for(i = 0; i < M1->r; i++)
	{
		for(j = 0; j < M1->c; j++)
		{
			M2.val[i][j] = scale * M1->val[i][j];
		}
	}
	
	if(opt_free)
	{
		FreeMat(M1);
	}
	
	return (M2);
}

void scaleMatrixRef(struct matrix *M1, 	struct matrix * M2, float scale, char opt_free)
{
	int i, j;
	
	for(i = 0; i < M1->r; i++)
	{
		for(j = 0; j < M1->c; j++)
		{
			M2->val[i][j] = scale * M1->val[i][j];
		}
	}
	
	if(opt_free)
	{
		FreeMat(M1);
	}
	
}

/*
	Do the SVD decomposition of any A matrix --dimension (m,n)--
	A = U.w.Vt
	U is orthogonal and dimension (m,n)
	w is a diagonal matrix of size (n,n)
	V is a square orthogonal matrix (n,n)
	
	The code has been adapted from numrec 92 -- reset mat from 0 to n-1 as in std C --
	
 */

#define CONVERGENCE 100

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))


static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))


/*
	Apparently clever way to do (a^2 + b^2)^1/2
 */
static float pythag(float a, float b)
{
	float absa=fabs(a),
	absb=fabs(b);
	
	if (absa > absb)
		return absa*sqrt(1.0+DSQR(absb/absa));
	else
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

static void svdcmp(double **A, int m, int n, double **u, double *w, double **v)
{
	int flag,i,its,j,jj,k,l,nm=0;
	
	double anorm,c,f,g,h,s,scale,x,y,z,
	*rv1;
	
	if(u==NULL || w==NULL || v==NULL)
		fprintf(stderr, "you should get the memory before using this function, bye\n"),exit(3);
	
	
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			u[i][j]= A[i][j];
	
	rv1=(double *)malloc( (size_t) n*sizeof(double) );
	if(!rv1)fprintf(stderr, "cannot allocate rv1, bye");
	
	g=scale=anorm=0.0;
	
	for (i=0;i<n;i++) {
		
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		
		if (i < m)
		{
			
			for (k=i;k<m;k++)
				scale += fabs(u[k][i]);
			
			if (scale)
			{
				
				for (k=i;k<m;k++)
				{
					u[k][i] /= scale;
					s += u[k][i]*u[k][i];
				}
				
				f=u[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i][i]=f-g;
				
				for (j=l;j<n;j++)
				{
					for (s=0.0,k=i;k<m;k++)
						s += u[k][i]*u[k][j];
					f=s/h;
					
					for (k=i;k<m;k++)
						u[k][j] += f*u[k][i];
				}
				
				for (k=i;k<m;k++)
					u[k][i] *= scale;
			}
		}
		
		
		w[i]=scale *g;
		g=s=scale=0.0;
		
		if (i < m && i != n-1)
		{
			for (k=l;k<n;k++)
				scale += fabs(u[i][k]);
			
			if (scale)
			{
				
				for (k=l;k<n;k++)
				{
					u[i][k] /= scale;
					s += u[i][k]*u[i][k];
				}
				
				f=u[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i][l]=f-g;
				
				for (k=l;k<n;k++)
					rv1[k]=u[i][k]/h;
				
				for (j=l;j<m;j++)
				{
					for (s=0.0,k=l;k<n;k++)
						s += u[j][k]*u[i][k];
					
					for (k=l;k<n;k++)
						u[j][k] += s*rv1[k];
				}
				
				for (k=l;k<n;k++)
					u[i][k] *= scale;
			}
		}
		
		anorm = DMAX( anorm , (fabs(w[i])+fabs(rv1[i])) );
		
	}
	
	
	
	for (i=n-1;i>=0;i--)
	{
		if (i < n-1){
			
			if (g) {
				for (j=l;j<n;j++)
					v[j][i]=(u[i][j]/u[i][l])/g;
				
				for (j=l;j<n;j++)
				{
					for (s=0.0,k=l;k<n;k++)
						s += u[i][k]*v[k][j];
					
					for (k=l;k<n;k++)
						v[k][j] += s*v[k][i];
				}
			}
			
			for (j=l;j<n;j++)
				v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	
	
	for (i=IMIN(m-1,n-1);i>=0;i--) {
		
		l=i+1;
		g=w[i];
		
		for (j=l;j<n;j++)
			u[i][j]=0.0;
		
		if (g)
		{
			g=1.0/g;
			
			for (j=l;j<n;j++)
			{
				for (s=0.0,k=l;k<m;k++)
					s += u[k][i]*u[k][j];
				
				f=(s/u[i][i])*g;
				
				for (k=i;k<m;k++)
					u[k][j] += f*u[k][i];
			}
			
			for (j=i;j<m;j++)
				u[j][i] *= g;
		}
		else
			for (j=i;j<m;j++)
				u[j][i]=0.0;
		++u[i][i];
	}
	
	
	
	for (k=n-1;k>=0;k--) {
		
		for (its=1;its<=CONVERGENCE;its++) {
			
			flag=1;
			
			for (l=k;l>0;l--) {
				
				nm=l-1;
				if ( (fabs(rv1[l])+anorm) == anorm ){
					flag=0;
					break;
				}
				
				if ( (fabs(w[nm])+anorm) == anorm )
					break;
			}
			
			if (flag) {
				
				c=0.0;
				s=1.0;
				
				for (i=l;i<k;i++) {
					
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ( (fabs(f)+anorm) == anorm ) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					
					for (j=0;j<m;j++) {
						
						y=u[j][nm];
						z=u[j][i];
						u[j][nm]=y*c+z*s;
						u[j][i]=z*c-y*s;
					}
				}
			}
			
			z=w[k];
			
			if (l == k) {
				
				if (z < 0.0)
				{
					w[k] = -z;
					for (j=0;j<n;j++)
						v[j][k] = -v[j][k];
				}
				break;
			}
			
			if (its == CONVERGENCE)
				fprintf(stderr, "no convergence in %d svdcmp iterations, bye\n", CONVERGENCE), exit(5);
			
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			
			for (j=l;j<=nm;j++) {
				
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				
				for (jj=0;jj<n;jj++) {
					
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				
				z=pythag(f,h);
				w[j]=z;
				
				if (z)
				{
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				
				f=c*g+s*y;
				x=c*y-s*g;
				
				for (jj=0;jj<m;jj++)
				{
					y=u[jj][j];
					z=u[jj][i];
					u[jj][j]=y*c+z*s;
					u[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	
	free(rv1);
}



/*
	Calculating the inverse of a matrix M1 by SVD
*/

// :: FOR LARGE MATRICES THIS BECOMES NUMERICALLY UNSTABLE :: ANY BETTER OPTION AVAILABLE? :: ALTERNATIVELY ONE COULD JUST DO BRUTE FORCE MATRIX MULTIPLICATION :: WOULD BE OF COMPLEXITY O(x n^2) FOR A MATRIX OF DIMENSION n THAT IS RAISED TO THE POWER x ::

void inverseMatrix(struct matrix * M1, struct matrix * inverseM)
{
	
	struct matrix U,w,V;    /* misc , SVD decomposition (U,w,V) and matrix product (Tmp) */
	
	int i,j,k;                  /* counters */
	
	/*
		get memory for the matrices
	 */
	U = MemMat( M1->r, M1->c);
	w = MemMat( 1, M1->c );
	V = MemMat( M1->c, M1->c );

	
	/*
	 svd decomposition --adapted from numrec--
	 svdcmp(double **A, int m, int n, double **u, double *w, double **v)
	 */
	
	svdcmp(M1->val, M1->r, M1->c, U.val, w.val[0], V.val);

	/*
		Check the w's for 0 values
	 */
	for(k = 0; k < inverseM->r; k++)
	{
		if( w.val[0][k] < 1e-6 )
		{
			fprintf(stderr, "warning w[%d]=%e, need extra check\n", k, w.val[0][k]);
		}
	}
			
			
	/*
	 This is equivalent to
	 V . Diag(w)^-1 . Ut
	 but more efficient in computation time O(m n^2)
	 */
	for(i = 0; i < inverseM->r; i++ )
	{
		for(j = 0; j < inverseM->c; j++ )
		{
			inverseM->val[i][j]=0;
			
			for(k = 0; k < inverseM->r; k++ )
			{
				inverseM->val[i][j] += V.val[i][k]*U.val[j][k]/w.val[0][k];
			}
		}
	}
	
	
	FreeMat(&U);
	FreeMat(&V);
	FreeMat(&w);
}



struct matrix LeastSquare( struct matrix *X, struct matrix *Y ){


	struct matrix beta;          /* this will contains the estimates you aim at */
	
	struct matrix Tmp, U,w,V;    /* misc , SVD decomposition (U,w,V) and matrix product (Tmp) */
	
	int i,j,k;                  /* counters */


	/*
		get memory for the matrices
	*/
	beta = MemMat( X->c, 1);
	Tmp = MemMat( X->c, X->r);
	U = MemMat( X->r, X->c);
	w = MemMat( 1, X->c );
	V = MemMat( X->c, X->c );
	
	
	/*
	 	svd decomposition --adapted from numrec--
	 	svdcmp(double **A, int m, int n, double **u, double *w, double **v)
	*/

	svdcmp(X->val,  X->r, X->c, U.val, w.val[0], V.val);
	
	
	/*
		Check the w's for 0 values
	*/
	for(k=0; k<Tmp.r; k++)
		if( w.val[0][k] < 1e-6 )
			fprintf(stderr, "warning w[%d]=%e, need extra check\n", k, w.val[0][k]);
	
	
	
	/*
		This is equivalent to 
		V . Diag(w)^-1 . Vt
		but more efficient in computation time O(m n^2)
	*/
	for(i=0; i<Tmp.r; i++ )
		for(j=0; j<Tmp.c; j++ )
		{
			Tmp.val[i][j]=0;
			
			for(k=0; k<Tmp.r; k++ )
				Tmp.val[i][j] += V.val[i][k]*U.val[j][k]/w.val[0][k];

		}
	
	
	
	/*
		From there, estimate the beta's
	*/
	beta = MatrixProduct( &Tmp , Y, 0 );
	
	
	FreeMat(&U);
	FreeMat(&V);
	FreeMat(&w);
	FreeMat(&Tmp);

	return beta;

}

