
/*
 * USAGE. The function is called from Matlab as follows
 *
 *      [P, Pinv] = nojd_in(B, kmax, epsilon);
 *
 * Input:
 * 	 B     : is m-by-(m*n) matrix such that B = [B1 B2 ... Bn]
 *           where B1, ..., Bn are target matrices and expected to
 *           be of the form Bi \approx Q Li Q^{-1}, where Li is a
 *           diagonal matrix and Q is the common eigenvectors matrix
 *   kmax  : is the maximum number of sweeps
 *   epsilon : is the treshold for stopping criterion
 *
 * IMPORTANT: Current implementation works correctly iff P is initialized
 *            with the identity matrix! Modifications are necessary if 
 *            other initializations are needed.
 *
 * Output:
 *   P     : P = Q', where Q is the common eigenvectors matrix
 *   Pinv  : Pinv = (Q^{-1})', is the respective inverse
 *
 * Examples: "algorithms/helpers/nojd_in_cpp.m"
 */

/* DESCRIPTION. Given n matrices B1, B2, ..., Bn, each of size m-by-m,
 * the function looks for such non-orthogonal matrix P that all matrices
 *          A1 = P^{-1}*B1*P, ..., An = P^{-1}*Bn*P 
 * are (approximately) as diagonal as possible. For that, it 
 * (approximately)minimizes the sum of all squared off-diagonal elements, 
 * i.e. the objective
 *      obj(P) = \sum_{r=1^n} \sum_{1<=i \ne j<=m} ( Ar[i,j] )^2.
 * The minimization procedure is based on the itertive Jacobi (=Givens)
 * rotations and shear transformations. At each iteration, a closed form 
 * solution is computed for the respective Jacobi rotation and an 
 * approximate closed form solution is computed for the respecitive shear
 * transformation, making the algorithm fast.
 */

/*
 * ACKNOWLEDGMENTS. This is a C++/MEX-Matlab implementation of the 
 * non-orthogonal joint diagonalization by similarity algorithm (for
 * real matrices) as described in the paper:
 *   X. Luciani, L. Albera. Joint eigenvalue decomposition using polar 
 *     matrix factorization. In Proc. LVA/ICA, 2010.
 */

/* Copyright 2016, Anastasia Podosinnikova */

#include <math.h>
#include <mex.h>
#include <matrix.h>

void print_matrix(double** a, int m, int n){
  for (int i=0; i<m; i++){
    for (int j=0; j<n; j++)
      printf("%4.4f ", a[i][j]);
    printf("\n");
  }
}

/*////////////////////////////////////////////////////
/////////////////////// Updates //////////////////////
////////////////////////////////////////////////////*/
/* A (M-by-(N*M)) is N matrices each of size M-by-M */
void update_A_shear(double** A, double y, int M, int N, int p, int q){
  double c  = cosh(y), s  = sinh(y);
  int n,t;
  double apt, aqt, atp, atq;
  for (n=0; n<N; n++){
    for (t=0; t<M; t++){
      apt = A[p][M*n+t];
      aqt = A[q][M*n+t];
      A[p][M*n+t] =   c*apt - s*aqt;
      A[q][M*n+t] = - s*apt + c*aqt;
    }
  }
  for (n=0; n<N; n++){
    for (t=0; t<M; t++){
      atp = A[t][M*n+p];
      atq = A[t][M*n+q];
      A[t][M*n+p] = c*atp + s*atq;
      A[t][M*n+q] = s*atp + c*atq;
    }
  }
}

/* A (M-by-(N*M)) is N matrices each of size M-by-M */
void update_A_unitary(double** A, double theta, int M, int N, int p, int q){
  int n,t;
  double c=cos(theta), s = sin(theta);
  double apt, aqt, atp, atq;
  for (t=0; t<M; t++){
    for (n=0; n<N; n++){
      apt = A[p][M*n+t];
      aqt = A[q][M*n+t];
      A[p][M*n+t] = c*apt - s*aqt;
      A[q][M*n+t] = s*apt + c*aqt;
    }
  }
  for (n=0; n<N; n++){
    for (t=0; t<M; t++){
      atp = A[t][M*n+p];
      atq = A[t][M*n+q];
      A[t][M*n+p] = c*atp - s*atq;
      A[t][M*n+q] = s*atp + c*atq;
    }
  }
}

/* Q is M-by-M matrix */
void update_Q_shear(double** Q, int M, int p, int q, double y){
  double c = cosh(y), s = sinh(y);
  double Qtp, Qtq;
  for (int t=0; t<M; t++){
    Qtp = Q[t][p];
    Qtq = Q[t][q];
    Q[t][p] = c*Qtp + s*Qtq;
    Q[t][q] = s*Qtp + c*Qtq;
  }
}

/* Q is M-by-M matrix */
void update_Qinv_shear(double** Qinv, int M, int p, int q, double y){
  double c = cosh(y), s = sinh(y);
  double Qinvpt, Qinvqt;
  for (int t=0; t<M; t++){
    Qinvpt = Qinv[p][t];
    Qinvqt = Qinv[q][t];
    Qinv[p][t] =   c*Qinvpt - s*Qinvqt;
    Qinv[q][t] = - s*Qinvpt + c*Qinvqt;
  }
}

/* Q is M-by-M matrix */
void update_Q_unitary(double** Q, int M, int p, int q, double theta){
  double c = cos(theta), s = sin(theta);
  double Qtp,Qtq;
  for (int t=0; t<M; t++){
    Qtp = Q[t][p];
    Qtq = Q[t][q];
    Q[t][p] = c*Qtp - s*Qtq;
    Q[t][q] = s*Qtp + c*Qtq;
  }
}

/* Q is M-by-M matrix */
void update_Qinv_unitary(double** Qinv, int M, int p, int q, double theta){
  double c = cos(theta), s = sin(theta);
  double Qinvpt,Qinvqt;
  for (int t=0; t<M; t++){
    Qinvpt = Qinv[p][t];
    Qinvqt = Qinv[q][t];
    Qinv[p][t] = c*Qinvpt - s*Qinvqt;
    Qinv[q][t] = s*Qinvpt + c*Qinvqt;
  }
}

/*////////////////////////////////////////////////////
/////////////////////// Jacobi ///////////////////////
////////////////////////////////////////////////////*/
/* A (M-by-(N*M)) is N matrices each of size M-by-M */
double compute_theta_unitary(double** A, int M, int N, int p, int q){
  /* all these strange manipulations to actually compute atan correctly */
  /* computation of the Jacobi (=Givens) rotations */
  double g11 = 0, g22 = 0, g12 = 0, g21 = 0, g1p, g2p;
  for (int n=0; n<N; n++){
    g1p =   A[p][M*n+p] - A[q][M*n+q];
    g2p = - A[p][M*n+q] - A[q][M*n+p];
    g11 += g1p*g1p; g21 += g2p*g1p;
    g12 += g1p*g2p; g22 += g2p*g2p;
  }
  double ton = g11 - g22; 
  double toff = g12 + g21;
  double theta = 0.5 * atan2(toff, ton + sqrt(ton*ton+toff*toff));
  return theta;
}

/*////////////////////////////////////////////////////
////////////////////// JDTM //////////////////////////
////////////////////////////////////////////////////*/
/* A (M-by-(N*M)) is N matrices each of size M-by-M */
void compute_jdtm_coefficients(double** A, int M, int N, int p, int q, double* coefs){
  coefs[0] = 0; coefs[1] = 0; coefs[2] = 0;
  double d, o;
  for (int n=0; n<N; n++){
    d = A[p][M*n+p] - A[q][M*n+q]; o = A[p][M*n+q] - A[q][M*n+p];
    coefs[0] += d*d;
    coefs[1] += d*o;
    coefs[2] += o*o;
  }
}

double compute_y_shear_jdtm(double** A, int M, int N, int p, int q){
  double* coefs = new double[3];
  compute_jdtm_coefficients(A,M,N,p,q, coefs);
  double a=coefs[0], b=coefs[1], c=coefs[2];
  free(coefs);
  double lam = 0.5 * ( (c-a) + sqrt( (a+c)*(a+c) - 4*b*b ) );
  
  double y = 0.5 * atanh( (lam - c) / b );
  return y;
}


/*////////////////////////////////////////////////////
///////////////////// OBJS  //////////////////////////
////////////////////////////////////////////////////*/
void compute_frobenious_norm_coefficients(double** A, int M, int N, int p, int q, double* coefs){
  coefs[0]=0; coefs[1]=0; coefs[2]=0; coefs[3]=0;
  double api, aqi, aip, aiq, dd, oo;
  for (int n=0; n<N; n++){
    for (int i=0; i<M; i++){
      if ( (i!=p) && (i!=q) ){
        api=A[p][M*n+i]; aqi=A[q][M*n+i]; aip=A[i][M*n+p]; aiq=A[i][M*n+q];
        coefs[0] += api*api + aqi*aqi + aip*aip + aiq*aiq;
        coefs[1] += api*aqi - aip*aiq;
      }
    }
    dd = A[p][M*n+p]-A[q][M*n+q]; oo = A[p][M*n+q]-A[q][M*n+p];
    coefs[2] += oo*dd;
    coefs[3] += oo*oo + dd*dd;
  }
}

/* computes \sum_n ||A'_n||_F^2 - ||A_n||_F^2 (to insure descent for the shear transform) */
double compute_norm_descent(double** A, int M, int N, int p, int q, double y){
  double s = sinh(2*y), c = cosh(2*y);
  double* coefs = new double[4];
  compute_frobenious_norm_coefficients(A,M,N,p,q, coefs);
  double res = coefs[0]*(c-1) - coefs[1]*2*s + coefs[2]*2*c*s + coefs[3]*s*s;
  free(coefs);
  return res/N;
}

/*////////////////////////////////////////////////////
///////////////////// MAIN  //////////////////////////
////////////////////////////////////////////////////*/
/* A (M-by-(N*M)) is N matrices each of size M-by-M */
void diagonalization(double** A, int M, int N, double** Q, double** Qinv, double eps, int kmax){
/* CONVERGENCE TRACKING: , double* objs_fro, double* objs_off */
 
  double y, delta, theta, off, fro;
  bool encore = true;
  int k=0; /* sweep */
  //int iter=0; /* a small iteration, i.e. k,p,q */
  while (encore){ /* next sweep */
    encore = false;
    for (int p=0; p<M-1; p++){
      for (int q=p+1; q<M; q++){
        /* the shear transform */
        y = compute_y_shear_jdtm(A,M,N,p,q);
        delta = compute_norm_descent(A,M,N,p,q,y);
        if (delta < 0){
          update_A_shear(A,y,M,N,p,q);
          update_Q_shear(Q,M,p,q,y);
          update_Qinv_shear(Qinv,M,p,q,y);
        }
        
        /* the unitary transform */
        theta = compute_theta_unitary(A,M,N,p,q);
        update_A_unitary(A, theta, M, N, p, q);
        update_Q_unitary(Q,M,p,q,theta);
        update_Qinv_unitary(Qinv,M,p,q,theta);
        
        if (encore || (fabs(sin(theta))>eps))  {encore=true;}
        if (encore || (fabs(sinh(y))>eps))     {encore=true;}
      }
      
    }
    k++;
    if (k > kmax) {break;}
  }
  
}


/* Input: A, kmax, eps */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
   /* Check the number of inputs/ outputs. */
  if (nrhs!=3) {mexErrMsgTxt("Wrong number of inputs.");}
  if (nlhs>2)  {mexErrMsgTxt("Wrong number of outputs.");}
  
  double** A;
  int NM = mxGetN(prhs[0]); /* number of cols */
  int M  = mxGetM(prhs[0]); /* number of rows */
  A = (double**)mxCalloc(M,sizeof(double*));
  for (int m=0; m<M; m++)
    A[m] = (double*)mxCalloc(NM,sizeof(double));
  
  for (int m=0; m<M; m++){
    for (int nm=0; nm<NM; nm++)
      A[m][nm] = ((double*)mxGetData(prhs[0]))[m+nm*M];
  }
  
  int kmax = (int)mxGetScalar(prhs[1]);
  double eps = (double)mxGetScalar(prhs[2]);
  
  int N = NM/M;
  
  plhs[0] = mxCreateNumericMatrix(M, M, mxDOUBLE_CLASS, mxREAL);
  double* QPr = (double*)mxGetData(plhs[0]);
  
  /* !!! Note: this, actually, is the transpose of what we want. 
               as Matlab and C++ store matrices differently*/
  double** Q = new double*[M];
  for (int m=0; m<M; m++)
    Q[m] = QPr + M*m;
  for (int m=0; m<M; m++)
    Q[m][m] = 1;
  
  plhs[1] = mxCreateNumericMatrix(M, M, mxDOUBLE_CLASS, mxREAL);
  double* QinvPr = (double*)mxGetData(plhs[1]);
  double** Qinv = new double*[M];
  for (int m=0; m<M; m++)
    Qinv[m] = QinvPr + M*m;
  for (int m=0; m<M; m++)
    Qinv[m][m] = 1;
  
  diagonalization(A,M,N,Q,Qinv,eps,kmax);
  
  for (int m=0; m<M; m++)
    mxFree(A[m]);
  mxFree(A);
  
}
