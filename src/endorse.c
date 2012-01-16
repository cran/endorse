#include <string.h>
#include <stdio.h>      
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"

void R2endorse(/*   Data   */
		int *dY,          /* ordinal outcome variable: 0, 1,
				     ..., J-1. Length N * J */
		int *dT,          /* endorsement matrix. Length N * J */
		double *dZ,      /* covariate matirx, length N * M */
		/*   Data Structure   */
		int *n_samp,     /* # of obs: N */ 
		int *n_pol,     /* # of policies: J */
		int *n_dim,     /* dimension of covariates  */
		int *n_act,      /* # of actors: K */
		int *n_cat,      /* # of categories for each questions:
				    L_j. Length  J*/
		int *max_n_cat,  /* max # of categories */
		/*   Starting values   */
		double *X,      /* vector of ideal points. 
				    No need to include constant
				    Length N */
		double *dS,       /* matrix of support, s_{ij}(k)
				     Length (N * J) */
		double *sbeta,    /* (alpha_j, beta_j), length J * 2 */
		double *stau,     /* cut points for each question;
				     the first cut point is set to 0
				     and the last one set to tau_{L_j - 1}+1000:
				     Length J * max {L_j} */
		double *slambda,  /* lambda, length (J * M) * K  */
		double *somega2,   /* omega2: variance of s given lambda,
				     length J * K */
		double *dtheta,   /* theta, vector of length K * M  */
		double *phi2,     /* phi, vector of length K * M */
		double *delta,    /* delta, vector of length M  */
		/*  Prior parameters  */
		double *dmu_beta,   /* prior mean of factor loadings */
		double *dA0_beta,     /* prior precision of factor loadings,
				    length 2 * 2: can be set to zero to
				    induce improper prior for beta alone */
		double *dA0_x,     /* known prior precision of ideal points, identical for all i */
		double *dmu_theta, /* prior mean for theta, length K */
		double *dA0_theta, /* prior precision of theta, identical for all k */
		double *mu_delta,  /* prior mean of delta, length M */
		double *dA0_delta,  /* prior precision of delta, length M * M */
		double *s0_omega2,      /* prior scale for InvChi2 of omega */
		int *nu0_omega2,         /* prior d.f. for InvChi2 of omega */
		double *s0_phi2,    /* prior scale for InvChi2 of phi */
		int *nu0_phi2,      /* prior d.f. for InvChi2 of phi */
		/* MCMC settings */
		int *n_gen,       /* # of gibbs draws */
		int *burn,       /* # of burnin period */
		int *thin,       /* thinning */
		int *mda,        /* marginal data augmentation? */
		int *mh,         /* Metropolis-Hasting step? */
		double *prop,    /* proposal variance for MH step */
		int *accept,     /* acceptance counter */
		int *x_sd,       /* Output of x: sd if 1, sample if 0 */
		int *tau_out,
		double *betaStore,
		double *tauStore,
		double *xStore,
		double *lambdaStore,
		double *thetaStore,
		double *deltaStore
		 ){
  /* loop counters */
  int i, j, k, l, m, n, main_loop, itemp;
  int ibeta = 0, itau = 0, ix = 0, ilambda = 0, itheta = 0, idelta = 0, keep = 1;
  double varx;

  /* storage vectors */
  int *Y_j = intArray(*n_samp);
  double *beta = doubleArray(2);
  double *mu_beta = doubleArray(2);
  double *s = doubleArray(1);
  double *var_epsilon = doubleArray(1);
  var_epsilon[0] = 1;
  double *mu_s = doubleArray(1);
  double *x_i = doubleArray(1);
  double *mu_x = doubleArray(1);
  double *lambda_jk = doubleArray(*n_dim);
  double *omega2_jk = doubleArray(1);
  double *mu_lambda = doubleArray(*n_dim);
  double *stemp = doubleArray(*n_samp);
  double *ztemp = doubleArray(*n_samp * *n_dim);
  double *theta_km = doubleArray(*n_dim);
  double *phi2_km = doubleArray(1);
  double *mu_theta = doubleArray(1);
  double *sig2_x = doubleArray(1);

  /* storage matrices */
  /** data **/
  int **T = intMatrix(*n_samp, *n_pol); /* treatment */
  int **Y = intMatrix(*n_samp, *n_pol); /* response */
  double **Z = doubleMatrix(*n_samp, *n_dim); /* covariates */

  /** starting values **/
  double **S = doubleMatrix(*n_samp, *n_pol); /* support parameter */
  double **Beta = doubleMatrix(*n_pol, 2); /* alpha and beta */
  double **Tau = doubleMatrix(*n_pol, *max_n_cat); /* cut point */
  double ***Lambda = doubleMatrix3D(*n_pol, *n_act, *n_dim); /* lambda */
  double **Omega2 = doubleMatrix(*n_pol, *n_act); /* omega2 */
  double **Theta = doubleMatrix(*n_act, *n_dim); /* theta */

  /** prior mean and precision **/
  double **Mu_beta = doubleMatrix(*n_pol, 2);
  double **A0_beta = doubleMatrix(2, 2);
  double **A0_s = doubleMatrix(1, 1);
  double **A0_x = doubleMatrix(1, 1);
  double **A0_lambda = doubleMatrix(*n_dim, *n_dim);
  double **A0_theta = doubleMatrix(1, 1);
  double **A0_delta = doubleMatrix(*n_dim, *n_dim);

  /** matrices for regression **/
  double **Ystar = doubleMatrix(*n_samp, *n_pol); /* latent random utility */
  double **U = doubleMatrix(*n_samp+2, 3); /* x_{i} + s_{ij}(k) */
  double **D_s = doubleMatrix(2, 2);
  double **D_x = doubleMatrix(*n_pol+1, 2);
  double **D_theta = doubleMatrix(*n_pol+1, 2);
  double **D_delta = doubleMatrix(*n_samp + *n_dim, *n_dim + 1);

  /* get random seed */
  GetRNGstate();

  /* packing data */
  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for (i = 0; i < *n_samp; i++)
      Y[i][j] = dY[itemp++];

  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for(i = 0; i < *n_samp; i++)
      T[i][j] = dT[itemp++];

  itemp = 0;
  for (m = 0; m < *n_dim; m++)
    for (i = 0; i < *n_samp; i++)
      Z[i][m] = dZ[itemp++];

  /* packing starting values */
  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for (i = 0; i < *n_samp; i++)
      S[i][j] = dS[itemp++];

  itemp = 0;
  for (m = 0; m < 2; m++)
    for (j = 0; j < *n_pol; j++)
      Beta[j][m] = sbeta[itemp++];

  itemp = 0;
  for (l = 0; l < *max_n_cat; l++)
    for (j = 0; j < *n_pol; j++)
      Tau[j][l] = stau[itemp++];

  itemp = 0;
  for (k = 0; k < *n_act; k++)
    for (j = 0; j < *n_pol; j++)
      for (m = 0; m < *n_dim; m++)
      Lambda[j][k][m] = slambda[itemp++];

  itemp = 0;
  for (k = 0; k < *n_act; k++)
    for (j = 0; j < *n_pol; j++)
      Omega2[j][k] = somega2[itemp++];

  itemp = 0;
  for (m = 0; m < *n_dim; m++)
    for (k = 0; k < *n_act; k++)
      Theta[k][m] = dtheta[itemp++];

  /* packing prior mean */
  itemp = 0;
  for (m = 0; m < 2; m++)
    for (j = 0; j < *n_pol; j++)
      Mu_beta[j][m] = dmu_beta[itemp++];

  /* packing prior precision */
  itemp = 0;
  for (m = 0; m < 2; m++)
    for (n = 0; n < 2; n++)
      A0_beta[n][m] = dA0_beta[itemp++];

  for (m = 0; m < *n_dim; m++)
    for (n = 0; n < *n_dim; n++)
      A0_lambda[m][n] = 0;

  itemp = 0;
  for (m = 0; m < *n_dim; m++)
    for (n = 0; n < *n_dim; n++)
      A0_delta[n][m] = dA0_delta[itemp++];

  A0_x[0][0] = *dA0_x;
  sig2_x[0] = (1 / *dA0_x);

  /* Gibbs Sampler */
  for (main_loop = 1; main_loop <= *n_gen; main_loop++) {
    if (main_loop == 1)
      Rprintf("Start Gibbs Sampler\n");

/*     if (main_loop % 100 == 0) */
      Rprintf("%6d / %6d\n", main_loop, *n_gen);

    /** start sampling alpha and beta **/
    for (j = 0; j < *n_pol; j++) {

      /*** vectors ***/
      double *tau = doubleArray(n_cat[j]);
      double *MHprop = doubleArray(n_cat[j]-2);

      /*** proposal variance vector ***/
      for (l = 0; l < (n_cat[j]-2); l++)
	MHprop[l] = *prop;

      /*** response of each question ***/
      for (i = 0; i < *n_samp; i++)
	Y_j[i] = Y[i][j];
      
      /*** systematic component of utility ***/
      for (i = 0; i < *n_samp; i++) {
	U[i][0] = -1;
	U[i][1] = X[i] + S[i][j];
      }

      /*** starting values of alpha and beta ***/
      for (m = 0; m < 2; m++)
	beta[m] = Beta[j][m];

      /*** starting values of tau ***/
      for (l = 0; l < n_cat[j]; l++)
	tau[l] = Tau[j][l];

      /*** prior mean ***/
      for (m = 0; m < 2; m++)
	mu_beta[m] = Mu_beta[j][m];
      
      bsupportoprobitMCMC(Y_j, U, beta, tau, *n_samp, 2, n_cat[j], 1,
			  mu_beta, A0_beta, *mda, *mh, MHprop, accept, 1);

      /*** packing sampled alpha and beta ***/
      for (m = 0; m < 2; m++)
	Beta[j][m] = beta[m];

      /*** packing sampled tau ***/
      for (l = 0; l < n_cat[j]; l++)
	Tau[j][l] = tau[l];
      
      /*** packing sampled latent random utility ***/
      for (i = 0; i < *n_samp; i++)
	Ystar[i][j] = U[i][2];

      /*** storing alpha and beta ***/
      if(main_loop > *burn) {
	if(keep == *thin) {

	  for (m = 0; m < 2; m++)
	    betaStore[ibeta++] = Beta[j][m];

	  if (*tau_out) {
	    for (l = 0; l < *max_n_cat-1; l++)
	      tauStore[itau++] = Tau[j][l];
	  }

	}
      }

      R_FlushConsole();
      R_CheckUserInterrupt();
    } /** end of sampling alpha and beta **/


    /** start sampling s **/
    for (i = 0; i < *n_samp; i++){
      for (j = 0; j < *n_pol; j++){

	k = T[i][j];
	
	/*** if not control, sample s ***/
	if (k > 0) {

	  D_s[0][0] = Beta[j][1];
	  D_s[0][1] = Ystar[i][j] + Beta[j][0] - Beta[j][1]*X[i];

	  mu_s[0] = 0;
	  for (m = 0; m < *n_dim; m++)
	    mu_s[0] += Lambda[j][k-1][m] * Z[i][m];

	  A0_s[0][0] = (1 / Omega2[j][k-1]);

	  bNormalReg(D_s, s, var_epsilon, 1, 1, 1, 1, mu_s, A0_s,
		     1, 1, 1, 1);

	  /*** packing sampled s ***/
	  S[i][j] = s[0];
	}

	R_FlushConsole();
	R_CheckUserInterrupt();
      }
    }/** end of sampling s **/

    /** start sampling x except the first and the last respondents **/
    for (i = 1; i < (*n_samp-1); i++) {

      for (j = 0; j < *n_pol; j++) {
	D_x[j][0] = Beta[j][1];
	D_x[j][1] = Ystar[i][j] + Beta[j][0] - Beta[j][1]*S[i][j];
      }
      
      /*** prior mean ***/
      mu_x[0] = 0;
      for (m = 0; m < *n_dim; m++)
	mu_x[0] += Z[i][m] * delta[m];
      
      bNormalReg(D_x, x_i, var_epsilon, *n_pol, 1, 1, 1, mu_x, A0_x,
		 1, 1, 1, 1);

      /*** packing sampled x ***/
      X[i] = x_i[0];

      R_FlushConsole();
      R_CheckUserInterrupt();
    }

    /*** storing sampled x ***/
    if(main_loop > *burn) {
      if(keep == *thin) {
	if (*x_sd){
	  varx = var(X, *n_samp, 1);
	  xStore[ix++] = sqrt(varx);
	} else {
	  for (i = 0; i < *n_samp; i++)
	    xStore[ix++] = X[i];
	}
      }
    }
    /** end of sampling x **/


    /** start sampling lambda and omega2 **/
    for (j = 0; j < *n_pol; j++) {
      for (k = 0; k < *n_act; k++) {

	/*** # of observations with treatment k for question j ***/
	n = 0;
	itemp = 0;
	for (i = 0; i < *n_samp; i++) {
	  if (k+1 == T[i][j]) {
	    stemp[n] = S[i][j];

	    for (m = 0; m < *n_dim; m++)
	      ztemp[itemp++] = Z[i][m];

	    n++;
	  }
	}
	
	if (n == 0)
	  continue;
	
	double **D_lambda = doubleMatrix(n + *n_dim, *n_dim+1);

	itemp = 0;
	for (l = 0; l < n; l++) {

	  for (m = 0; m < *n_dim; m++)
	    D_lambda[l][m] = ztemp[itemp++];

	  D_lambda[l][*n_dim] = stemp[l];
	}

	for (m = 0; m < *n_dim; m++) {
	  mu_lambda[m] = Theta[k][m];
	  A0_lambda[m][m] = (1 / phi2[(k * *n_dim) + m]);
	}
	
	bNormalReg(D_lambda, lambda_jk, omega2_jk, n, *n_dim, 1, 1, mu_lambda,
		   A0_lambda, 1, *s0_omega2, *nu0_omega2, 0);

	for (m = 0; m < *n_dim; m++)
	  Lambda[j][k][m] = lambda_jk[m];

	Omega2[j][k] = omega2_jk[0];

	FreeMatrix(D_lambda, n + *n_dim);

	R_FlushConsole();
	R_CheckUserInterrupt();
      }
    } 

    /*** storing sampled lambda ***/
    if (main_loop > *burn) {
      if (keep == *thin) {
	for (j = 0; j < *n_pol; j++)
	  for (k = 0; k < *n_act; k++)
	    for (m = 0; m < *n_dim; m++)
	      lambdaStore[ilambda++] = Lambda[j][k][m];
      }
    }
    /** end of sampling lambda and omega2 **/

    /** start sampling theta and phi2 **/
    for (k = 0; k < *n_act; k++) {
      for (m = 0; m < *n_dim; m++) {

	for (j = 0; j < *n_pol; j++) {
	  D_theta[j][0] = 1;
	  D_theta[j][1] = Lambda[j][k][m];
	}

	/*** prior mean for theta_{k} ***/
	mu_theta[0] = dmu_theta[m];

	A0_theta[0][0] = dA0_theta[m];

	bNormalReg(D_theta, theta_km, phi2_km, *n_pol, 1, 1, 1,
		   mu_theta, A0_theta, 1, *s0_phi2, *nu0_phi2, 0);

	/*** packing sampled theta ***/
	Theta[k][m] = theta_km[0];

	/*** packing sampled phi2 ***/
	phi2[(k * *n_dim) + m] = phi2_km[0];

	/*** storing sampled theta ***/
	if (main_loop > *burn) {
	  if (keep == *thin) {
	    thetaStore[itheta++] = Theta[k][m];
	  }
	}

	R_FlushConsole();
	R_CheckUserInterrupt();
      }
    } /** end of sampling theta and phi2 **/

    /** Start sampling delta **/
    for (i = 0; i < *n_samp; i++) {
      for (m = 0; m < *n_dim; m++)
	D_delta[i][m] = Z[i][m];
      
      D_delta[i][*n_dim] = X[i];
    }

    bNormalReg(D_delta, delta, sig2_x, *n_samp, *n_dim, 1, 1,
	       mu_delta, A0_delta, 1, 1, 1, 1);

    if (main_loop > *burn) {
      if (keep == *thin) {
	for (m = 0; m < *n_dim; m++)
	  deltaStore[idelta++] = delta[m];
      }
    }
    /** end of sampling delta **/

    /** update thinning counter **/
    if (keep == *thin) {
      keep = 1;
    } else {
      keep++;
    }

    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  Rprintf("End Gibbs Sampler\n");

  PutRNGstate();


  /* freeing memory */
  free(Y_j);
  free(beta);
  free(mu_beta);
  free(s);
  free(var_epsilon);
  free(mu_s);
  free(x_i);
  free(mu_x);
  free(lambda_jk);
  free(omega2_jk);
  free(mu_lambda);
  free(stemp);
  free(ztemp);
  free(theta_km);
  free(phi2_km);
  free(mu_theta);
  free(sig2_x);

  FreeintMatrix(T, *n_samp);
  FreeintMatrix(Y, *n_samp);
  FreeMatrix(Z, *n_samp);
  FreeMatrix(S, *n_samp);
  FreeMatrix(Beta, *n_pol);
  FreeMatrix(Tau, *n_pol);
  FreeMatrix(Theta, *n_act);
  Free3DMatrix(Lambda, *n_pol, *n_act);
  FreeMatrix(Omega2, *n_pol);
  FreeMatrix(Mu_beta, *n_pol);
  FreeMatrix(A0_beta, 2);
  FreeMatrix(Ystar, *n_samp);
  FreeMatrix(U, *n_samp+2);
  FreeMatrix(D_s, 2);
  FreeMatrix(A0_s, 1);
  FreeMatrix(D_x, *n_pol+1);
  FreeMatrix(A0_x, 1);
  FreeMatrix(A0_lambda, 1);
  FreeMatrix(D_theta, *n_pol+1);
  FreeMatrix(A0_theta, 1);
  FreeMatrix(A0_delta, *n_dim);
  FreeMatrix(D_delta, *n_samp + *n_dim);
}


void R2endorseNoCov(/*   Data   */
		    int *dY,          /* ordinal outcome variable: 0, 1,
					 ..., J-1. Length N * J */
		    int *dT,          /* endorsement matrix. Length N * J */
		    /*   Data Structure   */
		    int *n_samp,     /* # of obs: N */ 
		    int *n_pol,     /* # of policies: J */
		    int *n_act,      /* # of actors: K */
		    int *n_cat,      /* # of categories for each questions:
					L_j. Length  J*/
		    int *max_n_cat,  /* max # of categories */
		    /*   Starting values   */
		    double *X,      /* vector of ideal points. 
				       No need to include constant
				       Length N */
		    double *dS,       /* matrix of support, s_{ij}(k)
					 Length (N * J) */
		    double *sbeta,    /* (alpha_j, beta_j), length J * 2 */
		    double *stau,     /* cut points for each question;
					 the first cut point is set to 0
					 and the last one set to tau_{J-1}+1000:
					 Length J * max {L_j} */
		    double *slambda,  /* lambda, length J * K  */
		    double *somega2,   /* omega2: variance of s given lambda,
					  length J * K */
		    double *theta,   /* theta, vector of length K  */
		    double *phi2,     /* phi, vector of length K */
		    /*  Prior parameters  */
		    double *dmu_beta,   /* prior mean of factor loadings */
		    double *dA0_beta,     /* prior precision of factor loadings,
					     length 2 * 2: can be set to zero to
					     induce improper prior for beta alone */
		    double *dmu_x,     /* prior mean of ideal points. length N */
		    double *dA0_x,     /* prior precision of ideal points, identical for all i */
		    double *dmu_theta, /* prior mean for theta, length K */
		    double *dA0_theta, /* prior precision of theta, identical for all k */
		    double *s0_omega2,      /* prior scale for InvChi2 of omega */
		    int *nu0_omega2,         /* prior d.f. for InvChi2 of omega */
		    double *s0_phi2,    /* prior scale for InvChi2 of phi */
		    int *nu0_phi2,      /* prior d.f. for InvChi2 of phi */
		    /* MCMC settings */
		    int *n_gen,       /* # of gibbs draws */
		    int *burn,       /* # of burnin period */
		    int *thin,       /* thinning */
		    int *mda,        /* marginal data augmentation? */
		    int *mh,         /* Metropolis-Hasting step? */
		    double *prop,    /* proposal variance for MH step */
		    int *accept,     /* acceptance counter */
		    int *x_sd,
		    int *tau_out,
		    double *betaStore,
		    double *tauStore,
		    double *xStore,
		    double *lambdaStore,
		    double *thetaStore
		    ){
  /* loop counters */
  int i, j, k, l, m, n, main_loop, itemp;
  int ibeta = 0, itau = 0, ix = 0, ilambda = 0, itheta = 0, keep = 1;
  double varx;

  /* storage vectors */
  int *Y_j = intArray(*n_samp);
  double *beta = doubleArray(2);
  double *mu_beta = doubleArray(2);
  double *s = doubleArray(1);
  double *var_epsilon = doubleArray(1);
  var_epsilon[0] = 1;
  double *mu_s = doubleArray(1);
  double *x_i = doubleArray(1);
  double *mu_x = doubleArray(1);
  double *lambda_jk = doubleArray(1);
  double *omega2_jk = doubleArray(1);
  double *mu_lambda = doubleArray(1);
  double *stemp = doubleArray(*n_samp);
  double *theta_k = doubleArray(1);
  double *phi2_k = doubleArray(1);
  double *mu_theta = doubleArray(1);

  /* storage matrices */
  /** data **/
  int **T = intMatrix(*n_samp, *n_pol); /* treatment */
  int **Y = intMatrix(*n_samp, *n_pol); /* response */

  /** starting values **/
  double **S = doubleMatrix(*n_samp, *n_pol); /* support parameter */
  double **Beta = doubleMatrix(*n_pol, 2); /* alpha and beta */
  double **Tau = doubleMatrix(*n_pol, *max_n_cat); /* cut point */
  double **Lambda = doubleMatrix(*n_pol, *n_act); /* lambda */
  double **Omega2 = doubleMatrix(*n_pol, *n_act); /* omega2 */

  /** prior mean and precision **/
  double **Mu_beta = doubleMatrix(*n_pol, 2);
  double **A0_beta = doubleMatrix(2, 2);
  double **A0_s = doubleMatrix(1, 1);
  double **A0_x = doubleMatrix(1, 1);
  double **A0_lambda = doubleMatrix(1, 1);
  double **A0_theta = doubleMatrix(1, 1);

  /** matrices for regression **/
  double **Ystar = doubleMatrix(*n_samp, *n_pol); /* latent random utility */
  double **U = doubleMatrix(*n_samp+2, 3); /* x_{i} + s_{ij}(k) */
  double **D_s = doubleMatrix(2, 2);
  double **D_x = doubleMatrix(*n_pol+1, 2);
  double **D_theta = doubleMatrix(*n_pol+1, 2);

  /* get random seed */
  GetRNGstate();

  /* packing data */
  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for (i = 0; i < *n_samp; i++)
      Y[i][j] = dY[itemp++];

  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for(i = 0; i < *n_samp; i++)
      T[i][j] = dT[itemp++];

  /* packing starting values */
  itemp = 0;
  for (j = 0; j < *n_pol; j++)
    for (i = 0; i < *n_samp; i++)
      S[i][j] = dS[itemp++];

  itemp = 0;
  for (m = 0; m < 2; m++)
    for (j = 0; j < *n_pol; j++)
      Beta[j][m] = sbeta[itemp++];

  itemp = 0;
  for (l = 0; l < *max_n_cat; l++)
    for (j = 0; j < *n_pol; j++)
      Tau[j][l] = stau[itemp++];

  itemp = 0;
  for (k = 0; k < *n_act; k++)
    for (j = 0; j < *n_pol; j++)
      Lambda[j][k] = slambda[itemp++];

  itemp = 0;
  for (k = 0; k < *n_act; k++)
    for (j = 0; j < *n_pol; j++)
      Omega2[j][k] = somega2[itemp++];

  /* packing prior mean */
  itemp = 0;
  for (m = 0; m < 2; m++)
    for (j = 0; j < *n_pol; j++)
      Mu_beta[j][m] = dmu_beta[itemp++];

  /* packing prior precision */
  itemp = 0;
  for (m = 0; m < 2; m++)
    for (n = 0; n < 2; n++)
      A0_beta[n][m] = dA0_beta[itemp++];

  A0_x[0][0] = *dA0_x;

  A0_theta[0][0] = *dA0_theta;


  /* Gibbs Sampler */
  for (main_loop = 1; main_loop <= *n_gen; main_loop++) {
    if (main_loop == 1)
      Rprintf("Start Gibbs Sampler\n");

    /*     if (main_loop % 100 == 0) */
    Rprintf("%6d / %6d\n", main_loop, *n_gen);

    /** start sampling alpha and beta **/
    for (j = 0; j < *n_pol; j++) {

      /*** vectors ***/
      double *tau = doubleArray(n_cat[j]);
      double *MHprop = doubleArray(n_cat[j]-2);

      /*** proposal variance vector ***/
      for (l = 0; l < n_cat[j]-2; l++)
	MHprop[l] = *prop;

      /*** response of each question ***/
      for (i = 0; i < *n_samp; i++)
	Y_j[i] = Y[i][j];
      
      /*** systematic component of utility ***/
      for (i = 0; i < *n_samp; i++) {
	U[i][0] = -1;
	U[i][1] = X[i] + S[i][j];
      }

      /*** starting values of alpha and beta ***/
      for (m = 0; m < 2; m++)
	beta[m] = Beta[j][m];

      /*** starting values of tau ***/
      for (l = 0; l < n_cat[j]; l++)
	tau[l] = Tau[j][l];

      /*** prior mean ***/
      for (m = 0; m < 2; m++)
	mu_beta[m] = Mu_beta[j][m];
      
      bsupportoprobitMCMC(Y_j, U, beta, tau, *n_samp, 2, n_cat[j], 1,
			  mu_beta, A0_beta, *mda, *mh, MHprop, accept, 1);

      /*** packing sampled alpha and beta ***/
      for (m = 0; m < 2; m++)
	Beta[j][m] = beta[m];

      /*** packing sampled tau ***/
      for (l = 0; l < n_cat[j]; l++)
	Tau[j][l] = tau[l];
      
      /*** packing sampled latent random utility ***/
      for (i = 0; i < *n_samp; i++)
	Ystar[i][j] = U[i][2];

      /*** storing alpha and beta ***/
      if(main_loop > *burn) {
	if(keep == *thin) {

	  for (m = 0; m < 2; m++)
	    betaStore[ibeta++] = Beta[j][m];

	  if (*tau_out) {
	    for (l = 0; l < *max_n_cat-1; l++)
	      tauStore[itau++] = Tau[j][l];
	  }

	}
      }

      R_FlushConsole();
      R_CheckUserInterrupt();
    } /** end of sampling alpha and beta **/


    /** start sampling s **/
    for (i = 0; i < *n_samp; i++){
      for (j = 0; j < *n_pol; j++){

	k = T[i][j];
	
	/*** if not control, sample s ***/
	if (k > 0) {

	  D_s[0][0] = Beta[j][1];
	  D_s[0][1] = Ystar[i][j] + Beta[j][0] - Beta[j][1]*X[i];
	  
	  mu_s[0] = Lambda[j][k-1];
	  A0_s[0][0] = (1 / Omega2[j][k-1]);

	  bNormalReg(D_s, s, var_epsilon, 1, 1, 1, 1, mu_s, A0_s,
		     1, 1, 1, 1);

	  /*** packing sampled s ***/
	  S[i][j] = s[0];
	}

	R_FlushConsole();
	R_CheckUserInterrupt();
      }
    }/** end of sampling s **/

    /** start sampling x except the first and the last respondents **/
    for (i = 1; i < (*n_samp-1); i++) {

      for (j = 0; j < *n_pol; j++) {
	D_x[j][0] = Beta[j][1];
	D_x[j][1] = Ystar[i][j] + Beta[j][0] - Beta[j][1]*S[i][j];
      }
      
      /*** prior mean ***/
      mu_x[0] = dmu_x[i];
      
      bNormalReg(D_x, x_i, var_epsilon, *n_pol, 1, 1, 1, mu_x, A0_x,
		 1, 1, 1, 1);

      /*** packing sampled x ***/
      X[i] = x_i[0];

      R_FlushConsole();
      R_CheckUserInterrupt();
    }

    /*** storing sampled x ***/
    if(main_loop > *burn) {
      if(keep == *thin) {
	if (*x_sd) {
	  varx = var(X, *n_samp, 1);
	  xStore[ix++] = sqrt(varx);
	} else {
	  for (i = 0; i < *n_samp; i++)
	    xStore[ix++] = X[i];
	}
      }
    }
    /** end of sampling x **/


    /** start sampling lambda omega2 **/
    for (j = 0; j < *n_pol; j++) {
      for (k = 0; k < *n_act; k++) {

	/*** # of observations with treatment k for question j ***/
	n = 0;
	for (i = 0; i < *n_samp; i++) {
	  if (k+1 == T[i][j]) {
	    stemp[n] = S[i][j];
	    n++;
	  }
	}
	
	if (n == 0)
	  continue;
	
	double **D_lambda = doubleMatrix(n+1, 2);

	for (l = 0; l < n; l++) {
	  D_lambda[l][0] = 1;
	  D_lambda[l][1] = stemp[l];
	}

	mu_lambda[0] = theta[k];
	A0_lambda[0][0] = (1 / phi2[k]);
	
	bNormalReg(D_lambda, lambda_jk, omega2_jk, n, 1, 1, 1, mu_lambda,
		   A0_lambda, 1, *s0_omega2, *nu0_omega2, 0);

	Lambda[j][k] = lambda_jk[0];
	Omega2[j][k] = omega2_jk[0];

	FreeMatrix(D_lambda, n+1);

	R_FlushConsole();
	R_CheckUserInterrupt();
      }
    } 

    /*** storing sampled lambda ***/
    if (main_loop > *burn) {
      if (keep == *thin) {
	for (j = 0; j < *n_pol; j++)
	  for (k = 0; k < *n_act; k++)
	    lambdaStore[ilambda++] = Lambda[j][k];
      }
    }
    /** end of sampling lambda and omega2 **/


    /** start sampling theta and phi2 **/
    for (k = 0; k < *n_act; k++) {

      for (j = 0; j < *n_pol; j++) {
	D_theta[j][0] = 1;
	D_theta[j][1] = Lambda[j][k];
      }

      /*** prior mean for theta_{k} ***/
      mu_theta[0] = dmu_theta[k];

      bNormalReg(D_theta, theta_k, phi2_k, *n_pol, 1, 1, 1,
		 mu_theta, A0_theta, 1, *s0_phi2, *nu0_phi2, 0);

      /*** packing sampled theta ***/
      theta[k] = theta_k[0];

      /*** packing sampled phi2 ***/
      phi2[k] = phi2_k[0];

      /*** storing sampled theta ***/
      if (main_loop > *burn) {
	if (keep == *thin) {
	  thetaStore[itheta++] = theta[k];
	}
      }

      R_FlushConsole();
      R_CheckUserInterrupt();
    } /** end of sampling theta and phi2 **/


    /** update thinning counter **/
    if (keep == *thin) {
      keep = 1;
    } else {
      keep++;
    }


    R_FlushConsole();
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */
  Rprintf("End Gibbs Sampler\n");

  PutRNGstate();


  /* freeing memory */
  free(Y_j);
  free(beta);
  free(mu_beta);
  free(s);
  free(var_epsilon);
  free(mu_s);
  free(x_i);
  free(mu_x);
  free(lambda_jk);
  free(omega2_jk);
  free(mu_lambda);
  free(stemp);
  free(theta_k);
  free(phi2_k);
  free(mu_theta);

  FreeintMatrix(T, *n_samp);
  FreeintMatrix(Y, *n_samp);
  FreeMatrix(S, *n_samp);
  FreeMatrix(Beta, *n_pol);
  FreeMatrix(Tau, *n_pol);
  FreeMatrix(Lambda, *n_pol);
  FreeMatrix(Omega2, *n_pol);
  FreeMatrix(Mu_beta, *n_pol);
  FreeMatrix(A0_beta, 2);
  FreeMatrix(Ystar, *n_samp);
  FreeMatrix(U, *n_samp+2);
  FreeMatrix(D_s, 2);
  FreeMatrix(A0_s, 1);
  FreeMatrix(D_x, *n_pol+1);
  FreeMatrix(A0_x, 1);
  FreeMatrix(A0_lambda, 1);
  FreeMatrix(D_theta, *n_pol+1);
  FreeMatrix(A0_theta, 1);
}
