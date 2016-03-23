/* Statistical routines for Gumbel (type I extreme value) distributions.
 * 
 * Contents:
 *   1. Routine for evaluating densities and distributions
 *   2. Generic API routines: for general interface w/ histogram module
 *   3. Routines for dumping plots to files
 *   4. Routines for sampling (requires random module)
 *   5. Maximum likelihood fitting to data (requires minimizer module)
 *   6. Stats driver
 *   7. Unit tests
 *   8. Test driver
 *   9. Example
 *  10. Copyright and license information
 * 
 * Note: SRE, Mon Aug  6 13:42:09 2007
 * ML fitting routines will be prone to over/underfitting 
 * problems for scores outside a "normal" range, because
 * of exp(-lambda * x) calls. The Lawless ML estimation
 * may eventually need to be recast in log space.
 */
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "extern/easel/easel.h"
#include "extern/easel/esl_stats.h"
#include "extern/easel/esl_gumbel.h"


/*****************************************************************
 * 5. Routines for maximum likelihood fitting Gumbels to data
 * (fitting truncated distributions requires augmentation w/ minimizer module)
 *****************************************************************/ 

/*****************************************************************
 * Complete data, maximum a posteriori parameters
 *****************************************************************/ 

/* lawless416()
 * SRE, Thu Nov 13 11:48:50 1997 [St. Louis]
 * 
 * Purpose:  Equation 4.1.6 from [Lawless82], pg. 143, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to Gumbel lambda parameter.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *           
 * Args:     x      - array of sample values 
 *           n      - number of samples 
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.1.6 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.1.6 evaluated at lambda
 *           
 * Return:   (void)
 */ 
static void
lawless416(double *x, int n, double lambda, double *ret_f, double *ret_df)
{
  double esum;			/* \sum e^(-lambda xi)      */
  double xesum;			/* \sum xi e^(-lambda xi)   */
  double xxesum;		/* \sum xi^2 e^(-lambda xi) */
  double xsum;			/* \sum xi                  */
  int i;

  esum = xesum = xsum  = xxesum = 0.;
  for (i = 0; i < n; i++)
    {
      xsum   += x[i];
      xesum  += x[i] * exp(-1. * lambda * x[i]);
      xxesum += x[i] * x[i] * exp(-1. * lambda * x[i]);
      esum   += exp(-1. * lambda * x[i]);
    }
  *ret_f  = (1./lambda) - (xsum / n)  + (xesum / esum);
  *ret_df = ((xesum / esum) * (xesum / esum))
    - (xxesum / esum)
    - (1. / (lambda * lambda));
}

/* Function: esl_gumbel_FitComplete()
 * Synopsis: Estimates $\mu$, $\lambda$ from complete data.
 * Date:     SRE, Fri Nov 14 07:56:29 1997 [St. Louis] - HMMER's EVDMaxLikelyFit()
 * 
 * Purpose:  Given an array of Gumbel-distributed samples <x[0]..x[n-1]>,
 *           find maximum likelihood parameters <mu> and <lambda>.
 *           
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *            for lambda using Newton/Raphson iterations,
 *            then substitutes lambda into Lawless' equation 4.1.5
 *            to get mu. 
 *           
 * Args:     x          - list of Gumbel distributed samples
 *           n          - number of samples
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *           
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslENOHALT> if the fit doesn't converge.
 */
int
esl_gumbel_FitComplete(double *x, int n, double *ret_mu, double *ret_lambda)
{
  double  variance;
  double  lambda, mu;
  double  fx;			/* f(x)  */
  double  dfx;			/* f'(x) */
  double  esum;                 /* \sum e^(-lambda xi) */ 
  double  tol = 1e-5;
  int     i;

  /* 1. Find an initial guess at lambda
   *    (Evans/Hastings/Peacock, Statistical Distributions, 2000, p.86)
   */
  esl_stats_DMean(x, n, NULL, &variance);
  lambda = eslCONST_PI / sqrt(6.*variance);

  /* 2. Use Newton/Raphson to solve Lawless 4.1.6 and find ML lambda
   */
  for (i = 0; i < 100; i++)
    {
      lawless416(x, n, lambda, &fx, &dfx);
      if (fabs(fx) < tol) break;             /* success */
      lambda = lambda - fx / dfx;	     /* Newton/Raphson is simple */
      if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

  /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
   *      Resort to a bisection search. Worse convergence speed
   *      but guaranteed to converge (unlike Newton/Raphson).
   *      We assume that fx is a monotonically decreasing function of x;
   *      i.e. fx > 0 if we are left of the root, fx < 0 if we
   *      are right of the root.
   */ 
  if (i == 100)
    {
      double left, right, mid;
      ESL_DPRINTF1(("esl_gumbel_FitComplete(): Newton/Raphson failed; switchover to bisection"));

      /* First bracket the root */
      left  = 0.;	                 	/* for sure */
      right = eslCONST_PI / sqrt(6.*variance);  /* an initial guess */
      lawless416(x, n, lambda, &fx, &dfx);
      while (fx > 0.) 
	{		
	  right *= 2.;		/* arbitrary leap to the right */
	  if (right > 100.) /* no reasonable lambda should be > 100, we assert */
	    ESL_EXCEPTION(eslENOHALT, "Failed to bracket root in esl_gumbel_FitComplete().");
	  lawless416(x, n, right, &fx, &dfx);
	}

      /* Now, bisection search in left/right interval */
      for (i = 0; i < 100; i++)
	{
	  mid = (left + right) / 2.; 
	  lawless416(x, n, mid, &fx, &dfx);
	  if (fabs(fx) < tol) break;             /* success */
	  if (fx > 0.)	left = mid;
	  else          right = mid;
	}
      if (i == 100) 
	ESL_EXCEPTION(eslENOHALT, "Even bisection search failed in esl_gumbel_FitComplete().");

      lambda = mid;
    }

  /* 3. Substitute into Lawless 4.1.5 to find mu
   */
  esum = 0.;
  for (i = 0; i < n; i++)
    esum  += exp(-lambda * x[i]);
  mu = -log(esum / n) / lambda;

  *ret_lambda = lambda;
  *ret_mu     = mu;   
  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id: esl_gumbel.c 727 2011-10-24 17:17:32Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/trunk/esl_gumbel.c $
 *****************************************************************/
