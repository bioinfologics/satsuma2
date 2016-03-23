/* Foundation for the statistics modules.
 * 
 * Contents:
 *   1. The stats API.
 *   2. Unit tests.
 *   3. Test driver.
 *   4. Example.
 *   5. License and copyright information.
 * 
 */
#include <math.h>

#include "extern/easel/easel.h"
#include "extern/easel/esl_stats.h"


/* Function:  esl_stats_DMean()
 * Synopsis:  Calculates mean and $\sigma^2$ for samples $x_i$.
 * Incept:    SRE, Tue Jul 19 11:04:00 2005 [St. Louis]
 *
 * Purpose:   Calculates the sample mean and $s^2$, the unbiased
 *            estimator of the population variance, for a
 *            sample of <n> numbers <x[0]..x[n-1]>, and optionally
 *            returns either or both through <ret_mean> and
 *            <ret_var>.
 *            
 *            <esl_stats_FMean()> and <esl_stats_IMean()> do the same,
 *            for float and integer vectors.
 *
 * Args:      x        - samples x[0]..x[n-1]
 *            n        - number of samples
 *            opt_mean - optRETURN: mean
 *            opt_var  - optRETURN: estimate of population variance       
 *
 * Returns:   <eslOK> on success.
 */
int
esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var)
{
  double sum   = 0.;
  double sqsum = 0.;
  int i;

  for (i = 0; i < n; i++) 
    { 
      sum   += x[i];
      sqsum += x[i]*x[i];
    }
  if (opt_mean != NULL)  *opt_mean = sum / (double) n;
  if (opt_var  != NULL)  *opt_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
  return eslOK;
}
int
esl_stats_FMean(const float *x, int n, double *opt_mean, double *opt_var)
{
  double sum   = 0.;
  double sqsum = 0.;
  int i;

  for (i = 0; i < n; i++) 
    { 
      sum   += x[i];
      sqsum += x[i]*x[i];
    }
  if (opt_mean != NULL)  *opt_mean = sum / (double) n;
  if (opt_var  != NULL)  *opt_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
  return eslOK;
}
int
esl_stats_IMean(const int *x, int n, double *opt_mean, double *opt_var)
{
  double sum   = 0.;
  double sqsum = 0.;
  int i;

  for (i = 0; i < n; i++) 
    { 
      sum   += x[i];
      sqsum += x[i]*x[i];
    }
  if (opt_mean != NULL)  *opt_mean = sum / (double) n;
  if (opt_var  != NULL)  *opt_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id: esl_stats.c 685 2011-05-23 14:27:52Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/trunk/esl_stats.c $
 *****************************************************************/
