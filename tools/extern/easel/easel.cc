/* Easel's foundation.
 * 
 * Contents:
 *    1. Exception and fatal error handling.
 *    2. Memory allocation/deallocation conventions.
 *    3. Standard banner for Easel miniapplications.
 *    4. Improved replacements for some C library functions.
 *    5. Portable drop-in replacements for nonstandard C functions.
 *    6. Additional string functions, esl_str*()
 *    7. File path/name manipulation, including tmpfiles.
 *    8. Typed comparison functions.
 *    9. Commonly used background composition (iid) frequencies.
 *   10. Unit tests.
 *   11. Test driver.
 *   12. Examples. 
 *   13. Copyright and license. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <unistd.h>		
#include <sys/stat.h>
#include <sys/types.h>

#include "extern/easel/easel.h"


/*****************************************************************
 * 1. Exception and fatal error handling.
 *****************************************************************/
static esl_exception_handler_f esl_exception_handler = NULL;

/* Function:  esl_exception()
 * Synopsis:  Throw an exception.
 *
 * Purpose:   Throw an exception. An "exception" is defined by Easel
 *            as an internal error that shouldn't happen and/or is 
 *            outside the user's control; as opposed to "failures", that       
 *            are to be expected, and within user control, and
 *            therefore normal. By default, exceptions are fatal.
 *            A program that wishes to be more robust can register
 *            a non-fatal exception handler.
 *            
 *            Easel programs normally call one of the exception-handling
 *            wrappers <ESL_EXCEPTION()> or <ESL_XEXCEPTION()>, rather
 *            than calling  <esl_exception> directly.
 *            
 *            If no custom exception handler has been registered, the
 *            default behavior is to print a brief message to <stderr>
 *            then <abort()>, resulting in a nonzero exit code from the
 *            program.  Depending on what <errcode>, <sourcefile>,
 *            <sourceline>, and the <sprintf()>-formatted <format>
 *            are, this output looks like:
 *            
 *            Fatal exception (source file foo.c, line 42):
 *            Something wicked this way came.
 *            
 * Args:      errcode     - Easel error code, such as eslEINVAL. See easel.h.
 *            use_errno   - if TRUE, also use perror() to report POSIX errno message.
 *            sourcefile  - Name of offending source file; normally __FILE__.
 *            sourceline  - Name of offending source line; normally __LINE__.
 *            format      - <sprintf()> formatted exception message, followed
 *                          by any additional necessary arguments for that 
 *                          message.
 *                          
 * Returns:   void. 
 *
 * Throws:    No abnormal error conditions. (Who watches the watchers?)
 */
void
esl_exception(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, ...)
{
  va_list argp;

  if (esl_exception_handler != NULL) 
    {
      va_start(argp, format);
      (*esl_exception_handler)(errcode, use_errno, sourcefile, sourceline, format, argp);
      va_end(argp);
      return;
    } 
  else 
    {
      fprintf(stderr, "Fatal exception (source file %s, line %d):\n", sourcefile, sourceline);
      va_start(argp, format);
      vfprintf(stderr, format, argp);
      va_end(argp);
      fprintf(stderr, "\n");
      if (use_errno && errno) perror("system error");
      fflush(stderr);
      abort();
    }
}

/* Function:  esl_exception_SetHandler()
 * Synopsis:  Register a different exception handling function.
 *
 * Purpose:   Register a different exception handling function,
 *            <handler>. When an exception occurs, the handler
 *            receives at least four arguments: <errcode>, <sourcefile>,
 *            <sourceline>, and <format>. 
 * 
 *            <errcode> is an Easel error code, such as
 *            <eslEINVAL>. See <easel.h> for a list of all codes.
 * 
 *            <use_errno> is TRUE for POSIX system call failures. The
 *            handler may then use POSIX <errno> to format/print an
 *            additional message, using <perror()> or <strerror_r()>.
 *           
 *            <sourcefile> is the name of the Easel source code file
 *            in which the exception occurred, and <sourceline> is 
 *            the line number.
 *            
 *            <format> is a <vprintf()>-formatted string, followed by
 *            a <va_list> containing any additional arguments that
 *            formatted message needs.  Your custom exception handler
 *            will probably use <vfprintf()> or <vsnprintf()> to format
 *            its error message.
 *            
 * Args:      handler -  ptr to your custom exception handler.
 *
 * Returns:   void.
 *
 * Throws:    (no abnormal error conditions)
 */
void
esl_exception_SetHandler(void (*handler)(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, va_list argp))
{ 
  esl_exception_handler = handler; 
}


/* Function:  esl_exception_ResetDefaultHandler()
 * Synopsis:  Restore default exception handling.
 *
 * Purpose:   Restore default exception handling, which is to print
 *            a simple error message to <stderr> then <abort()> (see
 *            <esl_exception()>. 
 *      
 *            An example where this might be useful is in a program
 *            that only temporarily wants to catch one or more types
 *            of normally fatal exceptions.
 *            
 *            If the default handler is already in effect, this 
 *            call has no effect (is a no-op).
 *
 * Args:      (void)
 *
 * Returns:   (void)
 *
 * Throws:    (no abnormal error conditions)
 */
void
esl_exception_ResetDefaultHandler(void)
{
  esl_exception_handler = NULL; 
}


/* Function: esl_nonfatal_handler()
 * Synopsis: A trivial example of a nonfatal exception handler.
 * 
 * Purpose:  This serves two purposes. First, it is the simplest
 *           example of a nondefault exception handler. Second, this
 *           is used in test harnesses, when they have
 *           <eslTEST_THROWING> turned on to test that thrown errors
 *           are handled properly when a nonfatal error handler is
 *           registered by the application.
 *           
 * Args:      errcode     - Easel error code, such as eslEINVAL. See easel.h.
 *            use_errno   - TRUE on POSIX system call failures; use <errno> 
 *            sourcefile  - Name of offending source file; normally __FILE__.
 *            sourceline  - Name of offending source line; normally __LINE__.
 *            format      - <sprintf()> formatted exception message.
 *            argp        - <va_list> containing any additional necessary arguments for 
 *                          the <format> message.
 *                          
 * Returns:   void. 
 *
 * Throws:    (no abnormal error conditions)
 */
void
esl_nonfatal_handler(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, va_list argp)
{
  return; 
}


/* Function:  esl_fatal()
 * Synopsis:  Kill a program immediately, for a "violation".
 *
 * Purpose:   Kill a program for a "violation". In general this should only be used
 *            in development or testing code, not in production
 *            code. The main use of <esl_fatal()> is in unit tests.
 *            Another use is in assertions used in dev code.
 *            
 *            The only other case (and the only case that should be allowed in
 *            production code) is in a true "function" (a function that returns
 *            its answer, rather than an Easel error code), where Easel error
 *            conventions can't be used (because it can't return an error code),
 *            AND the error is guaranteed to be a coding error. For an example,
 *            see <esl_opt_IsOn()>, which triggers a violation if the code
 *            checks for an option that isn't in the code.
 * 
 * Args:      format  - <sprintf()> formatted exception message, followed
 *                      by any additional necessary arguments for that 
 *                      message.
 *
 * Returns:   (void)
 *
 * Throws:    (no abnormal error conditions)
 */
void
esl_fatal(const char *format, ...)
{
  va_list argp;

  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
  exit(1);
}
/*---------------- end, error handling conventions --------------*/




/*****************************************************************
 * 2. Memory allocation/deallocation conventions.
 *****************************************************************/

/* Function:  esl_Free2D()
 *
 * Purpose:   Free a 2D pointer array <p>, where first dimension is
 *            <dim1>. (That is, the array is <p[0..dim1-1][]>.)
 *            Tolerates any of the pointers being NULL, to allow
 *            sparse arrays.
 *
 * Returns:   void.
 */
void
esl_Free2D(void **p, int dim1)
{
  int i;
  if (p != NULL) {
    for (i = 0; i < dim1; i++)
      if (p[i] != NULL) free(p[i]);
    free(p);
  }
  return;
}

/* Function:  esl_Free3D()
 *
 * Purpose:   Free a 3D pointer array <p>, where first and second
 *            dimensions are <dim1>,<dim2>. (That is, the array is
 *            <p[0..dim1-1][0..dim2-1][]>.) Tolerates any of the
 *            pointers being NULL, to allow sparse arrays.
 *
 * Returns:   void.
 */
void
esl_Free3D(void ***p, int dim1, int dim2)
{
  int i, j;

  if (p != NULL) {
    for (i = 0; i < dim1; i++)
      if (p[i] != NULL) {
        for (j = 0; j < dim2; j++)
          if (p[i][j] != NULL) free(p[i][j]);
        free(p[i]);
      }
    free(p);
  }
}
/*------------- end, memory allocation conventions --------------*/

/*****************************************************************
 * 8. Typed comparison routines.
 *****************************************************************/

/* Function:  esl_DCompare()
 *
 * Purpose:   Compare two floating point scalars <a> and <b> for approximate equality.
 *            Return <eslOK> if equal, <eslFAIL> if not.
 *            
 *            Equality is defined by being within a relative
 *            epsilon <tol>, as <2*fabs(a-b)/(a+b)> $\leq$ <tol>.
 *            Additionally, we catch the special cases where <a>
 *            and/or <b> are 0 or -0. If both are, return <eslOK>; if
 *            one is, check that the absolute value of the other is
 *            $\leq$ <tol>.
 *            
 *            <esl_DCompare()> and <esl_FCompare()> work on <double> and <float>
 *            scalars, respectively.
 */
int
esl_DCompare(double a, double b, double tol)
{
  if (isinf(a) && isinf(b))                 return eslOK;
  if (isnan(a) && isnan(b))                 return eslOK;
  if (!isfinite(a) || !isfinite(b))         return eslFAIL;
  if (a == b)                               return eslOK;
  if (fabs(a) == 0. && fabs(b) <= tol)      return eslOK;
  if (fabs(b) == 0. && fabs(a) <= tol)      return eslOK;
  if (2.*fabs(a-b) / fabs(a+b) <= tol)      return eslOK;
  return eslFAIL;
}
int
esl_FCompare(float a, float b, float tol)
{ 
  if (isinf(a) && isinf(b))                 return eslOK;
  if (isnan(a) && isnan(b))                 return eslOK;
  if (!isfinite(a) || !isfinite(b))         return eslFAIL;
  if (a == b)                               return eslOK;
  if (fabs(a) == 0. && fabs(b) <= tol)      return eslOK;
  if (fabs(b) == 0. && fabs(a) <= tol)      return eslOK;
  if (2.*fabs(a-b) / fabs(a+b) <= tol)      return eslOK;
  return eslFAIL;
}

/* Function:  esl_DCompareAbs()
 *
 * Purpose:   Compare two floating point scalars <a> and <b> for
 *            approximate equality, by absolute difference.  Return
 *            <eslOK> if equal, <eslFAIL> if not.
 *            
 *            Equality is defined as <fabs(a-b) \leq tol> for finite
 *            <a,b>; or <inf=inf>, <NaN=NaN> when either value is not
 *            finite.
 *            
 *            Generally it is preferable to compare floating point
 *            numbers for equality using relative difference: see
 *            <esl_{DF}Compare()>, and also Knuth's Seminumerical
 *            Algorithms. However, cases arise where absolute
 *            difference comparison is preferred. One such case is in
 *            comparing the log probability values of DP matrices,
 *            where numerical error tends to accumulate on an absolute
 *            scale, dependent more on the number of terms than on
 *            their magnitudes. DP cells with values that happen to be
 *            very close to zero can have high relative differences.
 */
int
esl_DCompareAbs(double a, double b, double tol)
{
  if (isinf(a) && isinf(b))            return eslOK;
  if (isnan(a) && isnan(b))            return eslOK;
  if (!isfinite(a) || !isfinite(b))    return eslFAIL;
  if (fabs(a-b) <= tol)                return eslOK;
  return eslFAIL;
}
int
esl_FCompareAbs(float a, float b, float tol)
{ 
  if (isinf(a) && isinf(b))            return eslOK;
  if (isnan(a) && isnan(b))            return eslOK;
  if (!isfinite(a) || !isfinite(b))    return eslFAIL;
  if (fabs(a-b) <= tol)                return eslOK;
  return eslFAIL;
}





/* Function:  esl_CCompare()
 * Synopsis:  Compare two optional strings for equality.
 *
 * Purpose:   Compare two optional strings <s1> and <s2>
 *            for equality. 
 *            
 *            If they're non-<NULL> and identical up to their
 *            <NUL>-terminator, return <eslOK>.
 *            
 *            If they're both <NULL> (unset), return <eslOK>.
 *            
 *            Otherwise, they're not identical; return <eslFAIL>.
 */
int
esl_CCompare(char *s1, char *s2)
{
  if (s1 == NULL && s2 == NULL) return eslOK;
  if (s1 == NULL || s2 == NULL) return eslFAIL;
  if (strcmp(s1, s2) != 0)      return eslFAIL;
  return eslOK;
}


/*-------------- end, typed comparison routines --------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id: easel.c 748 2012-02-14 21:23:06Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/trunk/easel.c $
 *****************************************************************/  
