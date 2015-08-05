#ifndef _PROBTABLE_H_
#define _PROBTABLE_H_
#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <math.h>
#include "analysis/DNAVector.h"

// For testing 
#include "base/RandomStuff.h"


class ProbTable
{
  public:
    ProbTable();
    ProbTable(double targetSize, double cutoff);
    void Setup(double targetSize, double cutoff);
    bool IsGood(int length, double ident, double ident_expect);
    double GetMatchProbabilityRaw(int length, double ident, double ident_expect, double targetSize);
    double GetMatchProbability(double & ident,
        const DNAVector & target, 
        const DNAVector & query, 
        int startTarget,
        int startQuery,
        int length);
  private:
    int ExpectToIndex(double ident_expect);

    std::vector< std::vector < double > > m_table;

    double CDF(double x, double m, double s);
    double Sigma(double p, int N);
    double GCAdjustExpect(double gc, int N, double gc_target);

    double m_size;
    double m_cutoff;
    double dna_id_table[128][128];
    double dna_gc_table[128];

};

#endif
