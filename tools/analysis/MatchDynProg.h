#ifndef MATCHDYNPROG_H_
#define MATCHDYNPROG_H_

#include "SequenceMatch.h"


bool RunMatchDynProg(MultiMatches & out, const MultiMatches & in, bool bSelf = false);
bool RunMatchDynProgMult(MultiMatches & out, const MultiMatches & in);

 


#endif //MATCHDYNPROG_H_


