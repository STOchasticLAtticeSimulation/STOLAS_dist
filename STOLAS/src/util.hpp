#ifndef INCLUDED_util_hpp_
#define INCLUDED_util_hpp_

// useful macro
#define LOOP for(int i = 0; i < NLnoise; i++) for(int j = 0; j < NLnoise; j++) for(int k = 0; k < NLnoise; k++)

inline int index(int i, int j, int k) {
  return i * NLnoise * NLnoise + j * NLnoise + k;
}

#endif