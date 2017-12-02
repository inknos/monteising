#ifndef INFO_H
#define INFO_H
typedef unsigned int uint;

typedef struct{
  double _tempmin;
  double _tempmax;
  uint _N;
  uint _dim;
  uint _dim_vector;
  uint _I0;
  uint _iter;
  uint _tempstep;
  uint _datime_t;
} INFO;

#endif
