//
//  Block.h
//  
//
//  Created by Federico on 24/11/17.
//

#ifndef BLOCK_H
#define BLOCK_H

#include "TObject.h"

typedef unsigned int uint;

struct BLock : public TObject{
  
  double T;
  double M;
  double S;
  int E;
  uint I;
  
  Block() : TObject() { }
  
  Block(const uint& Id, const double& _T, const int& _E, const double& _M, const double& _S, const uint& _I): TObject(), T(_T), E(_E), M(_M), S(_S), I(_I){
    SetUniqueID(Id);
  }
  
  
  
  
};

#endif /* Block_h */
