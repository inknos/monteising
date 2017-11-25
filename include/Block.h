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

struct Block : public TObject{
  
  double T;
  double M;
  double S;
  int E;
  uint I;
  
  Block();
  
  Block(const uint& Id, const double& _T, const int& _E, const double& _M, const double& _S, const uint& _I);

ClassDef(Block, 1)
};
#endif /* Block_h */
