//
//  Block.cpp
//  
//
//  Created by Federico on 25/11/17.
//

#include "Block.h"


ClassImp(Block)

Block::Block() : TObject() { }
Block::Block(const uint& Id, const double& _T, const int& _E, const double& _M, const double& _S, const uint& _I): TObject(), T(_T), E(_E), M(_M), S(_S), I(_I){
  SetUniqueID(Id);
}
