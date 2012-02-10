/*
 * NNVM.h
 *
 *  Created on: Feb 7, 2012
 *      Author: stephan
 */


/**
 * include guard and includes
 */

#ifndef NNVM_H_
#define NNVM_H_
#include "NNVMTree.h"
namespace NNVM {

template <class T> class NNVM {
public:
	 NNVM(NNVMTree) {};
	~NNVM() {};
};
}
#endif
