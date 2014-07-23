/*!
 * parallel_hoc_command.cpp
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#include "parallel_hoc_command.h"

using namespace ug::membrane_potential_mapping::decorator;

/////////////////////////////////////////////////////////
/// perform
/////////////////////////////////////////////////////////
bool perform() {
	std::cout << "ParallelCommand invoked" << std::endl;
	return true;
}
