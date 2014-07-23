/*!
 * parallel_hoc_command_impl.cpp
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#include "parallel_hoc_command_impl.h"

using namespace ug::membrane_potential_mapping::decorator;

bool ParallelCommand::perform() {
	std::cout << "Parallel Command invoked" << std::endl;
	return true;
}


