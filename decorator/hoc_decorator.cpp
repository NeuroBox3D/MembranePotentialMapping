/*
 * hoc_decorator.cpp
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#include "hoc_command.h"

using namespace ug::membrane_potential_mapping::decorator;

bool BaseDecorator::perform() {
	return m_baseCommand->perform();
}


