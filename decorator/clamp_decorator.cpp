/*
 * clamp_decorator.cpp
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#include "clamp_decorator.h"

using namespace ug::membrane_potential_mapping::decorator;

bool ClampDecorator::perform() {
	return m_impl->perform();
}
