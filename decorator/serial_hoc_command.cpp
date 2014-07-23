/*!
 * serial_hoc_command.cpp
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */
#include "serial_hoc_command.h"

using namespace ug::membrane_potential_mapping::decorator;

/////////////////////////////////////////////////////////
/// perform
/////////////////////////////////////////////////////////
bool SerialCommand::perform() {
	return m_impl->perform();
}



