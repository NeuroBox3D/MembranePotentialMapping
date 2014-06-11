/*!
 * \file hoc_command_impl.h
 *  Created on: May 30, 2014
 *      Author: stephan
 */

#ifndef __H__UG__PLASMA_MEMBRANE__HOC_COMMAND_IMPL__
#define __H__UG__PLASMA_MEMBRANE__HOC_COMMAND_IMPL__

#include "hoc_command.h"

namespace ug {
	namespace membrane_potential_mapping {
		template <typename T>
		HocCommand::HocCommand(const std::vector<T>& commands) {
			for (std::vector<T>::const_iterator it = commands.begin(); it != commands.end(); ++it) {
				m_command << *it;
			}
		}

		template <typename T>
		HocCommand::HocCommand(const T& val) {
			m_command << val;
		}

		template <typename T>
		void HocCommand::build(const T& command) {
			try {
				m_command << command;
			} UG_CATCH_THROW("HocCommand.build() failed" << std::endl);
		}
	}
}
#endif //  __H__UG__PLASMA_MEMBRANE__HOC_COMMAND_IMPL__
