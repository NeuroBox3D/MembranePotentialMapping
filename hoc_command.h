/*!
 * \file hoc_command.h
 * \brief execute a hoc command
 *
 *  Created on: May 30, 2014
 *      Author: stephan
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMAND_BUILDER__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMAND_BUILDER__

// extern C functions (NEURON)
extern bool hoc_valid_stmt(const char* stmt, Object* ob);
extern int ivocmain(int, char**, char**);
extern double hoc_ac_;

#include "common/log.h"
#include "common/error.h"

namespace ug {
	namespace membrane_potential_mapping {
		/*!
		 * \class HocCommand
		 */
		template <typename T>
		class HocCommand {
		private:
			// number of arguments to hoc interpreter
			static int ARGC;
			// arguments to hoc interpreter
			static char* ARGV[];
			// environment of the hoc interpreter
			static char* ENV[];
			// command to be build
			std::stringstream m_command;

		public:
			HocCommand(const std::vector<T>& commands) {
				for (std::vector<T>::const_iterator it = commands.begin(); it != commands.end(); ++it) {
					m_command << *it;
				}
			}

			HocCommand(const T& val) {
				m_command << val;
			}

			HocCommand() {
			}

			SmartPtr<HocCommand> build(const T& command) {
				try {
					m_command << command;
				} UG_CATCH_THROW("HocCommand.build() failed" << std::endl);
				return SmartPtr<this>;
			}

			bool execute() const {
				return hoc_valid_stmt(m_command.str().c_str(), NULL);
			}

			number result() const {
				return hoc_ac_;
			}

			void clear() {
				m_command.clear();
			}
		};
	}
}

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMAND_BUILDER__
