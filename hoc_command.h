/*!
 * \file hoc_command.h
 * \brief execute a hoc command
 *
 *  Created on: May 30, 2014
 *      Author: stephan
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMAND_BUILDER__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMAND_BUILDER__


#include "common/log.h"
#include "common/error.h"

namespace ug {
	namespace membrane_potential_mapping {
		/*!
		 * \class HocCommand
		 */
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
			/*!
			 * \brief default ctor
			 */
			HocCommand() {
			}

			/*!
			 * \brief default dtor
			 */
			~HocCommand() {
			}

			/*!
			 * \brief ctor
			 * \param[in] val
			 */
			template <typename T>
			HocCommand(const T& val);

			/*!
			 * \brief ctor
			 * \param[in] commands
			 */
			template <typename T>
			HocCommand(const std::vector<T>& commands);

			/*!
			 * \brief build the comamnd string
			 * \param[in] command
			 */
			template <typename T>
			void build(const T& command);

			/*!
			 * \brief execute hoc statement
			 */
			bool execute() const;

			/*!
			 * \brief get result of hoc statement execution
			 */
			number result() const;

			/*!
			 * \brief clear the hoc command string
			 */
			void clear();
		};
	}
}

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMAND_BUILDER__
