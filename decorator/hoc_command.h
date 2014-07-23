/*!
 * hoc_command.h
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#ifndef HOC_COMMAND_H_
#define HOC_COMMAND_H_

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			/*!
			 * \brief base hoc command
			 */
			class BaseCommand {
			/// public methods
			public:
				/*!
				 * \brief virtual dtor
				 */
				virtual ~BaseCommand() = 0;

				/*!
				 * \brief execute commands
				 */
				virtual bool perform() = 0;
			};
		}
	}
}

#endif /// HOC_COMMAND_H_
