/*!
 * parallel_hoc_command.h
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#ifndef PARALLEL_HOC_COMMAND_H_
#define PARALLEL_HOC_COMMAND_H_

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			/*!
			 * \brief parallel command class
			 */
			class ParallelCommand : BaseCommand {
			public:
				/*!
				 * \brief default ctor
				 */
				ParallelCommand();

				/*!
				 * \brief perform command
				 */
				bool perform();

			};

		}
	}
}

#endif /// PARALLEL_HOC_COMMAND_H_
