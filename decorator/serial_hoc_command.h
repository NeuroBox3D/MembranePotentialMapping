/*!
 * serial_hoc_command.h
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */
#ifndef SERIAL_HOC_COMMAND_H
#define SERIAL_HOC_COMMAND_H

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			/*!
			 * \brief serial command class
			 */
			class SerialCommand : BaseCommand {
			private:
				std::unique_ptr<SerialCommandImpl> m_impl;
			public:
				/*!
				 * \brief default ctor
				 */
				SerialCommand();

				/*!
				 * \brief perform command
				 */
				bool perform() {
					return m_impl->perform();
				}
			};
		}
	}
}

#endif /// SERIAL_HOC_COMMAND_H
