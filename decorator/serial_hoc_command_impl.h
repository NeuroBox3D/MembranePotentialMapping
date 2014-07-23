/*!
 * serial_hoc_command_impl.h
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#ifndef SERIAL_HOC_COMMAND_IMPL_H_
#define SERIAL_HOC_COMMAND_IMPL_H_

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			struct SerialCommandImpl {
			/// private data structure to perform the command
			private:

			public:
				SerialCommandImpl();

				bool perform();
			};
		}
	}
}



#endif /* SERIAL_HOC_COMMAND_IMPL_H_ */
