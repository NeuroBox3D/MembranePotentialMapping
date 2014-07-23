/*!
 * parallel_hoc_command_impl.h
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#ifndef PARALLEL_HOC_COMMAND_IMPL_H_
#define PARALLEL_HOC_COMMAND_IMPL_H_

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			struct ParallelCommandImpl {
			/// private data structure to perform the command
			private:

			public:
				ParallelCommandImpl();

				bool perform();
			};
		}
	}
}

#endif /* PARALLEL_HOC_COMMAND_IMPL_H_ */
