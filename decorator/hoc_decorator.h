/*!
 * hoc_decorator.h
 * TODO: PIMPL classes and decorators can be refactored to structs!
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */
#ifndef HOC_DECORATOR_H
#define HOC_DECORATOR_H

#include "hoc_command.h"

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			/*!
			 * \brief base decorator
			 */
			struct BaseDecorator : BaseCommand {
			private:
				std::unique_ptr<BaseCommand> m_baseCommand;

			public:
				/*!
				 * \brief default ctor
				 */
				BaseDecorator(std::unique_ptr<BaseCommand> cmd);

				bool perform();
			};
		}
	}
}

#endif /// HOC_DECORATOR_H_



