/*!
 * plain_decorator.h
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#ifndef PLAIN_DECORATOR_H_
#define PLAIN_DECORATOR_H_

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			class PlainDecorator : BaseDecorator {
			public:
				/*!
				 * \brief default ctor
				 */
				PlainDecorator(std::unique_ptr<BaseCommand> cmd);

			};
		}
	}
}

#endif /// PLAIN_DECORATOR_H
