/*!
 * clamp_decorator.h
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#ifndef CLAMP_DECORATOR_H_
#define CLAMP_DECORATOR_H_

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			class ClampDecorator : BaseDecorator {
				/*!
				 * \brief default ctor
				 */
				ClampDecorator();

				bool perform();
			};
		}
	}
}

#endif /// CLAMP_DECORATOR_H_
