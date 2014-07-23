/*
 * clamp_decorator_impl.h
 *
 *  Created on: Jul 11, 2014
 *      Author: stephan
 */

#ifndef CLAMP_DECORATOR_IMPL_H_
#define CLAMP_DECORATOR_IMPL_H_

namespace ug {
	namespace membrane_potential_mapping {
		namespace decorator {
			struct ClampDecoratorImpl {
			public:
				bool perform();
			};
		}
	}
}



#endif /* CLAMP_DECORATOR_IMPL_H_ */
