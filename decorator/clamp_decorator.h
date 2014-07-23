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
			struct ClampDecorator : BaseDecorator {
			public:
				/*!
				 * \brief default ctor
				 */
				ClampDecorator(std::unique_ptr<BaseCommand> cmd);

				bool perform();

			private:
				std::unique_ptr<BaseCommand> m_baseCommand;
				std::unique_ptr<ClampDecoratorImpl> m_impl;

			};

			/// TODO: does this work? would be handy if it would for VRL
			class ClampDecoratedCommand : ClampDecorator {
			private:
				const std::unique_ptr<BaseCommand> cmd = new SerialCommand();
			public:
				ClampDecoratedCommand() {
					this->ClampDecorator(cmd);
				}

				bool perform() {
					return this->ClampDecorator::perform();
				}
			};
		}
	}
}

#endif /// CLAMP_DECORATOR_H_
