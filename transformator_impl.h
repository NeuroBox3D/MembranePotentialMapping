/*
 * transformator_impl.h
 *
 *  Created on: Mar 17, 2014
 *      Author: stephan
 */

#ifndef TRANSFORMATOR_IMPL_H_
#define TRANSFORMATOR_IMPL_H_

#include "transformator.h"

#include <lib_grid/common_attachments.h>
#include <lib_grid/grid/grid_util.h>

#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug {
	namespace membrane_potential_mapping {
	    /////////////////////////////////////////////////////////
		/// set_grid
	    /////////////////////////////////////////////////////////
		template <typename TGridFunction, typename TDomain>
		void Transformator::set_grid(TGridFunction& u, const char* uCmp) {
			SmartPtr<Domain<3> > spDom = u.domain();
			SmartPtr<ug::Grid> spGrid = spDom.get()->grid();
			mGrid = spGrid.get();
		}

	    /////////////////////////////////////////////////////////
		/// integrate_feedback
	    /////////////////////////////////////////////////////////
		template <typename TGridFunction, typename TDomain>
		number Transformator::integrate_feedback(TGridFunction& u, const char*cmp, const char* subsets, int quadOrder)  {
			m_saved_neuron_point_integral = Integral(u, cmp, subsets, quadOrder);
			return m_saved_neuron_point_integral;
		}

	    /////////////////////////////////////////////////////////
		/// integrate_feedbacks
	    /////////////////////////////////////////////////////////
		template <typename TGridFunction, typename TDomain>
		void Transformator::integrate_feedbacks(TGridFunction& u, const char*cmp, const char* subsets, int quadOrder)  {
			mIntegrals .push_back( Integral(u, cmp, subsets, quadOrder));
		}



	}
}

#endif /* TRANSFORMATOR_IMPL_H_ */
