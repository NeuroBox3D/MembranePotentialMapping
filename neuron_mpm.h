/*!
 * \file neuron_mpm.h
 * \brief maps membrane potentials with NEURON to grid points
 *
 *  Created on: Jun 24, 2015
 *      Author: stephan
 */

/// guard
#ifndef  __H__UG__MEMBRANE_POTENTIAL_MAPPING__NEURON_MPM__
#define  __H__UG__MEMBRANE_POTENTIAL_MAPPING__NEURON_MPM__

/// includes
#include "vm2ug_rework.h"
#include "transformator.h"

/*! \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace membrane_potential_mapping {
		/*!
		 * \brief
		 */
		class NeuronMPM {
			private:
				/// private non-static members
				SmartPtr<Mapper<3, number> > m_mapper;
				SmartPtr<Transformator> m_transformator;

			public:
				/*!
				 * \brief dtor
				 */
				NeuronMPM() : m_mapper(make_sp(new Mapper<3, number>())),
							  m_transformator(SPNULL) {
				}

				/*!
				 * \brief set's the NEURON interpreter
				 */
				inline void set_transformator(SmartPtr<Transformator> transformator) {
					m_transformator = transformator;
				}

				/*!
				 * \brief get's the Mapper
				 */
				inline SmartPtr<Mapper<3, number> > get_mapper() const {
					return m_mapper;
				}

				/*!
				 * \brief builds the kd-tree on demand
				 */
				inline void build_tree() {
					m_mapper.get()->build_tree(m_transformator.get()->get_vms().front());
				}

				/*!
				 * \brief get membrane potential for given coordinates
				 */
				inline number get_vm(number x, number y, number z) {
					return m_mapper.get()->get_data_from_nearest_neighbor(MathVector<3>(x, y, z));
				}
		};
	} /// end namespace mpm
} /// end namespace ug
//<! \}

#endif ///  __H__UG__MEMBRANE_POTENTIAL_MAPPING__NEURON_MPM__
