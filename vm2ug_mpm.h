/*!
 * \file vm2ug_mapper.h
 * \brief maps membrane potentials without NEURON to grid points
 *
 *  Created on: Jun 24, 2015
 *      Author: stephan
 */

/// guard
#ifndef  __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_MPM__
#define  __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_MPM__

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
		class Vm2uGMPM {
			/// private non-static members
			SmartPtr<Mapper<3, number> > m_mapper;

			public:
			/*!
			 * \brief dtor
			 */
			Vm2uGMPM() : m_mapper(make_sp(new Mapper<3, number>())) {

			}

			/*!
			 * \brief builds the kd-tree on demand with a file
			 * \param[in] filename
			 * \param[in] delimiter
			 */
			inline void build_tree(const std::string& filename, const std::string& delimiter) {
				m_mapper->build_tree(filename, delimiter);
			}

			/*!
			 * \brief builds the kd-tree on demand without a file
			 * \param[in] points
			 */
			inline void build_tree(const std::vector<std::pair<MathVector<3>, number> >& points) {
				m_mapper->build_tree(points);
			}

			/*!
			 * \brief get's membrane potential for given coordinates x, y and z
			 * \param[in] x
			 * \param[in] y
			 * \param[in] z
			 *
			 * \return \c membrane potential
			 */
			inline number get_vm(number x, number y, number z) {
				return get_vm(MathVector<3>(x, y, z));
			}

			/*!
			 * \brief get membrane potential for given vector
			 * \param[in] vec
			 *
			 * \return \c membrane potential
			 */
			inline number get_vm(const MathVector<3>& vec) {
				return m_mapper->get_data_from_nearest_neighbor(vec);
			}
		};
	} /// end namespace mpm
} /// end namespace ug
//<! \}

#endif ///  __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_MPM__
