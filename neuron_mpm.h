/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2015-06-24
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/// guard
#ifndef  UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__NEURON_MPM_H
#define  UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__NEURON_MPM_H

/// includes
#include "transformator.h"
#include "vm2ug.h"

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
		 * \brief get's the NEURON interpreter
		 */
		inline SmartPtr<Transformator> get_transformator() const {
			return m_transformator;
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

} // end namespace mpm
} // end namespace ug
//<! \}

#endif // UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__NEURON_MPM_H
