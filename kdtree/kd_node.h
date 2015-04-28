/*!
 * \file kd_node.h
 *
 *  Created on: Apr 28, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_NODE__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_NODE__

/* \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace membrane_potential_mapping {
		/*!
		 * \brief kd node with meta data
		 */
		template <size_t dim, typename M>
        struct kd_node {
           M m_coords[dim];
           M m_meta;
           kd_node<dim, M>* left;
           kd_node<dim, M>* right;

           /*!
            * \brief pretty print point coordinates with meta
            * \function operator<<
            */
           friend std::ostream& operator<<(std::ostream &out, const kd_node<dim, M>& node) {
        	   for (size_t i = 0; i < dim-1; i++) out << node.m_coords[i] << ", ";
        	   return (out << node.m_coords[dim] << "(" << node.m_meta << ")");
           }
        };
	} // namespace membrane_potential_mapping
} // namespace ug
//<! \}

#endif ///__H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_NODE__
