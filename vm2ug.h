/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2015-06-19
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

#ifndef UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__VM2UG_H
#define UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__VM2UG_H

#include <cstddef>
#include <vector>
#include <utility>
#include "common/types.h"  // for number
#include "common/math/math_vector_matrix/math_vector.h"  // for MathVector
#include "kdtree/kd_tree.h" // for kdtree

/*! \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
  namespace membrane_potential_mapping {

  template <size_t dim, typename M> class Mapper {
  private:
    //// non-static members
    kd_tree<dim, M> m_kdtree;
    size_t m_delimIndex;

    /// common typedefs
    typedef typename std::vector<std::pair<MathVector<dim, number>, M> >::const_iterator CITVPMNM;
    typedef typename std::vector<std::pair<std::vector<number>, M> >::const_iterator CITVPVNM;

  public:
    /*!
    * \brief ctor
    */
    Mapper()
    : m_kdtree(kd_tree<dim, M>()),
      m_delimIndex(0) { }

   /*!
    * \brief build an empty tree
    */
   void build_tree();

   /*!
    * \brief build populated tree from given vector of pairs
    * \param[in] points to be used for the tree construction
    */
   void build_tree
   (
      const std::vector<std::pair<MathVector<dim, number>, M> >& points
   );

   /*!
    * \brief build populated tree from given file
    * \param[in] filename where the points are stored (with meta data)
    * \param[in] delimiter in the file to separate values
    */
   void build_tree
   (
     const std::string& filename
   );

    /*!
    * \brief build populated tree from given vector of pairs
    * \param[in] points to be used for the tree construction
    */
   void build_tree
   (
     const std::vector<std::pair<std::vector<number>, M> >& points
   );
 
   /*!
    * \brief add a single node
    * \param[in] node a single point with meta data
    */
   void add_node
   (
     const std::pair<MathVector<dim, number>, M>& node
   );

   /*!
    * \brief add a single node
    * \param[in] node a single point with meta data
    */
   void add_node
   (
     const std::pair<std::vector<number>, M>& node
   );

  /*!
   * \brief add a single node with a value
   * \param[in] node a single point with meta data
   * \param[in] value the meta data
   */
  inline void add_node
  (
    const std::vector<number>&, 
    const M& value
  );

   /*!
    * \brief query the tree for the data of the very nearest neighbor
    * \param[in] query coordinates of a given point
    * \return \c data
    */
   M get_data_from_nearest_neighbor
   (
     const MathVector<dim, number>& query
   ) const;
    

   /*!
    * \brief query the tree for the data of the very nearest neighbor
    * \param[in] query coordinates of a given point
    * \return \c data
    */
   M get_data_from_nearest_neighbor
   (
     const std::vector<number>& query
   ) const;
   };
  } // end namespace mpm
} // end namespace ug
//<! \}

#include "vm2ug_impl.h"
#endif // UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__VM2UG_H
