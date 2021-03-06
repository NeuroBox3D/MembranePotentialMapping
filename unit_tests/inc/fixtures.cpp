/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2012-05-27
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

// includes
#include <boost/test/included/unit_test.hpp>

#include "mpm_config.h"  // for project-specific defines

#if (defined(MPMVGCC) && MPMVGCC == 1)
#include "../../bg_simple/bg.h"
#elif (defined(MPMVGCC) && MPMVGCC == 2)
#include "../../bg_default/bg.h"
#endif


// using directives
using namespace boost::unit_test;


/**
 * \brief fixture for class Vm2uG
template <class T> struct FixtureVUG {
	ug::membrane_potential_mapping::Vm2uG<T>* vm2ug;
	FixtureVUG() : vm2ug(new ug::membrane_potential_mapping::Vm2uG<T>("","")) { BOOST_TEST_MESSAGE("setup fixture >>vm2ugG<<");}
	~FixtureVUG() { BOOST_TEST_MESSAGE("teardown fixture >>vm2ug<<"); }

};
*/


#if (defined(MPMVGCC) && (MPMVGCC == 1 || MPMVGCC == 2))
/**
 * \brief fixture for class BG
 */
struct FixtureBG {
	ug::membrane_potential_mapping::bg::BG* bg;
	static size_t count;
	/**
	 * \brief construct a fixture for class bg
	 */
	FixtureBG() : bg(new ug::membrane_potential_mapping::bg::BG()) {  BOOST_TEST_MESSAGE("setup fixture >>bg<< no. " << ++count);}
	~FixtureBG() { BOOST_TEST_MESSAGE("teardown fixture >>bg<<"); }
};

// initialization of static members
size_t FixtureBG::count = 0;
#endif
