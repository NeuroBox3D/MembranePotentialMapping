/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2013-11-06
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

#ifndef UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__TRANSFORMATOR_H
#define UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__TRANSFORMATOR_H

#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <iostream>

#include "mpm_config.h"  // for project-specific defines

#include "common/types.h"
#define nil NULL

#include <lib_grid/common_attachments.h>
#include <lib_grid/grid/grid_util.h>

#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug {
namespace membrane_potential_mapping {

class Transformator {
  /// variables
  private:
    // the integrals
    std::vector<number> mIntegrals;
    // the points
    std::vector<ug::vector3> mNeuronPoints;
    // saved point
    ug::vector3 m_saved_neuron_point;
    // saved integral
    number m_saved_neuron_point_integral;
    // a grid
    ug::Grid* mGrid;
    // statements accumulator
    std::vector<std::string> m_stmts;
    // number of sections in hoc file
    size_t m_sections;
    // number of total points
    size_t m_totalPoints;
    // finitialize value of points
    number m_finitialize;
    // tstart
    number m_tstart;
    // tstop
    number m_tstop;
    // dt
    number m_dt;
    // t
    number m_t;
    // number of arguments to hoc interpreter
    static int ARGC;
    // arguments to hoc interpreter
    static char* ARGV[];
    // environment of the hoc interpreter
    static char* ENV[];
    // holds membrane potentials for limit timesteps
    std::vector<std::vector<std::pair<std::vector<double>, double> > > m_vms;
    // timestep
    size_t m_limit;
    // holds defining hoc stimulation protocol file
    std::string m_stimulation;
    // holds defining hoc geometry file
    std::string m_geometry;


  /// functions
  public:
    /*!
     * \brief adjust intracellular resistivity
     *
     * Sets the intracellular resisitivity for an obstacle compartment
     *
     * \param[in] uCmp the grid function name
     * \param[in] subset the obstacle subset
     * \param[in] Rai the new intracellular resistivity
     */
    void adjust_resistivity_obstacle(const char* uCmp, const char* subset, double Rai);

    /*!
     * \brief feedback to NEURON from UG
     *
     * Feedback coupling to NEURON from UG
     *
     * \param[in] u the grid function
     * \param[in] uCmp the grid function name
     * \param[in] subset the input subset for calcium flux
     * \param[in] subset_vol the volume subset inside of the surrounding surface subset for calcium influx
     * \param[in] density the density to be set
     */
    void feedback(const char* uCmp, const char* subset, const char* subset_vol, double density);

    /*!
     * \brief feedback to NEURON from UG
     *
     * See above, does the same except it allows multiple compartment points per subset
     */
    void feedbacks(const char* uCmp, const char* subset, const char* subset_vol, double density);

    /*!
     * \brief inits the feedback, i. e. insert the channel mechanism for subset
     */
    void init_feedback(const char* uCmp, const char* subset, const char* subset_vol, double density, const char* channel);

    /*!
     * \brief inits the feedback, i.e. inserts the channel mechanism for subset
     */
    void init_feedbacks(const char* uCmp, const char* subset, const char* subset_vol, double density, const char* channel);

    /*!
     * \brief set the actual feedback
     *
     * Sets the feedback then, helper function is necessary, because using NEURON with templates is very cumbersome...
     */
    void set_feedback();

    /*!
     * \brief sets the actual feedback
     *
     * See above, does the same for many compartment points per subset
     */
    void set_feedbacks();

    /*!
     * \brief sets the grid
     *
     * \param[in] u the grid function
     * \param[in] uCmp the grid function name
     */
    template <typename TGridFunction, typename TDomain>
    void set_grid(TGridFunction&u, const char* uCmp);

    /*!
     * \brief  integrate the calcium concentration
     *
     * \param[in] u the grid functionm
     * \param[in] uCmp the grid function name
     * \param[in] subsets the subset to be considered
     * \param[in] quadOrder which quad order to be used
     */
    template <typename TGridFunction, typename TDomain>
    double integrate_feedback(TGridFunction& u, const char*cmp, const char* subsets, int quadOrder);


    /*!
     * \brief integrate the calcium concentration
     *
     * See: above, does the same except it works for more than one compartment per subset
     */
    template <typename TGridFunction, typename TDomain>
    void integrate_feedbacks(TGridFunction& u, const char*cmp, const char* subsets, int quadOrder);

    /*!
     * \brief default ctor
     *
     * Initializes the hoc interpreter with default values,
     * set's up the hoc environment for subsequent execution
     * of hoc statements with respect to a preparation for the
     * membrane potential mapping plugin.
     */
    Transformator();

    /*!
     * \brief enhanced ctor I
     *
     * \param[in] argc number of command line arguments
     * \param[in] argv commanad line arguments as "strings"
     * \param[in] env hoc environment as "strings"
     *
     * In principle the same as the default ctor, but you
     * may supply own command line arguments and an own
     * environment for the hoc interpreter. This should only
     * be necessary in some really rare cases.
     */
    Transformator(int argc, char* argv[], char* env[]);

    /*!
     * \brief enhanced ctor II - supply mod files as argument
     *
     * Example: supply "-dll path/to/foo.so -dll path/to/bar.so"
     *
     */
    Transformator(const char* modFiles);

    /*!
     * \brief dtor
     */
    ~Transformator();

    /*!
     * \brief loads a hoc/NEURON geometry and initializes the environment by calling prepare()
     *
     * \param[in] file path to the hoc file (with OS specific delimiters)
     */
    void load_geom(const std::string& file);

    /*!
     * \brief loads a stimulation protocol, nothing more.
     */
    void load_stim(const std::string& file);

    /*!
     * \brief get's number of sections in hoc file
     */
    size_t get_sections() const;

    /*!
     * \brief get number of all points
     */
    size_t get_total_points() const;

    /*!
     * \brief get tstart
     */
    number get_tstart();

    /*!
     * \brief get tstart
     */
    number get_tstop();

    /*!
     * \brief get tstart
     */
    number get_dt();

    /*!
     * \brief get tstart
     */
    number get_finitialize() const;

    /*!
     * \brief get t
     */
    number get_t();

    /*!
     * \brief get timestep
     */
    size_t get_limit() const;

#ifdef MPMNEURON_REVISION
    /*!
     * \brief get all sections of current file
     */
    std::vector<std::string> get_all_sections();
    std::string get_all_sections_as_string();
#endif

    /*!
     * \brief adjust hoc run protocol / setup
     */
    void setup_hoc(number tstart, number tstop, number dt, number finitialize=-75.0);

    /*!
     * \brief execute one single hoc statement
     *
     * Please note, if any value is returned from hoc, we will return it as double
     */
    number execute_hoc_stmt(const std::string& stmt) const;

    /*!
     * Set a hoc variable
     */
    bool set_hoc_variable(const std::string& var, number value);


    /*!
     * Get a hoc variable
     */
    number get_hoc_variable(const std::string& var);

    /*!
     * Set a hoc variable for a section
     */
    bool set_hoc_variable_sec(const std::string& var, number value, const std::string& sec);

    /*!
     * Get a hoc variable for a section
     */
    number get_hoc_variable_sec(const std::string& var, const std::string& sec);

    /*!
     * Next timestep
     */
    bool fadvance();

    /*!
     * get the transformator
     */
    SmartPtr<Transformator> get_transformator() {
      SmartPtr<Transformator> temp(this);
      return temp;
    }

    /*!
     * \brief execute multiple single hoc statements
     *
     * Please note, if any value is returned from hoc, we will return it as a vector of doubles
     */
    std::vector<number> execute_hoc_stmts(const std::vector<std::string>& stmts) const;


    /*!
     * \brief get on the fly the membrane potential for limit next timestep
     *
     * \param[in] limit the number of timesteps to be extracted
     * \param[in] the number of timesteps which should be advanced between in NEURON
     *
     * \return \c vms for the next limit timesteps
     */
    void extract_vms(size_t limit=1, size_t steps=1);


    /*!
     * \brief returns the vms for limit next timesteps
     */
    std::vector<std::vector<std::pair<std::vector<number>, number > > > get_vms(size_t limit=1) const;

    /*!
     * \brief purge the hoc interpreter environment
     *
     * \param[in] reinit if true reinits the environment for the transformator to initial state
     */
    void purge();

    /*!
     * print setup information
     */
    void print_setup(bool verbose=true);

  private:
    /*!
     * \brief initializes the hoc interpreter in the default case
     *
     * \return \c success status
     */
    bool init();

    /*!
     * \brief initializes the hoc interpreter in the enhanced case
     *
     * \return \c success status
     */
    bool init(int argc, char* argv[], char* env[]);

    /*!
     * \brief initializes teh hoc interpreter with the mod files
     */
    bool init(const char* modFiles);

    /*!
     * \brief prepare the hoc interpreter for command execution
     */
    void prepare();

};

} // namespace membrane_potential_mapping
} // namespace ug

#include "transformator_impl.h"

#endif // UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__TRANSFORMATOR_H
