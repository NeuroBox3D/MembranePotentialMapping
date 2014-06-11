/*!
 * transformator.h
 *
 * TODO: optimize code in general
 * TODO: namings to be optimized
 * TODO: cleanup code in general again
 *
 * TODO: use mathvector of ug and introduce template arg
 * TODO: eliminate code duplication as a next step
 *
 *  Created on: Nov 6, 2013
 *      Author: stephangrein
 */

#ifndef __H__UG__PLASMA_MEMBRANE__TRANSFORMATOR__
#define __H__UG__PLASMA_MEMBRANE__TRANSFORMATOR__

#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <iostream>

#include "common/types.h"
#define nil NULL

#include <lib_grid/common_attachments.h>
#include <lib_grid/grid/grid_util.h>

#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/function_spaces/grid_function.h"
/*
// define nil (how to circumvent that define?)
#define nil NULL
// necessary NEURON includes
#include "/home/stephan/code/hg/nrn/nrn-7.3/src/ivoc/oc2iv.h"
#include "/home/stephan/code/hg/nrn/nrn-7.3/src/ivoc/ocjump.cpp"
#include "/home/stephan/code/hg/nrn/nrn-7.3/src/ivoc/ivocmain.cpp"

// extern C functions (NEURON)
extern bool hoc_valid_stmt(const char* stmt, Object* ob);
extern int ivocmain(int, char**, char**);
//void ivoc_cleanup() { return; } // UNDEF this afterwards TODO (when fixing includes in ivoc.cpp it becomes available!!!)
extern double hoc_ac_;
*/

// namespace ug
namespace ug {
	// namespace pm
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
				 * \brief enhanced ctor
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
				 * ���\brief get tstart
				 */
				number get_tstart();

				/*!
				 * ���brief get tstart
				 */
				number get_tstop();

				/*!
				 * ���\brief get tstart
				 */
				number get_dt();

				/*!
				 * ���\brief get tstart
				 */
				number get_finitialize() const;

				/*!
				 * \brief get t
				 */
				number get_t();

				/*!
				 * ���\brief get timestep
				 */
				size_t get_limit() const;

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
				 * ���brief returns the vms for limit next timesteps
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
				 * \brief prepare the hoc interpreter for command execution
				 */
				void prepare();

		};
	}
}
#include "transformator_impl.h"

#endif /* __H__UG__PLASMA_MEMBRANE__TRANSFORMATOR__ */
