/*!
 * \file transform.h
 * \brief header for the preprocessing step
 * transform (.hoc -> .obj) and
 * extract the associated membrane potential for each timestep
 * of the NEURON simulation setup (specified in .hoc)
 * \addtogroup mpm_plugin
 *
 * \date created on Aug 3, 2012
 * \author Stephan Grein
 *
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__

#include <string>

#include "/home/stephan/code/hg/nrn/nrn-7.3/src/ivoc/oc2iv.h"
#include "/home/stephan/code/hg/nrn/nrn-7.3/src/ivoc/ocjump.cpp"
#include "/home/stephan/code/hg/nrn/nrn-7.3/src/ivoc/ivocmain.cpp"

// begin namespace ug
namespace ug {
	/* begin namespace mpm */
	namespace membrane_potential_mapping {
		class Transform {
			public:
				/*!
				 * \brief default constructor for the Transform class
				 * \addtogroup mpm_plugin
				 *
				 * \param[in] hocfile the hoc file which shall be transformed
				 * \param[in] timestep_directory location to which the timesteps should be written
				 * \param[in] dt the delta_t (0.1) [ms]
				 * \param[in] steps number of steps (100) [#]
				 * \param[in] vinit the initial membrane potential (= resting potential) (-75.0) [mV]
				 *
s.erase(s.find_last_of("."), string::npos);
				 */
				Transform(const std::string& hocfile, const std::string& timestep_directory, number dt=0.1, long steps=100, number vinit=-75.0) : m_hocfile(hocfile), m_objfile(std::string(hocfile).erase(hocfile.find_last_of("."), std::string::npos)), m_xmlfile(std::string(hocfile).erase(hocfile.find_last_of("."), std::string::npos)), m_timestepdirectory(timestep_directory), m_dt(dt), m_steps(steps), m_vinit(vinit) { }
				~Transform() { };

				/*!
				 * \brief modifies the hoc setup
				 *
				 * \param[in] dt the delta_t [ms]
				 * \param[in] steps number of steps [#]
				 * \param[in] vinit the initial membrane potential (= resting potential) (-75.0) [mV]
				 *
				 * \return \c bool true if preparation succeeded
				 */
				void modify_hoc_setup(number dt, long steps, number vinit);

				/*!
				 * \brief extracts the timesteps and the object file (.obj)
				 *
				 * \param[in] gen_objfile perform second optional step iff True
				 * \param[in] neugen_executable the location of the neugen executable
				 * \param[in] neutria_executable the location of the neutria executable
				 *
				 * first step: extract membrane potentials associated with every timestep in the simulation setup
				 * second (optional) step: generate .obj file out of .hoc file iff gen_objfile iff True
				 *
				 * \return \c return the path to the .obj file
				 */
				void extract_timesteps_and_obj(bool gen_objfile=false, const std::string& neugen_executable="NeuGen3D", const std::string& neutria_executable="neutria");

				/*!
				 * \brief simple getter for the hoc file
				 *
				 * \return \c the hoc file as string
				 */
				 inline std::string get_hocfile() const {
					return m_hocfile;
				}

				/*!
				 * \brief simple getter for the obj file
				 *
				 * \return \c the obj file as string
				 */
				inline std::string get_objfile() const {
					return m_objfile;
				}

				/*!
				 * \brief simple getter for the xml file
				 *
				 * \return \c the xml file as string
				 */
				inline std::string get_xmlfile() const {
					return m_xmlfile;
				}

				/*!
				 * \brief simple getter for the timestep directory
				 *
				 * \return \c the timestep directory as string
				 */
				inline std::string get_timestepdirectory() const {
					return m_timestepdirectory;
				}

				/*!
				 * \brief simple getter for the timestep width
				 *
				 * \return \c the timestep width as double
				 */
				inline double get_dt() const {
					return m_dt;
				}


			private:
				std::string m_hocfile;
				std::string m_objfile;
				std::string m_xmlfile;
				std::string m_timestepdirectory;

				double m_dt;
				long m_steps;
				double m_vinit;
		};
	/* end namespace mpm */
	}
/* end namespace ug */
}

#endif /* __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__ */
