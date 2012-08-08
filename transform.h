/**
 * \file transform.h
 * \brief header for the preprocessing step: transform (.hoc -> .obj) and extract timesteps
 *
 * \date created on Aug 3, 2012
 * \author Stephan Grein
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__


// includes
#include <string>
#include <boost/filesystem.hpp>


// begin namespace ug (ug)
namespace ug {
	// begin namespace membrane_potential_mapping (mpm)
	namespace membrane_potential_mapping {
		class Transform {
			public:
				/**
				 * \brief default constructor for the Transform class
				 *
				 * \param[in] the hoc file which shall be transformed
				 * \param[in] timestep_directory location to which the timesteps should be written
				 * \param[in] dt the delta_t (0.1) [ms]
				 * \param[in] number of steps (100) [#]
				 * \param[in] vinit the initial membrane potential (= resting potential) (-75.0) [mV]
				 *
				 */
				Transform(const std::string& hocfile, const std::string& timestep_directory, const double dt=0.1, const long steps=100, const double vinit=-75.0) : m_hocfile(hocfile), m_objfile(boost::filesystem::path(hocfile).replace_extension(".obj").string()), m_xmlfile(boost::filesystem::path(hocfile).replace_extension(".xml").string()), m_timestepdirectory(timestep_directory), m_dt(dt), m_steps(steps), m_vinit(vinit) { }
				~Transform() { };

				/**
				 * \brief modifies the hoc setup
				 *
				 * \param[in] dt the delta_t [ms]
				 * \param[in] number of steps [#]
				 * \param[in] vinit the initial membrane potential (= resting potential) (-75.0) [mV]
				 *
				 * \return \c bool true if preparation succeeded
				 */
				void modify_hoc_setup(const double dt, const long steps, const double vinit);

				/**
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
				void extract_timesteps_and_obj(const bool gen_objfile=false, const std::string& neugen_executable="NeuGen3D", const std::string& neutria_executable="neutria");

				// getter
				inline const std::string& get_hocfile() const {
					return m_hocfile;
				}

				inline const std::string& get_objfile() const {
					return m_objfile;
				}

				inline const std::string& get_xmlfile() const {
					return m_xmlfile;
				}

				inline const std::string& get_timestepdirectory() const {
					return m_timestepdirectory;
				}

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
	}
}

#endif /* __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__ */
