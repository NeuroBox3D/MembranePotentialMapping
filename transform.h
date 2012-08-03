/**
 * \file transform.h
 * \brief TODO: preprocessing: transform (.hoc -> .obj) and extract timesteps
 *
 * \date created on Aug 3, 2012
 * \author Stephan Grein
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__

// includes
#include <string>


// begin namespace ug (ug)
namespace ug {
	// begin namespace membrane_potential_mapping (mpm)
	namespace membrane_potential_mapping {
		class Transform {
			public:
				/**
				 * \brief default constructor for the Transform class
				 *
				 * \param[in] dt the delta_t (0.1) [ms]
				 * \param[in] number of steps (100) [#]
				 * \param[in] vinit the initial membrane potential (= resting potential) (-75.0) [mV]
				 *
				 */
				Transform(const std::string hocfile, const double dt=0.1, const long steps=100, const double vinit=-75.0) : m_hocfile(hocfile), m_objfile(""), m_xmlfile(""), m_timestepfile(""), m_dt(dt), m_steps(steps), m_vinit(vinit) { }
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
				inline void modify_hoc_setup(const double dt, const long steps, const double vinit)  {
					m_dt = dt;
					m_steps = steps;
					m_vinit = vinit;
				}

				/**
				 * \brief extracts the timesteps and the object file (.obj)
				 *
				 * \return \c return the path to the .obj file
				 */
				void extract_timesteps_and_obj(const bool gen_objfile=false);

				// getter
				inline std::string& get_hocfile() {
					return m_hocfile;
				}

				inline std::string& get_objfile() {
					return m_objfile;
				}

				inline std::string& get_xmlfile() {
					return m_xmlfile;
				}

				inline std::string& get_timestepfile() {
					return m_timestepfile;
				}

				inline double get_dt() {
					return m_dt;
				}


			private:
				std::string m_hocfile;
				std::string m_objfile;
				std::string m_xmlfile;
				std::string m_timestepfile;

				double m_dt;
				long m_steps;
				double m_vinit;
		};
	}
}

#endif /* __H__UG__MEMBRANE_POTENTIAL_MAPPING__TRANSFORM__ */
