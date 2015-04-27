/*!
 * transformator.cpp
 *
 * TODO: use ForAllSections(sec) construct to replace the cumbersome for loop constructs in the other methods if necessary!*
 * NOTE: Section* sec = chk_access() and secname(sec) will give the currently accessed session, this could be handy in future, e. g. for the clamp constructs
 * TODO: find out why purge does not work as expected!!! exit() is called on quit() hoc command, this is not what we want in advance
 *
 *  Created on: Nov 6, 2013
 *      Author: stephangrein
 */

#include "transformator.h"

// define nil (how to circumvent that define?)
#define nil NULL
// necessary NEURON includes
#include "oc2iv.h"
#include "ocjump.cpp"
#include "ivocmain.cpp"
#include "section.h"
// #include "cabcode.c" -> depends on python probably
// extern void hoc_final_exit();
// extern void ivoc_final_exit();

// necessary stdlib includes
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <iostream>

// ug includes
#include <common/log.h>
#include <common/error.h>
#include <lib_disc/function_spaces/integrate.h>

// extern C functions (NEURON)
extern bool hoc_valid_stmt(const char* stmt, Object* ob);
#ifndef MPMNEURON_REVISION
extern int ivocmain(int, char**, char**);
#else
extern int ivocmain(int, const char**, const char**);
extern const char* secname(Section*);
extern Section* chk_access(void);
// extern Section* nrn_noerr_access(); -> depends on python probably
#endif
//void ivoc_cleanup() { return; } // UNDEF this afterwards TODO (when fixing includes in ivoc.cpp it becomes available!!!)
extern double hoc_ac_;

// using namespace
using namespace ug::membrane_potential_mapping;

// initialize static members
int Transformator::ARGC = 0;
char* Transformator::ARGV[] = {0};
char* Transformator::ENV[] = {0};

/////////////////////////////////////////////////////////
// feedback
/////////////////////////////////////////////////////////
void Transformator::feedback(const char* uCmp, const char* subset, const char* subset_vol, double density) {
			ug::GridSubsetHandler sh(*mGrid);
			// size_t si = sh.get_subset_index(subset_vol);
			mGrid->attach_to_vertices(ug::aPosition);
			ug::Grid::VertexAttachmentAccessor<ug::APosition3> aaPos(*mGrid, ug::aPosition);

			// get NEURON parameters
						size_t no_sections = m_sections;
						for (size_t k = 0; k <= no_sections; k++) {
							hoc_valid_stmt("j = 0", 0);
							hoc_valid_stmt("current_section = new SectionList()", 0);
							std::stringstream cmd;
							cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
							hoc_valid_stmt(cmd.str().c_str(), 0);
							cmd.str("");
							bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
							size_t no_points = static_cast<size_t>(hoc_ac_);
								if (ret) {
									for (size_t j=0; j <= no_points; j++) {
										// stmt string
										std::stringstream stmt;

										// X
										stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number x3d = hoc_ac_;

										// Y
										stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number y3d = hoc_ac_;

										// Z
										stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number z3d = hoc_ac_;
										// iterate over all triangles in subset 0
										for(TetrahedronIterator iter = sh.begin<Tetrahedron>(0); iter != sh.end<Tetrahedron>(0); ++iter) {
											Tetrahedron* h = *iter;
											ug::vector3 neuron_point(x3d, y3d, z3d);
											if (ContainsPoint(h, neuron_point, aaPos)) {
												m_saved_neuron_point = neuron_point;

												//double Ca_i = 0.;
												// TODO: integrate on subset "subset" with integrate.h from Andreas Vogel
										//		 IIntegrand<number, 3> integrand;
											//	number Ca_i =  ug::Integrate<1,2, Tetrahedron>(sh.begin<Tetrahedron>(0), sh.end<Tetrahedron>(0), aaPos, integrand, 0, "");

											//	stmt << "forsec current_section { Cai = " << Ca_i << "}"; // for this to work specify channel model with NMDL languag of NEURON TODO
												break;
											}
										}
									} // for points
								} // ret no_points
							// remove every covered section of current_section SectionList ...
							hoc_valid_stmt("{current_section.remove(current_section)}", 0);
						} // for no_sections
		}

void Transformator::feedbacks(const char* uCmp, const char* subset, const char* subset_vol, double density) {
			ug::GridSubsetHandler sh(*mGrid);
			// size_t si = sh.get_subset_index(subset_vol);
			mGrid->attach_to_vertices(ug::aPosition);
			ug::Grid::VertexAttachmentAccessor<ug::APosition3> aaPos(*mGrid, ug::aPosition);

			// get NEURON parameters
						size_t no_sections = m_sections;
						for (size_t k = 0; k <= no_sections; k++) {
							hoc_valid_stmt("j = 0", 0);
							hoc_valid_stmt("current_section = new SectionList()", 0);
							std::stringstream cmd;
							cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
							hoc_valid_stmt(cmd.str().c_str(), 0);
							cmd.str("");
							bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
							size_t no_points = static_cast<size_t>(hoc_ac_);
								if (ret) {
									for (size_t j=0; j <= no_points; j++) {
										// stmt string
										std::stringstream stmt;

										// X
										stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number x3d = hoc_ac_;

										// Y
										stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number y3d = hoc_ac_;

										// Z
										stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number z3d = hoc_ac_;
										// iterate over all triangles in subset 0
										for(TetrahedronIterator iter = sh.begin<Tetrahedron>(0); iter != sh.end<Tetrahedron>(0); ++iter) {
											Tetrahedron* h = *iter;
											ug::vector3 neuron_point(x3d, y3d, z3d);
											if (ContainsPoint(h, neuron_point, aaPos)) {
	//											m_saved_neuron_point = neuron_point;
												mNeuronPoints.push_back(neuron_point);
												//double Ca_i = 0.;
												// TODO: integrate on subset "subset" with integrate.h from Andreas Vogel
										//		 IIntegrand<number, 3> integrand;
											//	number Ca_i =  ug::Integrate<1,2, Tetrahedron>(sh.begin<Tetrahedron>(0), sh.end<Tetrahedron>(0), aaPos, integrand, 0, "");

											//	stmt << "forsec current_section { Cai = " << Ca_i << "}"; // for this to work specify channel model with NMDL languag of NEURON TODO
												break;
											}
										}
									} // for points
								} // ret no_points
							// remove every covered section of current_section SectionList ...
							hoc_valid_stmt("{current_section.remove(current_section)}", 0);
						} // for no_sections
		}



/// set feedback
void Transformator::set_feedback() {
			// get NEURON parameters
						size_t no_sections = m_sections;
						for (size_t k = 0; k <= no_sections; k++) {
							hoc_valid_stmt("j = 0", 0);
							hoc_valid_stmt("current_section = new SectionList()", 0);
							std::stringstream cmd;
							cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
							hoc_valid_stmt(cmd.str().c_str(), 0);
							cmd.str("");
							bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
							size_t no_points = static_cast<size_t>(hoc_ac_);
								if (ret) {
									for (size_t j=0; j <= no_points; j++) {
										// stmt string
										std::stringstream stmt;

										// X
										stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number x3d = hoc_ac_;

										// Y
										stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number y3d = hoc_ac_;

										// Z
										stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number z3d = hoc_ac_;

										ug::vector3 neuron_point(x3d, y3d, z3d);
										double Ca_i = 0;
										if (m_saved_neuron_point.x() == neuron_point.x() && m_saved_neuron_point.y() == neuron_point.y() && m_saved_neuron_point.z() == neuron_point.z()) {
											Ca_i = m_saved_neuron_point_integral;
										}
										stmt << "forsec current_section { cai = " << Ca_i << "}"; // for this to work specify channel model with NMDL language of NEURON
										break;
									} // for points
							} // ret no_points
						// remove every covered section of current_section SectionList ...
					hoc_valid_stmt("{current_section.remove(current_section)}", 0);
			} // for no_sections
}

/// set feedbacks
void Transformator::set_feedbacks() {
			// get NEURON parameters
						size_t no_sections = m_sections;
						for (size_t k = 0; k <= no_sections; k++) {
							hoc_valid_stmt("j = 0", 0);
							hoc_valid_stmt("current_section = new SectionList()", 0);
							std::stringstream cmd;
							cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
							hoc_valid_stmt(cmd.str().c_str(), 0);
							cmd.str("");
							bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
							size_t no_points = static_cast<size_t>(hoc_ac_);
								if (ret) {
									for (size_t j=0; j <= no_points; j++) {
										// stmt string
										std::stringstream stmt;

										// X
										stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number x3d = hoc_ac_;

										// Y
										stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number y3d = hoc_ac_;

										// Z
										stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number z3d = hoc_ac_;

										ug::vector3 neuron_point(x3d, y3d, z3d);
										double Ca_i = 0;
										size_t index = 0;
										for (std::vector<ug::vector3>::const_iterator cit = mNeuronPoints.begin(); cit != mNeuronPoints.end(); cit++) {

											if (cit->x() == neuron_point.x() && cit->y() == neuron_point.y() && cit->z() == neuron_point.z()) {
												Ca_i = mIntegrals[index];
												stmt << "forsec current_section { cai = " << Ca_i << "}"; // for this to work specify channel model with NMDL language of NEURON
												break;
											}
											index++;
										}
									} // for points
							} // ret no_points
						// remove every covered section of current_section SectionList ...
					hoc_valid_stmt("{current_section.remove(current_section)}", 0);
			} // for no_sections
}



/////////////////////////////////////////////////////////
// init feedback
/////////////////////////////////////////////////////////
void Transformator::init_feedback(const char* uCmp, const char* subset, const char* subset_vol, double density, const char* channel) {
			ug::GridSubsetHandler sh(*mGrid);
			// size_t si = sh.get_subset_index(subset_vol);
			mGrid->attach_to_vertices(ug::aPosition);
			ug::Grid::VertexAttachmentAccessor<ug::APosition3> aaPos(*mGrid, ug::aPosition);

			// get NEURON parameters
						size_t no_sections = m_sections;
						for (size_t k = 0; k <= no_sections; k++) {
							hoc_valid_stmt("j = 0", 0);
							hoc_valid_stmt("current_section = new SectionList()", 0);
							std::stringstream cmd;
							cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
							hoc_valid_stmt(cmd.str().c_str(), 0);
							cmd.str("");
							bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
							size_t no_points = static_cast<size_t>(hoc_ac_);
								if (ret) {
									for (size_t j=0; j <= no_points; j++) {
										// stmt string
										std::stringstream stmt;

										// X
										stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number x3d = hoc_ac_;

										// Y
										stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number y3d = hoc_ac_;

										// Z
										stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number z3d = hoc_ac_;
										// iterate over all triangles in subset 0
										for(TetrahedronIterator iter = sh.begin<Tetrahedron>(0); iter != sh.end<Tetrahedron>(0); ++iter) {
											Tetrahedron* h = *iter;
											ug::vector3 neuron_point(x3d, y3d, z3d);
											if (ContainsPoint(h, neuron_point, aaPos)) {
												stmt.clear();
												stmt << "forsec current_section { insert " << channel << "}";
												hoc_valid_stmt(stmt.str().c_str(), 0);
												break;
											}
										}
									} // for points
								} // ret no_points
							// remove every covered section of current_section SectionList ...
							hoc_valid_stmt("{current_section.remove(current_section)}", 0);
						} // for no_sections
		}

// init feedbacks
void Transformator::init_feedbacks(const char* uCmp, const char* subset, const char* subset_vol, double density, const char* channel) {
			ug::GridSubsetHandler sh(*mGrid);
			// size_t si = sh.get_subset_index(subset_vol);
			mGrid->attach_to_vertices(ug::aPosition);
			ug::Grid::VertexAttachmentAccessor<ug::APosition3> aaPos(*mGrid, ug::aPosition);

			// get NEURON parameters
						size_t no_sections = m_sections;
						for (size_t k = 0; k <= no_sections; k++) {
							hoc_valid_stmt("j = 0", 0);
							hoc_valid_stmt("current_section = new SectionList()", 0);
							std::stringstream cmd;
							cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
							hoc_valid_stmt(cmd.str().c_str(), 0);
							cmd.str("");
							bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
							size_t no_points = static_cast<size_t>(hoc_ac_);
								if (ret) {
									for (size_t j=0; j <= no_points; j++) {
										// stmt string
										std::stringstream stmt;

										// X
										stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number x3d = hoc_ac_;

										// Y
										stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number y3d = hoc_ac_;

										// Z
										stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
										hoc_valid_stmt(stmt.str().c_str(), 0);
										stmt.str("");
										number z3d = hoc_ac_;
										// iterate over all triangles in subset 0
										for(TetrahedronIterator iter = sh.begin<Tetrahedron>(0); iter != sh.end<Tetrahedron>(0); ++iter) {
											// Tetrahedron* h = *iter;
											ug::vector3 neuron_point(x3d, y3d, z3d);
													for (std::vector<ug::vector3>::const_iterator cit = mNeuronPoints.begin(); cit != mNeuronPoints.end(); cit++) {

											if (cit->x() == neuron_point.x() && cit->y() == neuron_point.y() && cit->z() == neuron_point.z()) {
												stmt << "forsec current_section { insert  " << channel << "}"; // for this to work specify channel model with NMDL language of NEURON
												hoc_valid_stmt(stmt.str().c_str(), 0);
												stmt.clear();
											}
													}
										}
									} // for points
								} // ret no_points
							// remove every covered section of current_section SectionList ...
							hoc_valid_stmt("{current_section.remove(current_section)}", 0);
						} // for no_sections
		}



/////////////////////////////////////////////////////////
// adjust_resisitivity_obstacle
/////////////////////////////////////////////////////////
void Transformator::adjust_resistivity_obstacle(const char* uCmp, const char* subset, double Rai) {
				ug::GridSubsetHandler sh(*mGrid);
				//size_t si = sh.get_subset_index(subset);
				mGrid->attach_to_vertices(ug::aPosition);
				ug::Grid::VertexAttachmentAccessor<ug::APosition3> aaPos(*mGrid, ug::aPosition);

			// get NEURON parameters
				size_t no_sections = m_sections;
				for (size_t k = 0; k <= no_sections; k++) {
					hoc_valid_stmt("j = 0", 0);
					hoc_valid_stmt("current_section = new SectionList()", 0);
					std::stringstream cmd;
					cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
					hoc_valid_stmt(cmd.str().c_str(), 0);
					cmd.str("");
					bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
					size_t no_points = static_cast<size_t>(hoc_ac_);
						if (ret) {
							for (size_t j=0; j <= no_points; j++) {
								// stmt string
								std::stringstream stmt;

								// X
								stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
								hoc_valid_stmt(stmt.str().c_str(), 0);
								stmt.str("");
								number x3d = hoc_ac_;

								// Y
								stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
								hoc_valid_stmt(stmt.str().c_str(), 0);
								stmt.str("");
								number y3d = hoc_ac_;

								// Z
								stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
								hoc_valid_stmt(stmt.str().c_str(), 0);
								stmt.str("");
								number z3d = hoc_ac_;
								// iterate over all triangles in subset 0
								for(TetrahedronIterator iter = sh.begin<Tetrahedron>(0); iter != sh.end<Tetrahedron>(0); ++iter) {
									Tetrahedron* h = *iter;
									ug::vector3 neuron_point(x3d, y3d, z3d);
									if (ContainsPoint(h, neuron_point, aaPos)) {
										double Rai_new = Rai;
										stmt << "forsec current_section { Ra = " << Rai_new << "}";
										break;
									}
								}
							} // for points
						} // ret no_points
					// remove every covered section of current_section SectionList ...
					hoc_valid_stmt("{current_section.remove(current_section)}", 0);
				} // for no_sections
}


Transformator::Transformator() : m_sections(0), m_totalPoints(0), m_finitialize(0.0), m_tstart(0), m_tstop(0), m_dt(0), m_t(0), m_limit(0){
	init(ARGC, ARGV, ENV);
}

Transformator::Transformator(int argc, char* argv[], char* env[]) : m_sections(0), m_totalPoints(0), m_finitialize(0.0), m_tstart(0), m_tstop(0), m_dt(0), m_t(0), m_limit(0) {
	init();
}

Transformator::Transformator(const char* modFiles) :  m_sections(0), m_totalPoints(0), m_finitialize(0.0), m_tstart(0), m_tstop(0), m_dt(0), m_t(0), m_limit(0) {
	init(modFiles);
}

bool Transformator::init(const char* modFiles) {
	// init ivocmain (hoc) interpreter
   #ifdef MPMNEURON_REVISION
	 const char* args[] = { modFiles };
	 const char* env[] = {0};
	 int init = ivocmain(1, args, env);
   #else
	 std::stringstream ss;
	 ss << modFiles;
	 char* args[] = { const_cast<char*>(ss.str().c_str()) };
	 char* env[] = {0};
	 int init = ivocmain(1, ARGV, ENV);
   #endif

	  // check for success
	 if (init != 0) {
		 UG_THROW("ivocmain (hoc interpreter) could not be initialized by Transformator::init." << std::endl);
	 }

 return true;
}

Transformator::~Transformator() {
	purge();
}

bool Transformator::init() {
	 // init ivocmain (hoc) interpreter
   const char* a[] = {0};
   const char* e[] = {0};
   
   #ifdef MPMNEURON_REVISION
	 int init = ivocmain(0, a, e);
   #else
   int init = ivocmain(ARGC, ARGV, ENV);
   #endif

	 // check for success
	 if (init != 0) {
		 UG_THROW("ivocmain (hoc interpreter) could not be initialized by Transformator::init." << std::endl);
	 }

	 return true;
}

bool Transformator::init(int argc, char* argv[], char* env[]) {
   const char* a[] = {0};
   const char* e[] = {0};
   
   #ifdef MPMNEURON_REVISION
	 int init = ivocmain(0, a, e);
   #else
   int init = ivocmain(ARGC, ARGV, ENV);
   #endif
	 // check for success
	 if (init != 0) {
		 UG_THROW("ivocmain (hoc interpreter) could not be initialized by Transformator::init." << std::endl);
	 }
	 // return status
	 return true;
 }

void Transformator::prepare() {
	 // add initializing statements
	 m_stmts.clear();
	 m_stmts.push_back("objref current_section");
	 m_stmts.push_back("current_section_index = 0"); // TODO: was this correct here previously we excluded this!?
	 m_stmts.push_back("current_section = new SectionList()");
	 m_stmts.push_back("forall { current_section_index = current_section_index + 1 }");
	 m_stmts.push_back("hoc_ac_ = current_section_index");

	 // error indicator
	 bool ret = true;

	 // execute all statements in accumulator
	 for (std::vector<std::string>::const_iterator it = m_stmts.begin(); it != m_stmts.end(); it++) {
	       ret = (ret && hoc_valid_stmt((*it).c_str(), 0));
	 }

	 // clear accumulator
	 m_stmts.clear();

	 // assign number of retrieved sections to member
	 m_sections = static_cast<size_t>(hoc_ac_);

	 // total points
	 size_t total_points = 0;
	   for (size_t k = 0; k < m_sections+1; k++) { // was +1
		   hoc_valid_stmt("j = 0", 0);
		   hoc_valid_stmt("current_section = new SectionList()", 0);
		   std::stringstream cmd;
		   cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
		   hoc_valid_stmt(cmd.str().c_str(), 0);
		   cmd.str("");
		   ret = (ret && hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0));
		   total_points += static_cast<size_t>(hoc_ac_);
	   }
	 m_totalPoints = total_points;

	 // check for succcess of all prepare statements
	 if (!ret) {
		 UG_LOG("At least one of the statements in Transformator::prepare did not succeed." << std::endl);
	 }
}

void Transformator::load_geom(const std::string& file) {
	if (m_geometry.compare(file) == 0) {
		UG_WARNING("NEURON/hoc interpreter does not allow reloading the same geometry file. You may encounter undefined behaviour!")
	}
	// create hoc command
	std::stringstream ss;
	ss << "{load_file(1, \"" << file << "\")}";

	// check for success
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);
	if (!ret) {
		UG_LOG("Could not load hoc file: " << file << " with Transformator::load_hoc." << std::endl);
	}

	// cleanup stringstream
	ss.str("");
	ss.flush();

	// set defining geometry file
	m_geometry = file;
	// set number of sections
	prepare();
}

void Transformator::load_stim(const std::string& file) {
	if (m_stimulation.compare(file) == 0) {
		UG_WARNING("NEURON/hoc interpreter does not allow reloading the same stimulation file. You may encounter undefined behaviour!")
	}
	// create hoc command
		std::stringstream ss;
		ss << "{load_file(1, \"" << file << "\")}";

		// check for success
		bool ret = hoc_valid_stmt(ss.str().c_str(), 0);
		if (!ret) {
			UG_LOG("Could not load hoc file: " << file << " with Transformator::load_hoc." << std::endl);
		}

		// cleanup stringstream
		ss.str("");
		ss.flush();

		// set defining stimulation file
		m_stimulation = file;
}

void Transformator::setup_hoc(number tstart, number tstop, number dt, number finitialize) {
	// command string
	std::stringstream ss;
	// return value
	bool ret = true;

	// tstart
	ss << "{tstart = " << tstart << "}";
	ret = (ret && hoc_valid_stmt(ss.str().c_str(), 0));
	ss.str("");

	// tstop
	ss << "{tstop = " << tstop << "}";
	ret = (ret && hoc_valid_stmt(ss.str().c_str(), 0));
	ss.str("");

	// dt
	ss << "{dt = " << dt << "}";
	ret = (ret && hoc_valid_stmt(ss.str().c_str(), 0));
	ss.str("");

	// finitialize
	ss << "{finitialize(" << finitialize << ")}";
	ret = (ret && hoc_valid_stmt(ss.str().c_str(), 0));
	ss.str("");
	m_finitialize = finitialize;

	// setup initial values (i. e. -75.0 mV)
	m_vms.clear();
	size_t no_sections = m_sections;
	std::vector<std::pair<std::vector<number>, number> > local_vms;

	for (size_t k = 0; k < no_sections+1; k++) {
		hoc_valid_stmt("j = 0", 0);
		hoc_valid_stmt("current_section = new SectionList()", 0);
		std::stringstream cmd;
		cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
		hoc_valid_stmt(cmd.str().c_str(), 0);
		cmd.str("");
		bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
		size_t no_points = hoc_ac_;
		if (ret) {
			for (size_t j=0; j < no_points; j++) {
				// stmt string
				std::stringstream stmt;

				// VM
				stmt << "forsec current_section { hoc_ac_ = v(" << j << ")}";
				hoc_valid_stmt(stmt.str().c_str(), 0);
				number vm = hoc_ac_;
				stmt.str("");

				// X
				stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
				hoc_valid_stmt(stmt.str().c_str(), 0);
				stmt.str("");
				number x3d = hoc_ac_;

				// Y
				stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
				hoc_valid_stmt(stmt.str().c_str(), 0);
				stmt.str("");
				number y3d = hoc_ac_;

				// Z
				stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
				hoc_valid_stmt(stmt.str().c_str(), 0);
				stmt.str("");
				number z3d = hoc_ac_;

				// create vector with membrane potential
				std::vector<number> coords;
				coords.push_back(x3d);
				coords.push_back(y3d);
				coords.push_back(z3d);
				std::pair<std::vector<number>, number> pair = std::pair<std::vector<number>, number>(coords, vm);
				local_vms.push_back(pair);
			} // for points
		} // ret no_points
		// remove every covered section of current_section SectionList ...
		hoc_valid_stmt("{current_section.remove(current_section)}", 0);
		} // for no_sections
	m_vms.push_back(local_vms);
	local_vms.clear();
}

bool Transformator::set_hoc_variable(const std::string& var, number value) {
	std::stringstream ss;
	ss.precision(std::numeric_limits<number>::digits10);
	// TODO bug here: http://www.cplusplus.com/reference/sstream/stringstream/ uses std::
	//http://www.cplusplus.com/reference/iomanip/setprecision/
	ss << var << " = " << value;
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);
	if (ret != 1) {
		UG_LOG("Could not execute hoc stmt: " << "var = " << value << " in Transformator::set_hoc_variable" << std::endl);
	}
	return ret;
}

bool Transformator::set_hoc_variable_sec(const std::string& var, number value, const std::string& section) {
	std::stringstream ss;
	ss << "access" << " " << section;
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);

	if (ret != 1) {
		UG_LOG("Could not execute hoc stmt: " << "access " << section << " in Transformator::set_hoc_variable_sec (Probably section does not exist)" <<  std::endl);
		return ret;
	}

	return set_hoc_variable(var, value);
}

number Transformator::get_hoc_variable(const std::string& var) {
	std::stringstream ss;
	ss << "hoc_ac_" << " = " << var;
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);
	if (ret != 1) {
		UG_LOG("Could not execute hoc stmt: " << "hoc_ac_ = " << var << " in Transformator::get_hoc_variable" << std::endl);
	}
	return static_cast<number>(hoc_ac_);
}

number Transformator::get_hoc_variable_sec(const std::string& var, const std::string& section) {
	std::stringstream ss;
	ss << "access" << " " << section;
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);

	if (ret != 1) {
		UG_LOG("Could not execute hoc stmt: " << "access " << section << " in Transformator::set_hoc_variable_sec (Probably section does not exist)" <<  std::endl);
		return ret;
	}

	return get_hoc_variable(var);
}

number Transformator::execute_hoc_stmt(const std::string& stmt) const {
	bool ret = hoc_valid_stmt(stmt.c_str(), 0);
	number val = 0.;

	if (ret != 1) {
		UG_LOG("Could not execute hoc stmt: " << stmt << " in Transformator::execute_hoc_stmt." << std::endl);
	}

	val = static_cast<number>(hoc_ac_);

	return val;
}

std::vector<number> Transformator::execute_hoc_stmts(const std::vector<std::string>& stmts) const {
	bool ret = true;
	std::vector<number> val;

	for (std::vector<std::string>::const_iterator it = stmts.begin(); it != stmts.end(); it++) {
		ret = (ret && hoc_valid_stmt(it->c_str(), 0));
		val.push_back(static_cast<number>(hoc_ac_));
	}

	if (ret != 1) {
		UG_LOG("Could not execute at least one of the hoc stmts in Transformator::execute_hoc_stmts." << std::endl);
	}

	return val;
}

bool Transformator::fadvance() {
	bool ret = hoc_valid_stmt("fadvance()", 0);

	if (!ret) {
		UG_LOG("Could not execute the fadvance() command in Transformator::fadvance()." << std::endl);
	}

	return ret;
}

size_t Transformator::get_limit() const {
	return m_limit;
}

size_t Transformator::get_sections() const {
	return m_sections;
}

size_t Transformator::get_total_points() const {
	return m_totalPoints;
}

number Transformator::get_tstart() {
	// command string
	std::stringstream ss;
	ss << "hoc_ac_ = tstart";
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);
	if (ret != 1) {
		UG_LOG("Could not acquire value of tstart by Transformator::get_tstart." << std::endl);
	}

	number tstart = static_cast<number>(hoc_ac_);
	m_tstart = tstart;
	return tstart;
}

number Transformator::get_tstop() {
	// command string
	std::stringstream ss;
	ss << "hoc_ac_ = tstop";
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);
	if (ret != 1) {
		UG_LOG("Could not acquire value of tstop by Transformator::get_tstop." << std::endl);
	}

	number tstop = static_cast<number>(hoc_ac_);
	m_tstop = tstop;
	return tstop;
}

number Transformator::get_dt()  {
	// command string
	std::stringstream ss;
	ss << "hoc_ac_ = dt";
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);
	if (ret != 1) {
		UG_LOG("Could not acquire value of dt by Transformator::get_dt." << std::endl);
	}

	number dt = static_cast<number>(hoc_ac_);
	m_dt = dt;
	return dt;
}

number Transformator::get_finitialize() const {
	return m_finitialize;
}

number Transformator::get_t() {
	// command string
	std::stringstream ss;
	ss << "hoc_ac_ = t";
	bool ret = hoc_valid_stmt(ss.str().c_str(), 0);
	if (ret != 1) {
		UG_LOG("Could not acquire value of t by Transformator::get_t." << std::endl);
	}

	number t = static_cast<number>(hoc_ac_);
	m_t = t;
	return t;
}


void Transformator::extract_vms(size_t limit, size_t steps) {
	m_limit += limit;
	m_vms.clear();
	size_t no_sections = m_sections;
	size_t no_steps = steps;
	for (size_t i=1; i <= limit; i++) {
		std::vector<std::pair<std::vector<number>, number> > local_vms;

		// fadvance multiple times
		bool ret_outer = true;
		for (size_t l=1; l <= no_steps; l++) {
			ret_outer = (ret_outer && hoc_valid_stmt("{fadvance()}", 0));
		}

		if (ret_outer) {
			for (size_t k = 0; k <= no_sections; k++) {
				hoc_valid_stmt("j = 0", 0);
				hoc_valid_stmt("current_section = new SectionList()", 0);
				std::stringstream cmd;
				cmd << "forall{ if (j == " << k << ") { {current_section.append() } } { j = j + 1} }";
				hoc_valid_stmt(cmd.str().c_str(), 0);
				cmd.str("");
				bool ret = hoc_valid_stmt("forsec current_section { hoc_ac_ = n3d() - 1}", 0);
				size_t no_points = static_cast<size_t>(hoc_ac_);
				if (ret) {
					for (size_t j=0; j <= no_points; j++) {
						// stmt string
						std::stringstream stmt;

						// VM
						stmt << "forsec current_section { hoc_ac_ = v(" << j << ")}";
						hoc_valid_stmt(stmt.str().c_str(), 0);
						number vm = hoc_ac_;
 						stmt.str("");

						// X
						stmt << "forsec current_section { hoc_ac_ = x3d(" << j << ")}";
						hoc_valid_stmt(stmt.str().c_str(), 0);
						stmt.str("");
						number x3d = hoc_ac_;

						// Y
						stmt << "forsec current_section { hoc_ac_ = y3d(" << j << ")}";
						hoc_valid_stmt(stmt.str().c_str(), 0);
						stmt.str("");
						number y3d = hoc_ac_;

						// Z
						stmt << "forsec current_section { hoc_ac_ = z3d(" << j << ")}";
						hoc_valid_stmt(stmt.str().c_str(), 0);
						stmt.str("");
						number z3d = hoc_ac_;

						// create vector with membrane potential
						std::vector<number> coords;
						coords.push_back(x3d);
						coords.push_back(y3d);
						coords.push_back(z3d);
						std::pair<std::vector<number>, number> pair = std::pair<std::vector<number>, number>(coords, vm);
						local_vms.push_back(pair);
					} // for points
				} // ret no_points
			// remove every covered section of current_section SectionList ...
			hoc_valid_stmt("{current_section.remove(current_section)}", 0);
			} // for no_sections
		} // if ret_outer
	m_vms.push_back(local_vms);
	local_vms.clear();
	} // for limit
}

std::vector<std::vector<std::pair<std::vector<number>, number> > > Transformator::get_vms(size_t limit) const {
	std::vector<std::vector<std::pair<std::vector<number>, number> > >::const_iterator begin = m_vms.begin();
	std::vector<std::vector<std::pair<std::vector<number>, number> > >::const_iterator end = m_vms.begin() + limit;
	if (limit > m_vms.size()) {
		return m_vms;
	} else {
		return std::vector<std::vector<std::pair<std::vector<number>, number> > >(begin, end);
	}
}

void Transformator::purge() {
	m_vms.clear(); // TODO: was this necessary here? previously this was not in here ...


	hoc_valid_stmt("forall pt3dclear()", 0);
	hoc_valid_stmt("forall delete_section()", 0);
	m_sections = 0;
	m_totalPoints = 0;
  //m_limit = 0;
	m_finitialize = -75.0;

	bool ret = hoc_valid_stmt("current_section_index = 0", 0);
		if (ret != 1) {
			UG_LOG("Could not purge the hoc interpreter by Transformator::purge. " << std::endl);
		}
	//	hoc_final_exit();
 //ivoc_final_exit();
}

void Transformator::print_setup(bool verbose) {
	this->get_t();
	this->get_dt();
	this->get_tstart();
	this->get_tstop();
	if (verbose) {
		std::cout << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "NEURON (hoc) setup:" << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "Geometry" << std::endl;
		std::cout << "No. sections:                     " << m_sections << " [#]" << std::endl;
		std::cout << "No. points:                       " << m_totalPoints << " [#]" << std::endl;
		std::cout <<  std::endl;
		std::cout << "Simulation protocol" << std::endl;
		std::cout << "NEURON t:                         " << m_t << " [s]" << std::endl;
		std::cout << "NEURON dt:                        " << m_dt  << " [s]" << std::endl;
		std::cout << "NEURON tstart                     " << m_tstart  << " [s]" << std::endl;
		std::cout << "NEURON tstop:                     " << m_tstop << " [s]" << std::endl;
		std::cout << "NEURON finitialize:               "  << m_finitialize << " [mV]" << std::endl;
		std::cout <<  std::endl;
		std::cout << "Defining files" << std::endl;
		std::cout << "Stimulation protocol:             " << m_stimulation << std::endl;
		std::cout << "Geometry:                         " << m_geometry << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout <<  std::endl;
	} else {
		std::cout << "NEURON (hoc) setup:" << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "Defining files" << std::endl;
		std::cout << "Stimulation protocol:             " << m_stimulation << std::endl;
		std::cout << "Geometry:                         " << m_geometry << std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
	}
}


		/////////////////////////////////////////////////////////
		/// get_all_sections
	    /////////////////////////////////////////////////////////
		#ifdef MPMNEURON_REVISION
		std::vector<std::string> Transformator::get_all_sections() {
			// number of sections known to the current instance of the hoc interpreter
			size_t no_sections = m_sections;

			// store the section names we retrieve below
			std::vector<std::string> section_names;

			if (no_sections == 0) {
				UG_LOG("Transformator::get_all_sections: #sections is zero -> therefore we have no sections and no section names!");
				return std::vector<std::string>();
			}

			// qsec, sec are required somehow for the looping
			hoc_Item* qsec;
			Section* sec;
			ForAllSections(sec) // omit opening brace due to MACRO definition! {
				section_names.push_back(secname(sec));
			}

			// return all names of present sections
			return section_names;
		}

	std::string Transformator::get_all_sections_as_string() {
					// number of sections known to the current instance of the hoc interpreter
					size_t no_sections = m_sections;

					// store the section names we retrieve below
					std::vector<std::string> section_names;

					if (no_sections == 0) {
						UG_LOG("Transformator::get_all_sections: #sections is zero -> therefore we have no sections and no section names!");
						return "";
					}

					// qsec, sec are required somehow for the looping
					hoc_Item* qsec;
					Section* sec;
					ForAllSections(sec) // omit opening brace due to MACRO definition! {
						section_names.push_back(secname(sec));
					}

					std::stringstream ss;
					for (std::vector<std::string>::const_iterator it = section_names.begin(); it != section_names.end(); it++) {
						ss << *it << ";";
					}

					// return all names of present sections
					return ss.str();
				}



		#endif

