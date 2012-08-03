/**
 * \file transform.cpp
 * \brief implementation of preprocessing: transform (.hoc -> .obj) and extract timesteps

 * \date created on Aug 3, 2012
 * \author Stephan Grein
 */


// includes
extern "C" {
#include "Python.h"
}

#include <exception>
#include <sstream>
#include <string>

#include "transform.h"

#include <common/log.h>
#include <common/error.h>


// using directives
using namespace ug::membrane_potential_mapping;


void Transform::extract_timesteps_and_obj(const bool gen_objfile) {
	try {
		Py_Initialize();
		std::stringstream command;
		command << "h.load_file" << m_hocfile << ")";
		PyRun_SimpleString("from neuron import h");
		PyRun_SimpleString(command.str().c_str());
		command.clear();
		PyRun_SimpleString("h.load_file(mview.hoc)");
		PyRun_SimpleString("h.define_shape()");
		PyRun_SimpleString("modelView = h.ModelView(0)");
		PyRun_SimpleString("modelxml = h.ModelViewXML(modelView)");
		command << "modelxml.xportLevel1(" << m_xmlfile << ")";
		PyRun_SimpleString(command.str().c_str());
		command.clear();

		command << "h.dt =" << m_dt;
		PyRun_SimpleString(command.str().c_str());
		command.clear();
		command << "h.tstop =" << m_steps*m_dt;
		PyRun_SimpleString(command.str().c_str());
		command.clear();
		command << "h.finitialize(" << m_vinit << ")";
		PyRun_SimpleString(command.str().c_str());
		// TODO: check while loop!
		command.clear();
		command << "while (h.t <= h.tstop): \n\t" << "sprint(fname, timestep%f.csv, t/1000)\n\t";
		command << "open(fname, 'w') \n\t";
		command << "for i in h.n3d()-1:\n\t";
		command << "fprintf('%f %f %f %f\n', h.x3d(i), h.y3d(i), h.z3d(i), v(i))\n\t";
		command << "fname.close() \n\t" << "h.fadvance()\n\n";
		PyRun_SimpleString(command.str().c_str());
		command.clear();

		if (gen_objfile) {
			// TODO: implement .obj generation with NeuGen and NeuTria
		}
		Py_Finalize();
	} catch (const std::exception& exception) {
		UG_THROW("Fatal error in Transform::extract_timesteps_and_obj occured. Stopping execution.");
	}
}
