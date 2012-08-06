/**
 * \file transform.cpp
 * \brief implementation of preprocessing: uses high-level embedding of the Python language
 *
 * \see http://docs.python.org/extending/embedding.html
 *
 * \date created on Aug 3, 2012
 * \author Stephan Grein
 */


// includes
extern "C" {
#include "Python.h"
}
#include <boost/python.hpp>

#include <exception>
#include <sstream>
#include <string>

#include "transform.h"

#include <common/log.h>
#include <common/error.h>


// using directives
using namespace ug::membrane_potential_mapping;


void Transform::extract_timesteps_and_obj(const bool gen_objfile, const std::string neugen_executable, const std::string neutria_executable) {
	try {

		system("python foo.py");

		// initializes the Python interpreter
		Py_Initialize();
		/*PyRun_SimpleString("from time import ctime");
		PyRun_SimpleString("import sys");*/

		PyObject * mainModule = PyImport_AddModule("__main__");
		PyObject * hashlibModule = PyImport_ImportModule("neuron");
		PyModule_AddObject(mainModule, "neuron", hashlibModule);
		// fails because we compile with clang, whereas python modules are compiled with gcc!
		PyRun_SimpleString("print neuron");

		// import statements for python interface to NEURON (TODO: make sure neuron is in the user's path somehow...), check return value of PyRun_SimpleString!

		/*if (PyRun_SimpleString("from neuron import h") == -1)
			throw;

		if (PyRun_SimpleString("from neuron import hoc") == -1)
			throw;


		// TODO: use PyRun_SimpleFile for whole script (first create file into tmp dir with boost!)

		// command string buffer
		std::stringstream command;

		// 1st NEURON cmd
		command << "h.load_file"
				<< m_hocfile
				<< ")";
		PyRun_SimpleString(command.str().c_str());
		command.clear();

		// 2nd NEURON cmd
		PyRun_SimpleString("h.load_file(mview.hoc)");

		// 3rd NEURON cmd
		PyRun_SimpleString("h.define_shape()");

		// 4th NEURON cmd
		PyRun_SimpleString("modelView = h.ModelView(0)");

		// 5th NEURON cmd
		PyRun_SimpleString("modelxml = h.ModelViewXML(modelView)");

		// 6th NEURON cmd
		command << "modelxml.xportLevel1(" << m_xmlfile << ")";
		PyRun_SimpleString(command.str().c_str());
		command.clear();

		// 1st hoc cmd
		command << "hoc.execute('"
				<< "dt="
				<< m_dt
				<< ')\'\\n';
		PyRun_SimpleString(command.str().c_str());
		command.clear();

		// 2nd hoc cmd
		command << "hoc.execute('"
				<< "tstop"
				<< m_steps*m_dt
				<< ')\'\\n';
		PyRun_SimpleString(command.str().c_str());
		command.clear();

		// 3rd hoc cmd
		command << "hoc.execute('"
				<< "finitialize("
				<< m_vinit
				<< ')'
				<< ')\'\\n';
		PyRun_SimpleString(command.str().c_str());
		command.clear();

		// 4th hoc cmd
		command << "hoc.execute('"
				<< "while(t < tstop) { " << "\\n"
				<< "sprint(fname, timestep%f.csv, t/1000) \\n"
				<< "outfile = wopen(fname, \'w\') \\n"
				<< "forall for i=0,n3d()-1 fprint(\"%f %f %f %f\n\", x3d(i), y3d(i), z3d(i), v(i)) \\n"
				<< "wopen() \\n"
				<< "fadvance() \\n"
				<< "} \\n"
				<< '\')\\n';
		PyRun_SimpleString(command.str().c_str());
		command.clear();

		// shuts the Python interpreter down
		Py_Finalize();

		if (gen_objfile) {
			#ifdef __linux__
				command << "java -jar " << neugen_executable << " " << m_xmlfile << std::endl;
				system(command.str().c_str());
				command.clear();
				// TODO: use boost to create files and a temporary directory! write to file content with ofstream out("file")!
				command << "cat > neutria.lua <<NEUTRIA" << std::endl
						<< "neutria.set_hoc_params_low()" << std::endl
						<< "neutria.open_simple_hoc('test.shoc')" << std::endl
						<< "neutria.create_grid_hoc_spline()" << std::endl
						<< "neutria.save_grid_to_file('" << m_objfile << "')" << std::endl
						<< "NEUTRIA" << std::endl;
				system(command.str().c_str());
				command.clear();

				command << neutria_executable << " neutria.lua" << std::endl;
				system(command.str().c_str());
				command.clear();
			#elif __APPLE__
				std::cout << "APPLE" << std::endl;
			#elif __WIN32__
				std::cout << "WIN32" << std::endl;
			#endif
 			// TODO: implement .obj generation with NeuGen and NeuTria
		}*/
	} catch (const std::exception& exception) {
		UG_THROW("Fatal error in Transform::extract_timesteps_and_obj occured. Stopping execution.");
	}
}
