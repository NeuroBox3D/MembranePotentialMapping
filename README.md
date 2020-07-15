# MembranePotentialMapping
UG4 plugin providing functionality to map membrane potential values calculated 
with the NEURON simulator to 3D plasma membrane surfaces for simulations with ug4.

## Installation
Enable the plugin by `cmake -DMembranePotentialMapping=ON` during ug4 configuration
and provide NEURON installation directory by `-DMPMNEURON=/path/to/your/local/neuron/installation`.

The user can decide which version to use, for compability reasons a stable version (v7.3)
which works for all platforms is used, but the user can specify also to use the HEAD
revision of the NEURON development branch by `-DMPMNEURON_REVISION=ON` or the v7.3 
branch by `-DMPMNEURON_REVISION=OFF` (default).

### Enable the Boost testsuite
To enable the Boost testsuite issue the following command: `cmake -D MPMTestsuite=ON`.

### Issues
Static builds will prevent NEURON's `.mod` files from being loaded, thus use a dynamic build
by `-DSTATIC_BUILD=OFF`.

## How to use
There are two ways to do this:
1. Running Neuron in advance and writing the membrane potential values to file. They can then be read by the `Vm2uGMPM` class.
2. Running Neuron at the same time as the 3D simulation. The membrane potential values can then be mapped using the `NeuronMPM` class.
