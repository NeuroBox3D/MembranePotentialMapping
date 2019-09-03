# cable_neuron
UG4 plugin providing functionality to map membrane potential values calculated with Neuron
to 3D plasma membrane surfaces for simulations with UG 4.

There are two ways to do this:
Running Neuron in advance and writing the membrane potential values to file.
They can then be read by the Vm2uGMPM class.
Running Neuron at the same time as the 3D simulation. The membrane potential values can
then be mapped using the NeuronMPM class.