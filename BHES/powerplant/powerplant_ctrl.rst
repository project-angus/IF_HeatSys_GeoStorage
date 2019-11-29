Power plant control file
++++++++++++++++++++++++

The power plant control file contains the following information:

* The names of the power plant models (identifier for the respective model specifications)
	* hp: heat pump
	* ice: internal combustion engine
	* ccbp: combined cycle with back pressure steam turbine
	* eb: electrode boiler
	* plb: peak load boiler
	
The model specifications are:

* path: path to the power plant model
* T_ff: connection (in power plant model) with feed flow temperature information
* T_ff: connection with back flow temperature information
* m: connection with mass flow information
* heat_bus: bus with transferred heat information
* power_bus: bus with power generation/consumption information
