# Custom implementation of the Oseen and Rotne-Prager-Yamakawa method into the framework of ESPResSo
This work represents an extension of the Extensible Simulation Package for Research on Soft Matter Systems (ESPResSo). The [user guide](https://espressomd.github.io/doc/index.html) of ESPResSo will walk you through the basic usage of ESPResSo.

Users of this custom implementation have to be aware of modification to the UI in regards to the Brownian dynamics routine of ESPResSo.

The main changes are:

New features avaiable: OSEEN and OSEEN_RPY
The Brownian dynamics integrator now accepts the particle diameter as input parameter.
The build/test_task folder contains multiple scripts used for analysis.

## License
This is a free extension developed for the software ESPResSo.  Feel free to redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
option) any later version <https://www.gnu.org/licenses/>.
