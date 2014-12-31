![logo](parscale_logo.png)
======
A Compilation of Particle Scale Models.

ParScale is part of the [NanoSim Project](http://sintef.no/NanoSim)


Graz University of Technology and DCS Computing GmbH release ParScale
------------------
31st of Decemmber 2014

The Graz University of Technology (TU Graz) together with DCS Computing GmbH (DCS) are please to announced the release of the `1.0.1-beta` version of the tool `ParScale`. This version is distributed is licensed under the [Lesser General Public License](http://www.gnu.org/licenses/lgpl.html) by TU Graz and DCS.

###Features
Version 1.0.1-beta is the first public release of the ParScale library and is meant to introduce its features to a wide audience of users and possible developers of this library. Specifically, the current version of ParScale is able to

- simulate heat and mass transport in one or multiple porous particles.
- perform calculations on an equi-distant spherical one-dimensional coordinate system, or using a simplified model for reactive transport (i.e., the Shrinking Core model)
- predict the outcome of a single non-isothermal irreversible (homogeneous or heterogeneous) chemical reaction (full coupling of temperature and species concentration fields) involving an arbitrary number of solid and gas-phase species.
- consider various "grain-scale models" to picture the formation of a product layer on solids present in the porous of a particle.
- couple with a JSON file (i.e., read data from a JSON file to simulate transient boundary conditions or model parameters), as well as a LIGGGHTS simulation.

Test and verification cases, as well as documentation is supplied with ParScale. More test cases will be added step-by-step during the next years.

Note that ParScale is still under active development, and the success of ParScale relies to a large extend on contributions from users. Specifically, we plan to release

- modules to picture multiple reversible reactions that are specified by the user via CHEMKIN-II files,
- an improved version that uses orthogonal collocation methods for more efficient spatial discretization,
- refined sub-models for multi-component diffusion and temperature-dependent transport properties.

The developers of ParScale are thankful for any comments, ideas, or contributions to currently available or future features of ParScale that are posted in the ParScale user forum](http://www.cfdem.com/forum).

###Important Notes
'beta' means, that this version of ParScale has been tested internally to some extend, and now seeks for testing and application by external groups (i.e., non TUG or DCS). Thus, while TUG and DCS have carefully tested ParScale, it is likely that users might find significant bugs, run into unexpected errors, or might miss information on how to use ParScale. We are greatful for any users feedback, which should be directed to the [ParScale user forum](http://www.cfdem.com/forum).

###Acknowledgement
Parts of the code were developed in the frame of the [NanoSim project](http://www.sintef.no/Projectweb/NanoSim/) funded by the European Commission through FP7 Grant agreement no. 604656.


