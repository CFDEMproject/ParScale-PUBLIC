![logo](parscale_logo.png)
======
A Compilation of Particle Scale Models.

ParScale is part of the [NanoSim Project](http://sintef.no/NanoSim)


Graz University of Technology and DCS Computing GmbH releases ParScale 1.1.1 beta
------------------
24th of August 2014

###Features
Version 1.1.1-beta is a major upgrade of the ParScale library with extended features to study intra-particle transport phenomena including phase change phenomena (e.g., evaporation of water inside the pores of the particles). Specifically, the following features are new, or have been significantly improved in the latest version of ParScale:

- restructuring of code and documentation that was required for extended ParScale to handle multiple phases
- handling of up to three phases (i.e., solid, liquid, gas), including phase change phenomena. Two phase change models have been added (i.e., "evaporation", as well as "equilibrium sharp"). The volumetric phase change rate can be accounted for in the species transport equations.
- model equations can be used to not only compute concentration profiles, but also update the local phase fraction (e.g., to study the evolution of the local liquid saturation in a porous particle)
- output per model eqn. is now more flexible (it is possible to request more detailed output if required)
- the CVODE integrator can be fully parameterized at runtime (see the settings/integrator.json file)
- the convective flux (caused by phase change phenomena and consumption due to chemical reactions) can now be taken into account in the species transport equation.
- Knudsen diffusion effects have been added 
- new verification cases have been added
- a number of bugs has been fixed

Graz University of Technology and DCS Computing GmbH release ParScale
------------------
31st of December 2014

The Graz University of Technology (TU Graz) together with DCS Computing GmbH (DCS) are please to announced the release of the `1.0.1-beta` version of the tool `ParScale`. This version is distributed is licensed under the [Lesser General Public License](http://www.gnu.org/licenses/lgpl.html) by TU Graz and DCS.

###Features
Version 1.0.1-beta is the first public release of the ParScale library and is meant to introduce its features to a wide audience of users and possible developers of this library. Specifically, the current version of ParScale is able to

- simulate heat and mass transport in one or multiple porous particles.
- perform calculations on an equi-distant spherical one-dimensional coordinate system, or using a simplified model for reactive transport (i.e., the Shrinking Core model)
- predict the outcome of a single non-isothermal irreversible (homogeneous or heterogeneous) chemical reaction (full coupling of temperature and species concentration fields) involving an arbitrary number of solid and gas-phase species.
- consider various "grain-scale models" to picture the formation of a product layer on solids present in the porous of a particle.
- couple with a JSON file (i.e., read data from a JSON file to simulate transient boundary conditions or model parameters), as well as a LIGGGHTS simulation.

A number of test and verification cases, as well as documentation is supplied with ParScale

Note that ParScale is still under active development, and the success of ParScale relies to a large extend on contributions from users. Specifically, we plan to release

- modules to picture multiple reversible reactions that are specified by the user via CHEMKIN-II files,
- an improved version that uses orthogonal collocation methods for more efficient spatial discretization,
- refined sub-models for multi-component diffusion and temperature-dependent transport properties.

The developers of ParScale are thankful for any comments, ideas, or contributions to currently available or future features of ParScale that are posted in the ParScale user forum](http://www.cfdem.com/forum).

###Important Notes
'beta' means, that this version of ParScale has been tested internally to some extend, and now seeks for testing and application by external groups (i.e., non TUG or DCS). Thus, while TUG and DCS have carefully tested ParScale, it is likely that users might find significant bugs, run into unexpected errors, or might miss information on how to use ParScale. We are greatful for any users feedback, which should be directed to the [ParScale user forum](http://www.cfdem.com/forum).

###Acknowledgement
Parts of the code were developed in the frame of the [NanoSim project](http://www.sintef.no/Projectweb/NanoSim/) funded by the European Commission through FP7 Grant agreement no. 604656.
