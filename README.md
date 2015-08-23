![logo](parscale_logo.png)
======
A Compilation of Particle Scale Models.

ParScale is part of the [NanoSim Project](http://sintef.no/NanoSim)


Copyright Notice
------------------

- Copyright 2014 - Graz University of Technology (S. Radl, T. Forgber).
- Copyright 2014 - DCS Computing GmbH, Linz (C. Goniva, C. Kloss).

All rights reserved.

Note that some modules within the ParScale package are not under the copyright of the above copyright holders. These modules are included in the ParScale package with a valid open-source license. Please see the text in the source files, or the corresponding sub-directory, for more detailed information on copyright and licensing issues.

License
-----------------
See the [LICENSE.md](LICENSE.md) file for details.

Warranty
-----------------
ParScale is distributed in the hope that it will be applied and further developed by researchers and engineerings in academia or industry. However, ParScale comes without any warranty, without even the implied warranty of merchantability or fitness for a particular purpose. 

Scope
---------------------------------------
ParScale stands for “Particle Scale Models”, and is a publicly available library to simulate intra-particle transport phenomena, e.g., heat or mass transfer, as well as homogeneous and heterogeneous reactions in porous particles. For example, the rate of oxidation of porous metal particles suspended by a hot gas stream can be predicted by ParScale.

ParScale is available via https://github.com/CFDEMproject, and is designed as an add-on to the soft-sphere particle simulation tool `LIGGGHTS`. As such, ParScale is part of the [CFDEMproject](http://www.cfdem.com) initiated by Christoph Kloss and Christoph Goniva. However, ParScale can also be run in stand-alone mode for predicting, e.g., the chemical conversion of particles, fitting of model parameters to experimental data, or optimization studies of reactive particulate systems.

Important Hints for Usage
-----------------

- for installation instructions see the file src/INSTALL.md
- there are separate README.md files in each sub-directory, so you might want to check them
- for a detailed model and user documentation see the doc folder
- there is background information for developers in src/DEVELOPERS.md
- in case you plan to add more details to the documentation, Markdown should be used for this purpose. The tool "ReText" (in the OpenSUSE standard package) or gedit (with the appropriate [gedit-markdown plugin](http://www.jpfleury.net/en/software/gedit-markdown.php)) for previewing markdown documentation on linux systems should be used to edit and view documentation.


Credits
-------------------
The following persons have made significant contributions to ParScale (in alphabetic order)

- Andreas Aigner (DCS)
- Thomas Forgber (TU Graz)
- Christoph Kloss (DCS)
- Stefan Radl (TU Graz)
