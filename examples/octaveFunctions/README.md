Repository with octave post-processing functions
===================

Usage
------------

Users should add functions needed to postprocess LIGGGHTS simulations here. The user should then include the following lines to the octave-style post processing script:

```
PARSCALE_SRC_DIR = getenv('PASCAL_SRC_DIR');
if(isempty(PARSCALE_SRC_DIR))
    error('The user has not set the environment variable "PASCAL_SRC_DIR". Please do so, e.g., in your .bashrc file, to use this octave script.')
end
addpath([PARSCALE_SRC_DIR,'/../examples/octaveFunctions']);

```

 
