1D Non-Catalytic Heterogeneous Reaction Validation
=======================================================

This validation case validates the implementation of a single reaction according Wen, "Noncatalytic Heterogenous Solid Fluid Reaction Models", p. 42-44 and offers a comparison to the presented analytical solution for Thiele Modulus of ~ 3.16.

Parameters
---------------------

This test case is based on the seetings for copper oxidation presented by Noorman et al. (CEJ 167:297-307). Results fo this test case are shown in Figure 2 and 3 of this paper.
Key assumptions are
- isothermal case
- 1-order with respect to gas phase (i.e., "speciesA")
- 1 Cu + 1/2 O2 -> 1 CuO

How to run:
-----------

- Open a terminal at the `in.file_*`´s postion and just type `./Allrun_convective`
- In case `Allrun_*` and `Allclean` are not executable type `chmod +x All*`
- After the simulation is done type `Allclean` to delete all generated files
