# Examples

This directory contains a large suite of examples demonstrating most
of the features available in PyLith. The examples progress from vary
simple problems to more complex ones similar to those encountered in
research. The examples are broken up into the small sets listed below,
in which all of the simulations within a set use the same geometry.

All of the examples use HDF5 output. PyLith supports writing output as
VTK files, but this does not scale to large research level problems.
Additionally, PyLith supports HDF5 output with external data sets, and
this is the output recommended for large parallel runs.

## 2d/box

Simple axial and shear loading in a 2-D box with Dirichlet and Neumann
boundary conditions and isotropic, linear elasticty.

## 3d/box

Same simulations as those in 2d/box, but in 3-D.

## 2d/strike-slip

## 2d/reverse

## 2d/subduction

## 2d/magma-chamber

Considering for future release.

## 3d/magma-chamber

Considering for future release.

## 3d/strike-slip

Planned for future release.

## 3d/subduction

## bar_shearwave

Not yet updated for v3.0.


Add example of using the initial conditions.



3d/box
  tet mesh (default)
  hex mesh
  two materials (upper & lower crust)

  Step01: axial extension w/Dirichlet BC
    Dirichlet BC (UniformDB)
    IsotropicLinearElasticity (2x UniformDB)
  Step02: simple shear w/Dirichlet BC
    Dirichlet BC (SimpleDB)
    IsotropicLinearElasticity (1x SimpleDB)
  Step03: simple shear w/Dirichlet + Neumann BC
    Dirichlet BC (SimpleDB)
    Neumann BC (UniformDB)
    IsotropicLinearElasticity (1x SimpleDB)
  Step03: time-dependent shear w/Dirichlet and Neumenn BC
    Dirichlet BC (UniformDB)
    Neumann BC (UniformDB)
    IsotropicLinearElasticity (1x SimpleDB)
  Exercises
    Change to hex mesh
    Change material properties
    Change BC to give axial compression in the +y direction
    Change basis and quadrature order


2d/strike-slip w/throughgoing fault
  tri mesh
  quad mesh (default)
  two materials (25% contrast)

  Step01: static w/prescribed slip; fixed boundaries
  Step02: quasistatic: boundaries move; fault catches up w/prescribed slip
  Step03: quasistatic: multiple earthquake cycles w/presscribed slip

  Exercises
    Change to tri mesh

2d/reverse w/buried fault + splay fault
  tri mesh (default)
  quad mesh
  three materials: upper crust, wedge, lower crust

  Step01: gravitational body forces (no reference state)
  Step02: gravitational body forces (reference state)
    Reference state (SimpleDB)
  Step03: gravitational body forces + incompressible elasticity
  Step05: distributed surface load
    Neumann BC (SimpleGridDB)
  Step06: coseismic slip (main fault)
  Step07: coseismic slip (main + splay fault)
  Step08: coseismic slip + viscoelastic relaxation w/linear Maxwell
  Step09: coseismic slip + viscoelastic relaxation w/powerlaw

  Exercises
    Change to quad mesh
    Compare basis order 1 and 2 for Steps 1 and 2
    Change material properties

3d/strike-slip w/throughgoing fault :LATER:?
  match 2d/strike-slip w/throughgoing fault

2d/subduction :UPDATE:

2d/magmachamber :LATER:

3d/subduction :UPDATE:

3d/strike-slip :LATER:

debugging :UPDATE:
  quadrature orders don't match
  typos in component names

meshing [no changes]



Concepts:
  * Output
    + domain
      - solution
    + boundary
      - solution
    + fault
    + material
      - auxiliary subfields
      - derived field
      - solution
    + points
      - solution
    + Dirichlet BC
      - auxiliary field
    + Neumann BC
      - auxiliary field
    + Add field to HDF5 and update Xdmf
    + Projection of field to basis order 0/1
  * Spatial databases
    + Generate spatial database using Python
    + SimpleGridDB
    + UniformDB
  * Simulation
    + Neumann BC
      - side loads
      - surface loads
    + Dirichlet BC
      - roller BC
      - constant rate BC
    + Initial conditions
    + Elasticity
      - Referense state
      - Body forces
      - Gravitational body forces
      - Rheologies
        * IsotropicLinearElasticity
	* IsotropicLinearMaxwell
	* IsotropicPowerlaw
	* IsotropicDruckerPrager
    + Incompressible elasticity
      - gravitational body forces
      - initial conditions
    + Fault
      + Prescribed slip step function, single rupture
      + Prescribed slip step function + creep, multiple ruptures
  * Order of basis functions
    - gravity w/1st order vs 2nd order
