# Master thesis

## Modeling of Blood Flow in Zebrafish With Cancer

For decades, the improvement of cancer research has relied on in vivo representations
for analyzing cancer development and treatment. Particularly, one of
the main reasons why zebrafishes are suitable for studying a wide variety of cancer
types is the dynamic visualization of tumor growth in vivo. We can observe
tumor cells spreading locally into nearby healthy tissues, or globally through
lymph and blood vessels. In the blood vessels, abnormally high permeability
of the vessel walls can result in blood eruptions, consequently causing leakage
of tumor cells into the surrounding tissue. In this thesis, we aim to model
this phenomenon by means of partial differential equations. Different governing
equations for flow modeling are proposed for different domains. Flow in the
blood vessel (viscous domain) is modeled by simplified Navier-Stokes equations,
while Darcy’s law is representing the flow in the tissue (porous domain). Together,
coupled Darcy-Stokes equations will be used as an approximation of the
flow within and between the two domains. The advantage of this modeling is
the possibility to observe eruptions that emerge at the interface (endothelial
cells) between the domains. The results generated in this thesis illustrate the
evolution of the eruption from various values of permeability in a simplified geometry.
The accuracy of the simulated results have been verified by comparing
them with available experimental data.

### Tools used in the master thesis

* Python2-3
* FEniCS
  * Version 2019.1.0
* Gmsh
  * Used for generation of the mesh.
* ParaView
  * Used for visualization of the simulations.
* Docker 
  * FEniCS: prebuilt, high-performance Docker images)
* High Performance Computing
  * Sigma2
* Matlab
* ImageJ
  * Used for scaling.
* THIS PROJECT FOLDER IS UNDER CONSTRUCTION
