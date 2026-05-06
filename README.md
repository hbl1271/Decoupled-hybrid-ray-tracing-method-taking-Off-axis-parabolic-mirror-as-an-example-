Here is the professionally drafted README documentation for your project, reflecting the progressive logic from foundational verification to complex surface analysis.

DHRT Off-Axis Reflective System Simulation & Validation Framework
This repository provides a comprehensive MATLAB-based simulation suite for Off-Axis Parabolic (OAP) systems using the Discrete Hybrid Ray Tracing (DHRT) algorithm. 
The project is structured as a three-stage progressive study, moving from fundamental algorithm benchmarking to high-fidelity manufacturing error analysis, 
with rigorous cross-validation against industry-standard Zemax software.

🏗 Progressive Research Logic
The simulation framework is divided into three distinct scenarios. Each scenario consists of a Core Calculation Engine and a corresponding Visualization Script.

Phase 1: Baseline & Nominal System Verification
Engine (Scenario1_Calculation.m): Establishes an ideal OAP system. 
It calculates wavefront errors using the DHRT algorithm and compares them against Analytical and standard Geometric methods.

Visualization (Scenario1_Plotting.m): Focuses on establishing the precision baseline and confirming the numerical stability of the algorithm under ideal conditions.

Phase 2: Rigid-Body Displacements Sensitivity Analysis
Engine (Scenario2_Calculation.m): Introduces assembly-level errors, including lateral/longitudinal decenters and rotational tilts.

Visualization (Scenario2_Plotting.m): Validates the algorithm's robustness when the system's geometric topology is perturbed. 
Results are benchmarked against Zemax to prove the algorithm's reliability for sensitivity analysis.

Phase 3: Surface Deformation 
Engine (Scenario3_Calculation.m): Simulates realistic manufacturing defects. 
It integrates low-frequency Zernike figure errors and mid-to-high frequency sinusoidal modulations (S1 & S2) directly onto the OAP surface.

Visualization (Scenario3_Plotting.m): Analyzes complex wavefront degradation, Zernike term decomposition, and Ray Fan distortions. 
This stage demonstrates that DHRT maintains micron-level consistency with Zemax even under complex micro-perturbations.

Main Simulation Scripts:
Numerical_verification_of_the_nominal_optical_system_calculation.m,   Phase 1.1.
Numerical_verification_of_the_nominal_optical_system_visualization.m, Phase 1.2.
Robustness_under_rigid_body_displacements_calculation.m,     Phase 2.1.
Robustness_under_rigid_body_displacements_visualization.m,   Phase 2.2.
Computational_fidelity_for_complex_surface_deformations_calculation.m,     Phase 3.1.
Computational_fidelity_for_complex_surface_deformations_visualization.m,   Phase 3.2.

Core Function Library
calculate_reflection_DHRT.m: DHRT method for OAP
generate_oap_points.m: Generates OAP point clouds with integrated misalignment or surface sag(This point data is the same as FEA data).
fit_oap.m: Reconstructs geometric parameters (Vertex/Focus) from sampled point clouds.
matrix2zemax.m: Converts MATLAB matrices into .dat Grid Sag files for Zemax interoperability.
