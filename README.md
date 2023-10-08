# NKPOC
This repository contains an example code for simulating CO2 adsorption on NKPOC (Nanostructured Kerogen-Pillared Organic Cage) material at 298K. The simulation is performed using a Hybrid MD/GCMC (Molecular Dynamics/Grand Canonical Monte Carlo) sampling approach.

## Files and Directory Structure

- `run_gcmc.py`: This file contains the implementation of the Hybrid MD/GCMC sampling mechanism. It is the main script used to run the simulations.

- `cage.full.itp`, `cage_mol.itp`, `cage_nb_redefined.itp`: Input files containing parameters for the NKPOC structure, individual cage molecules, and redefined inter-cage interactions parameters, respectively.

- `co2.itp`: Input file containing parameters for CO2 molecules.

## Details of Simulation Parameters

The `run_gcmc.py` script allows you to customize various simulation parameters to control the behavior of the CO2 adsorption simulation. Here are the details of these parameters:

- `T = 298 * kelvin`: Sets the temperature of the simulation to 298 Kelvin.

- `dt = 1.0`: Specifies the time step for the simulation in femtoseconds (fs). In this case, it is set to 1.0 fs.

- `n_cycles = 10000000`: Sets the total number of simulation steps. This parameter determines the duration of the entire simulation.

- `n = 0`: Specifies the initial number of adsorbed molecules. The simulation starts with 0 adsorbed molecules and will incrementally add or remove them.

- `dn = 1.0`: Specifies the number of molecules to insert or delete in each GCMC move. In this case, it is set to 1.0, meaning one molecule is added or removed at a time.

- `P = 0.1 * bar`: Sets the CO2 pressure in the simulation to 0.1 bar.

- `P_ext = 1.0 * bar`: Specifies the external pressure applied to the system.

- `freq_mc = 50`: Determines the frequency of performing GCMC sampling within the MD simulation. In this case, GCMC sampling is performed every 50 MD steps.

You can adjust these parameters according to your specific simulation requirements by modifying the `run_gcmc.py` script.


## Usage

To run the CO2 adsorption simulation, follow these steps:

1. Ensure you have the GPU version of OpenMM installed.

2. Run the simulation using the following command:
   ```bash
   python run_gcmc.py
