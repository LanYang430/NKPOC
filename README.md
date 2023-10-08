# NKPOC
This repository contains an example code for simulating CO2 adsorption on NKPOC (Nanostructured Kerogen-Pillared Organic Cage) material at 298K. The simulation is performed using a Hybrid MD/GCMC (Molecular Dynamics/Grand Canonical Monte Carlo) sampling approach.

## Files and Directory Structure

- `run_gcmc.py`: This file contains the implementation of the Hybrid MD/GCMC sampling mechanism. It is the main script used to run the simulations.

- `cage.full.itp`,`cage_mol.itp`,`cage_nb_redefined.itp`,: Input file containing parameters for the NKPOC cages.

- `co2_parameters.txt`: Input file containing parameters for CO2 molecules.

- `structure.xyz`: Input file specifying the structure of the NKPOC material.

## Usage

To run the CO2 adsorption simulation, follow these steps:

1. Ensure you have Python and the required libraries installed.

2. Modify the `cage_parameters.txt`, `co2_parameters.txt`, and `structure.xyz` files to specify your system's parameters.

3. Open a terminal and navigate to the project directory.

4. Run the simulation using the following command:
   ```bash
   python run_gcmc.py
