# simulate-openmm

**simulate-openmm** is a Python library designed to simplify and enhance molecular dynamics (MD) simulations using [OpenMM](https://openmm.org/). It provides tools for setting up simulations, integrating Packmol for building initial configurations, and managing simulation workflows efficiently.

## Features

- **Packmol Integration**: Streamlined setup of molecular systems with `mdapackmol`.
- **Simulation Workflows**: Intuitive tools for running and managing MD simulations with OpenMM.
- **Custom Reporters**: Modular reporters for tracking simulation progress and results.
- **Utilities**: Additional tools for parsing inputs and handling common tasks in MD simulations.

## Directory Overview

- `mdapackmol`: Integrates Packmol for molecular system building.
  - Core logic in `mdapackmol.py`.
- `simulate`: Core functionality for running simulations with OpenMM.
  - **`__main__.py`**: Entry point for the package.
  - **`parse`**: Handles input file parsing.
  - **`reporters`**: Custom reporters for simulations.
  - **`run`**: Manages simulation execution workflows.
  - **`utils`**: Auxiliary tools for various tasks.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
