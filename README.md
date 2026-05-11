# Molecule-Bond-Orders
**An exploration of quantum chemistry bond order formulations.**

---

This work presents a computational exploration of bond order calculations through density matrices using a custom C++ implementation of an SCF-based quantum chemistry method.  The program was developed to generate overlap and density matrices and to evaluate chemical bonding using Mulliken, Wiberg, and Mayer bond order formulations. By comparing these methods across molecular systems, the study investigates how differing mathematical treatments of electron density and orbital overlap influence calculated bond order values. The project also demonstrates how quantum chemistry methods can be implemented computationally to provide insight into molecular structures and chemical bonding behavior.

This repository covers the following bond order formulations:
- Wiberg
- Mayer
- Mulliken

## Usage

```bash
docker build -t molecule_bond_order:latest . #1. Download image

bash interactive.sh                          #2. Enter the interactive session

bash run.sh                                  #3. Builds & runs the program, output goes to project_output/.
                                             #   Includes matrix stdout, matrix hdf5, and stats stdout

bash test.sh                                 #4. Runs pytests for this program
```





