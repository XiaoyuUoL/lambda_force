## Overview
This repository provides
1. Script to call QC calculation via G16/ORCA and then calculate electronic properties
2. Shell script to create a Docker image

## Electronic Properties Calculation (in `opt/` folder)
There are three folders in `opt` folder:
### 1. `opt/binary_folder`
Python codes (`.py` files) and input files (`.yml` files). The calculation system can be defined in the `molecule.yml` file. QC calculation options can be modified in the `calculator.yml` file. (e.g., `Software` for QC package, `ProcNumber` for parallel cores, `Functional` for DFT functional, etc.) Calculated electronic properties can also be modified in the `calculator.yml` file (`provides`), including
(1). `lambda_4p`: reorganization energy/HR factor via 4-point and displacement approaches for S$_1$. (J. Chem. Phys., 2001, 115, 9103.)

(2). `lambda_force`: reorganization energy/Huang-Rhys factor via force approach for S$_1$. (https://doi.org/10.1021/acs.jpclett.3c00749)

(3). `BOD`: bond order difference between closed-shell configuration (usually S$_0$) and HOMO->LUMO excitation configuration (usually S$_1$ state).

(4). `NAC`: nonadiabatic coupling between S$_0$ and S$_1$.

(5). `SOC`: spin-orbit coupling (only ORCA available now) between singlet and triplet states.

Here, (1), (2), and (3) are related to https://doi.org/10.1021/acs.jpclett.3c00749; (4) and (5) are related to [J. Chem. Phys. 2022, 157, 134106.]

### 2. `opt/software_folder`
Thirds party software is not in GitHub, considering the copyright issue. There are only examples of loading scripts (`.sh` files). Users should provide QC packages themselves, and the name of packages should be consistent with loading scripts. For example, Gaussian 16 (`opt/software_folder/G16` and `g16.sh`) and ORCA 5.0.4 (`opt/software_folder/orca_5_0_4_linux_x86-64_shared_openmpi411` and `orca.sh`).

### 3. `opt/result_folder`
An empty where we put the input/output files of QC packages (e.g., `.gjf`/`.log`/`.chk`/`.fchk` files of g16 calculation and `.inp`/`.log` files of the orca calculation) and output files of our python script (`result.yml`, `HR_factor_disp.dat`/`HR_factor_force.dat`, `BOD.dat`, `SOC.dat`, and `NAC.dat` files).


## Docker Image Creation
#### 1. create the `env.lock` file from env.yml for Dockfile
```shell
./creat_lock.sh
```

#### 2. create the docker image
```shell
docker build --tag lambda_force:vXXX .
```

#### 3. run the docker image
```shell
docker run lambda_force:vXXX
```

### Notice: if users want to use the Python scripts instead of creating the docker image, they should modify `.sh` files in `opt/software_folder` (value of `SOFTWARE`) and then run
```shell
cd opt/binary_folder
python main.py
```
