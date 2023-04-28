## creat docker image to calculation  via ORCA/g16
### 1. creat env.lock from env.yml for Dockfile
```shell
./creat_lock.sh
```

### 2. create docker image
```shell
docker build --tag lambda_force:vXXX .
```

### 3. run docker image
```shell
docker run lambda_force:vXXX
```

## What can be obtained (calculator 'provides')
### reorganization energy via 4-point approach: 'lambda-4p'
### reorganization energy via displacement approach: 'lambda-disp' (J. Chem. Phys., 2001, 115, 9103.)
### reorganization energy via force approach: 'lambda-force' (https://doi.org/10.1021/acs.jpclett.3c00749)
### bond order difference: 'BOD'
### spin-orbit coupling: 'SOC' (only for ORCA)
### nonadiabatic coupling: 'NAC'