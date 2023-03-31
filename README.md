## creat docker image to calculation reorganization energy via ORCA/g16

### 1. creat env.lock from env.yml for Dockfile
```shell
./creat_lock.sh
```

### 2. create docker image
```shell
docker build --tag lambda_force:v0.0.1 .
```

### 3. run docker image
```shell
docker run lambda_force:v0.0.1
```