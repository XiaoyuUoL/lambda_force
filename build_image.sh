#!/bin/bash -e
set -euxo pipefail

source ../project_config.sh

if [[ -z "$NAME" ]]; then
    echo "Please edit project_config.sh and specify a NAME".
    exit 1
fi

VERSION=$(git describe)

if [[ ! "$VERSION" =~ ^$NAME/v* ]]; then
    echo "Last tag did not confirm to naming spec $NAME/v1.1.1"
    exit 1
fi

docker build  --tag $VERSION .
