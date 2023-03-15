docker run -it --rm -v $(pwd):/tmp mambaorg/micromamba:1.3.1 \
   /bin/bash -c "micromamba create --yes --name new_env --file /tmp/env.yml && \
                 micromamba env export --name new_env --explicit > env.lock"
