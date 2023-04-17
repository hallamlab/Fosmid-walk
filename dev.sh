NAME=fosmid-walk
DOCKER_IMAGE=quay.io/hallam_lab/$NAME
echo image: $DOCKER_IMAGE
echo ""

HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

case $1 in
    --get-usearch)
        # pcc1 vector backbone seq. from https://novoprolabs.com/vector/Vgm4tanq 

        cd docker/lib
        file=usearch.gz
        wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz -O $file
        gunzip $file
        chmod +x *
    ;;
    --build|-b)
        # change the url in python if not txyliu
        # build the docker container locally *with the cog db* (see above)
        cp -r dist/* docker/load
        cd docker
        docker build -t $DOCKER_IMAGE .
    ;;
    --push|-p)
        # login and push image to quay.io, remember to change the python constants in src/
        # sudo docker login quay.io
	    docker push $DOCKER_IMAGE:latest
    ;;
    --sif)
        # test build singularity
        mkdir -p $HERE/scratch
        singularity build $HERE/scratch/$NAME.sif docker-daemon://$DOCKER_IMAGE:latest
    ;;
    --run|-r)
        # test run docker image
        docker run -it --rm \
            -e XDG_CACHE_HOME="/ws"\
            --mount type=bind,source="$HERE/docker/load/",target="/app" \
            --workdir="/ws" \
            -u $(id -u):$(id -g) \
            $DOCKER_IMAGE \
            /bin/bash
    ;;
    --pip-setup)
        # make an environment before hand
        # in that env, install these build tools
        pip install build
    ;;
    --pip-build|-l)
        # build the packge for container
        # this is NOT suitable for pypi because of poorly handled dependencies
        rm -r build && rm -r dist
        python -m build --wheel
    ;;
    --pip-remove|-x)
        pip uninstall $NAME -y
    ;;
    -t)
        shift
        # test run docker image
        docker run -it --rm \
            -e XDG_CACHE_HOME="/ws"\
            --mount type=bind,source="$HERE/src/fosmid_walk",target="/opt/conda/envs/main/lib/python3.11/site-packages/fosmid_walk" \
            --mount type=bind,source="$HERE/scratch/cache",target="/ws" \
            --workdir="/ws" \
            -u $(id -u):$(id -g) \
            $DOCKER_IMAGE \
                foswalk -r /ws/inputs/lake-test.fastq \
                -o /ws/testout \
                -t 12 -b /ws/inputs/bb.fasta
    ;;
    *)
        echo "bad option"
        echo $1
    ;;
esac
