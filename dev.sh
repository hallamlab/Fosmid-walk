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
        singularity build $2/$NAME.sif docker-daemon://$DOCKER_IMAGE:latest
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
    -t)
        echo "x"
    ;;
    *)
        echo "bad option"
        echo $1
    ;;
esac
