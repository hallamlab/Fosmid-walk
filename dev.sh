NAME=simple-metagenomics
DOCKER_IMAGE=quay.io/txyliu/$NAME
echo image: $DOCKER_IMAGE
echo ""

HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

case $1 in
    --get-cog)
        # get the cog db is needed for the pipeline and is packaged with the container
        wget -P ./docker/db/ https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz
    ;;
    --build|-b)
        # change the url in python if not txyliu
        # build the docker container locally *with the cog db* (see above)
        cd docker 
        sudo docker build -t $DOCKER_IMAGE .
    ;;
    --push|-p)
        # login and push image to quay.io, remember to change the python constants in src/
        # sudo docker login quay.io
	    sudo docker push $DOCKER_IMAGE:latest
    ;;
    --sif)
        # test build singularity
        sudo singularity build $2/$NAME.sif docker-daemon://$DOCKER_IMAGE:latest
    ;;
    --run|-r)
        # test run docker image
        docker run -it --rm \
            -e XDG_CACHE_HOME="/ws"\
            --mount type=bind,source="$HERE/scratch",target="/ws" \
            --mount type=bind,source="$HERE/scratch/res",target="/ref"\
            --mount type=bind,source="$HERE/scratch/res/.ncbi",target="/.ncbi" \
            --mount type=bind,source="$HERE/docker/load/",target="/app" \
            --workdir="/ws" \
            -u $(id -u):$(id -g) \
            $DOCKER_IMAGE \
            /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && conda run -n main ${*: 2:99}"
    ;;
    --pip-setup)
        # make an environment before hand
        # in that env, install these build tools
        pip install twine build
    ;;
    --pip-install|-n)
        # build and install the package locally
        python setup.py build \
        && python setup.py install
    ;;
    --pip-build|-l)
        # build the packge for upload to pypi
        rm -r build && rm -r dist
        python -m build --wheel
    ;;
    --pip-upload|-u)
        # upload to pypi
        # use testpypi for dev
        # PYPI=testpypi
        PYPI=pypi
        TOKEN=`cat secrets/${PYPI}`
        python -m twine upload --repository $PYPI dist/* -u __token__ -p $TOKEN
    ;;
    --pip-remove|-x)
        pip uninstall simple-metagenomics
    ;;
    -t)
        #
        # scratch space for testing stuff
        #
        cd $HERE/src
        # # python -m simple_meta setup -ref $HERE/scratch/test1/ref -c docker
        # python -m simple_meta run -ref $HERE/scratch/res -i SRR19573024 -o $HERE/scratch/test2/ws -t 16
        # # --mock

        # python -m simple_meta setup -ref $HERE/scratch/res -c singularity
        python -m simple_meta run -r $HERE/scratch/res -i SRR19573024 -o $HERE/scratch/test3/ws -t 14
    ;;
    *)
        echo "bad option"
        echo $1
    ;;
esac
