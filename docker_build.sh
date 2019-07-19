
# Update the Git and AWS strings in the Dockerfile!!

docker build --no-cache -t c4 -f DockerfileAWS

# docker run --rm -e STATE=va -e SEED=300 -e METHOD=POWER c4

docker tag c4:latest 808035620362.dkr.ecr.us-east-1.amazonaws.com/cluscious:latest

$(aws ecr get-login --no-include-email)
docker push 808035620362.dkr.ecr.us-east-1.amazonaws.com/cluscious:latest

