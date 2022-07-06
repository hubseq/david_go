NAME=david_go
TAG=00.00.00

all: build push

build:
	docker build -t $(NAME):$(TAG) -f Dockerfile ../
	docker tag $(NAME):$(TAG) $(REGISTRY):$(TAG)

push:
	docker push $(REGISTRY):$(TAG)
