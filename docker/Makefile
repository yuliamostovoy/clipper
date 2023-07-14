VERSION = 1.1
TAG1 = quay.io/ymostovoy/clipper:$(VERSION)
TAG2 = quay.io/ymostovoy/clipper:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
