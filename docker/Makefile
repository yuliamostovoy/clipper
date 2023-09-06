VERSION = 1.1
TAG1 = quay.io/ymostovoy/bigclipper:$(VERSION)
TAG2 = quay.io/ymostovoy/bigclipper:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
