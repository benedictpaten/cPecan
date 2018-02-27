# Definitions
build_tool = runtime-container.DONE
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = tpesout/cpecan
tag = ${git_commit}

# Steps
build:
	# so that our docker tag and build match the current git hash and repo
	-rm -rf ../../tmp.cPecan
	-rm -rf tmp.cPecan
	mkdir ../../tmp.cPecan
	-rm -rf tmp.cPecan
	cp -r ../* ../../tmp.cPecan
	-cd ../../tmp.cPecan && make clean
	mv ../../tmp.cPecan tmp.cPecan
	# do the docker build
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest
	rm -rf tmp.cPecan
	touch ${build_tool}

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest

clean:
	-rm ${build_tool}
