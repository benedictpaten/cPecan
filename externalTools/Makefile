rootPath = ../
include ../include.mk

all: 
	cd lastz-distrib-1.03.54 && make install
	mv ./lastz ${binPath}/cPecanLastz
	mv ./lastz_D ${binPath}/cPecanLastz_D

clean: 
	cd lastz-distrib-1.03.54 && make clean 
	rm -rf ${binPath}/cPecanLastz ${binPath}/cPecanLastz_D
 