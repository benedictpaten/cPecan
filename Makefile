rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

cPecanDependencies =  ${basicLibsDependencies}
cPecanLibs = ${basicLibs}

all : ${libPath}/cPecanLib.a ${binPath}/cPecanLibTests ${binPath}/cPecanRealign ${binPath}/cPecanEm ${binPath}/cPecanModifyHmm 
	cd externalTools && make all
	
clean : 
	rm -f ${binPath}/cPecanRealign ${binPath}/cPecanEm ${binPath}/cPecanLibTests  ${libPath}/cPecanLib.a
	cd externalTools && make clean
	
test : all
	python allTests.py

${binPath}/cPecanRealign : cPecanRealign.c ${libPath}/cPecanLib.a ${cPecanDependencies} 
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cPecanRealign cPecanRealign.c ${libPath}/cPecanLib.a ${cPecanLibs}
	
${binPath}/cPecanEm : cPecanEm.py
	cp cPecanEm.py ${binPath}/cPecanEm
	chmod +x ${binPath}/cPecanEm
	
${binPath}/cPecanModifyHmm : cPecanModifyHmm.py
	cp cPecanModifyHmm.py ${binPath}/cPecanModifyHmm
	chmod +x ${binPath}/cPecanModifyHmm

${binPath}/cPecanLibTests : ${libTests} tests/*.h ${libPath}/cPecanLib.a ${cPecanDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cPecanLibTests ${libTests} ${libPath}/cPecanLib.a ${cPecanLibs}
	
${libPath}/cPecanLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources} 
	ar rc cPecanLib.a *.o
	ranlib cPecanLib.a 
	rm *.o
	mv cPecanLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/
