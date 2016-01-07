rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

cPecanDependencies =  ${basicLibsDependencies}
cPecanLibs = ${basicLibs}

all : ${libPath}/cPecanLib.a ${binPath}/cPecanLibTests ${binPath}/vanillaAlign ${binPath}/trainModels ${binPath}/signalAlign ${sonLibrootPath}/nanoporelib.py
	# disabled right now so that we don't build Lastz every time I do an update
	#cd externalTools && make all
	
clean : 
	rm -f ${binPath}/cPecanRealign ${binPath}/cPecanEm ${binPath}/cPecanLibTests  ${libPath}/cPecanLib.a
	cd externalTools && make clean
	
test : all
	python allTests.py

${binPath}/cPecanRealign : cPecanRealign.c ${libPath}/cPecanLib.a ${cPecanDependencies} 
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cPecanRealign cPecanRealign.c ${libPath}/cPecanLib.a ${cPecanLibs}

${binPath}/vanillaAlign : vanillaAlign.c ${libPath}/cPecanLib.a ${cPecanDependencies} 
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/vanillaAlign vanillaAlign.c ${libPath}/cPecanLib.a ${cPecanLibs}

${binPath}/trainModels : ${rootPath}scripts/trainModels.py
	cp ${rootPath}scripts/trainModels.py ${binPath}/trainModels
	chmod +x ${binPath}/trainModels
	
${binPath}/signalAlign : ${rootPath}scripts/signalAlign.py
	cp ${rootPath}scripts/signalAlign.py ${binPath}/signalAlign
	chmod +x ${binPath}/signalAlign

${sonLibrootPath}/nanoporelib.py : ${rootPath}scripts/nanoporeLib.py
	cp ${rootPath}scripts/nanoporeLib.py ${sonLibRootPath}/nanoporeLib.py

${binPath}/cPecanEm : cPecanEm.py
	cp cPecanEm.py ${binPath}/cPecanEm
	chmod +x ${binPath}/cPecanEm

${binPath}/cPecanLibTests : ${libTests} tests/*.h ${libPath}/cPecanLib.a ${cPecanDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cPecanLibTests ${libTests} ${libPath}/cPecanLib.a ${cPecanLibs}
	
${libPath}/cPecanLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources} 
	ar rc cPecanLib.a *.o
	ranlib cPecanLib.a 
	rm *.o
	mv cPecanLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/
