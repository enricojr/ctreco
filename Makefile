currentDirectory=`pwd`
libraryDirectory=libraries
OBJS = ${libraryDirectory}/CTProjectorClass.o	${libraryDirectory}/frameClass.o ${libraryDirectory}/FourierClass.o
EXEC = makeSinogram sumFramesWithFormat normalizeToOB interpolateBadPixels filterImage reconstruct_FBP reconstruct_OSEM
EXTRACXXFLAGS=-O2 -W -Wall

opt:    EXTRACXXFLAGS=-O2 -Wall

dbg:    EXTRACXXFLAGS=-g3 -Wall

opt:    all

dbg:    all

all:	info makeSinogram sumFramesWithFormat normalizeToOB interpolateBadPixels filterImage reconstruct_FBP reconstruct_OSEM

clean:
	@rm -f $(OBJS) $(EXEC) *.so

info:
	@echo " --------------------------------------------------------"
	@echo " - Building ctreco"
	@echo " - 	currentDirectory=${currentDirectory}"
	@echo " - 	libraryDirectory=${libraryDirectory}"
	@echo " - make dbg   provides a debug build"
	@echo " --------------------------------------------------------"

${libraryDirectory}/CTProjectorClass.o:	${libraryDirectory}/CTProjectorClass.cc
	@echo "library CTProjectorClass ..."
	g++ ${EXTRACXXFLAGS} -std=c++0x -c ${libraryDirectory}/CTProjectorClass.cc
	mv CTProjectorClass.o ${libraryDirectory}

${libraryDirectory}/frameClass.o:	${libraryDirectory}/frameClass.cc
	@echo "library frameClass ..."
	g++ ${EXTRACXXFLAGS} -c ${libraryDirectory}/frameClass.cc
	mv frameClass.o ${libraryDirectory}

${libraryDirectory}/FourierClass.o:	${libraryDirectory}/FourierClass.cc
	@echo "library FourierClas ..."
	g++ ${EXTRACXXFLAGS} -c ${libraryDirectory}/FourierClass.cc -lfftw3 -lm
	mv FourierClass.o ${libraryDirectory}

makeSinogram:	${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o ${libraryDirectory}/FourierClass.o makeSinogram.cc
	@echo "makeSinogram ...."
	g++ ${EXTRACXXFLAGS} -o makeSinogram ${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o makeSinogram.cc -lfftw3 -lm

sumFramesWithFormat:	${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o sumFramesWithFormat.cc
	g++ ${EXTRACXXFLAGS} -o sumFramesWithFormat ${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o sumFramesWithFormat.cc -lfftw3 -lm

normalizeToOB:	${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o normalizeToOB.cc
	g++ ${EXTRACXXFLAGS} -o normalizeToOB ${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o normalizeToOB.cc -lfftw3 -lm

interpolateBadPixels:	${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o interpolateBadPixels.cc
	g++ ${EXTRACXXFLAGS} -o interpolateBadPixels ${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o interpolateBadPixels.cc -lfftw3 -lm

filterImage:	${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o filterImage.cc
	g++ ${EXTRACXXFLAGS} -o filterImage ${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o filterImage.cc -lfftw3 -lm

reconstruct_FBP:	${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o reconstruct_FBP.cc
	g++ ${EXTRACXXFLAGS} -o reconstruct_FBP ${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o reconstruct_FBP.cc -lfftw3 -lm

reconstruct_OSEM:	reconstruct_OSEM.o ${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o reconstruct_OSEM.cc
	g++ ${EXTRACXXFLAGS} -pthread -o reconstruct_OSEM ${libraryDirectory}/CTProjectorClass.o ${libraryDirectory}/frameClass.o  ${libraryDirectory}/FourierClass.o reconstruct_OSEM.o -lfftw3 -lm

reconstruct_OSEM.o:	reconstruct_OSEM.cc
	g++ ${EXTRACXXFLAGS} -std=c++0x -pthread -c reconstruct_OSEM.cc


