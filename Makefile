CC=g++
CFLAGS=-c -Wall -O2 -std=c++11
LDFLAGS= -std=c++11

vpath %.cpp src

SOURCES=Converter.cpp\
		Sequence.cpp\
		global.cpp\
		SVMtrain.cpp\
		SequenceNames.cpp\
		KTree.cpp\
		TreeLeafData.cpp


OBJECTS=$(patsubst %.cpp,obj/%.o,$(SOURCES)) 

SRC_GLKKERNEL=mainGLKkernel.cpp
OBJ_GLKKERNEL=obj/mainGLKkernel.o
GLKKERNEL=glk_kernel

SRC_SVMTRAIN=mainSVMtrain.cpp
OBJ_SVMTRAIN=obj/mainSVMtrain.o
SVMTRAIN=glk_train

SRC_SVMCLASSIFY=mainSVMclassify.cpp
OBJ_SVMCLASSIFY=obj/mainSVMclassify.o
SVMCLASSIFY=glk_classify

all: $(SOURCES) $(GLKKERNEL) $(SVMTRAIN) $(SVMCLASSIFY)

clean:
	rm -f $(OBJECTS) $(OBJ_GLKKERNEL) $(OBJ_SVMTRAIN) $(OBJ_SVMCLASSIFY) $(GLKKERNEL) $(SVMTRAIN) $(SVMCLASSIFY)

obj/%.o : %.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

$(GLKKERNEL): $(OBJECTS) $(OBJ_GLKKERNEL)
	$(CC) $(LDFLAGS) $(OBJECTS) $(OBJ_GLKKERNEL) -o $@

$(SVMTRAIN): $(OBJECTS) $(OBJ_SVMTRAIN)
	$(CC) $(LDFLAGS) $(OBJECTS) $(OBJ_SVMTRAIN) -o $@

$(SVMCLASSIFY): $(OBJECTS) $(OBJ_SVMCLASSIFY)
	$(CC) $(LDFLAGS) $(OBJECTS) $(OBJ_SVMCLASSIFY) -o $@

install:
	cp glk_kernel glk_train glk_classify /bin
