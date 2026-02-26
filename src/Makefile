# Makefile for surge
VERSION=2_0

# Define CC to be the name of the C compiler, and CCOPT to be
#  switches for good optimization in your environment.
CC=gcc
CCOPT=-O3 -march=native

# Define NAUTY to be the directory where files like nauty.h are
#  located. Define NAUTYLIB to be the directory where files like
#  nautyL1.a are located.
NAUTY=${HOME}/nauty
NAUTYLIB=${NAUTY}

# If you have zlib installed, define ZLIB to be "-DZLIB" and LZ
#  to be the linker command for zlib.
ZLIB=-DZLIB
LZ=-lz
# If you don't have zlib installed, make both of them null.
#ZLIB=
#LZ=

SURGE=-g -I ${NAUTY} ${CCOPT} -DWORDSIZE=64 -DMAXN=WORDSIZE \
  -DOUTPROC=surgeproc -DPREPRUNE=surgepreprune \
  -DPRUNE=surgeprune -DGENG_MAIN=geng_main

all : surge

surge : surge.c
	${CC} -o surge ${SURGE} ${ZLIB} \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyL1.a ${LZ}

samples : surge_0 surge_1 surge_2 surge_3

surge_0 : surge.c plugin0.c
	${CC} -o surge_0 ${SURGE} ${ZLIB} -DSURGEPLUGIN='"plugin0.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyL1.a ${LZ}

surge_1 : surge.c plugin1.c
	${CC} -o surge_1 ${SURGE} ${ZLIB} -DSURGEPLUGIN='"plugin1.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyL1.a ${LZ}

surge_2 : surge.c plugin2.c
	${CC} -o surge_2 ${SURGE} ${ZLIB} -DSURGEPLUGIN='"plugin2.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyL1.a ${LZ}

surge_3 : surge.c plugin3.c
	${CC} -o surge_3 ${SURGE} ${ZLIB} -DSURGEPLUGIN='"plugin3.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyL1.a ${LZ}

canonsdf : canonsdf.c
	${CC} -o canonsdf -I${NAUTY} ${CCOPT} canonsdf.c ${NAUTYLIB}/nautyL1.a

clean :
	rm -rf surge surge_0 surge_1 surge_2 surge_3 canonsdf

tarfile :
	rm -rf surge${VERSION}
	mkdir surge${VERSION}
	cp surge.c planarity.c geng.c makefile LICENSE-2.0 README.txt \
	  canonsdf.c plugin[0-3].c surge${VERSION}.pdf surge${VERSION}/
	touch surge${VERSION}/This_is_surge_${VERSION}.txt
	tar cvf surge${VERSION}.tar surge${VERSION}/
