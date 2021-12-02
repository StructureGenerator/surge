# Makefile for surge
VERSION=1_0

# Define CC to be the name of the C compiler, and CCOPT to be
#  switches for good optimization in your environment.
CC=gcc
CCOPT=-march=native -mtune=native -O3

# Define NAUTY to be the directory where files like nauty.h are
#  located. Define NAUTYLIB to be the directory where files like
#  nautyW1.a are located.
NAUTY=${HOME}/nauty
NAUTYLIB=${NAUTY}

# If you have zlib installed, define ZLIB to be "-DZLIB" and ZLIBLIB
#  to be the linker command for zlib.
ZLIB=-DZLIB
ZLIBLIB=-lz
# If you don't have zlib installed, make both of them null.
#ZLIB=
#ZLIBLIB=

SURGE=-g -I ${NAUTY} -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc ${CCOPT} \
  -DPREPRUNE=surgepreprune ${ZLIB} -DPRUNE=surgeprune -DGENG_MAIN=geng_main

all : surge surge_0 surge_1 surge_2 surge_3

surge : surge.c 
	${CC} -o surge ${SURGE} \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a ${ZLIBLIB}

surge_0 : surge.c plugin0.c
	${CC} -o surge_0 ${SURGE} -DSURGEPLUGIN='"plugin0.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a ${ZLIBLIB}

surge_1 : surge.c plugin1.c
	${CC} -o surge_1 ${SURGE} -DSURGEPLUGIN='"plugin1.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a ${ZLIBLIB}

surge_2 : surge.c plugin2.c
	${CC} -o surge_2 ${SURGE} -DSURGEPLUGIN='"plugin2.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a ${ZLIBLIB}

surge_3 : surge.c plugin3.c
	${CC} -o surge_3 ${SURGE} -DSURGEPLUGIN='"plugin3.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a ${ZLIBLIB}

clean : 
	rm -rf surge surge_0 surge_1 surge_2 surge_3 canonsdf

# This one is unoptimized for easier debugging
surge_g : 
	${CC} -o surge_g ${SURGE} -g -O0 \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a ${ZLIBLIB}

canonsdf : canonsdf.c
	${CC} -o canonsdf -I {NAUTY} canonsdf.c ${NAUTYLIB}/nautyL1.a

tarfile :
	rm -rf surge${VERSION}
	mkdir surge${VERSION}
	cp SDFformats.txt surge.c planarity.c geng.c makefile COPYRIGHT \
	   canonsdf.c plugin[0-3].c surge1_0.pdf surge${VERSION}/
	touch surge${VERSION}/This_is_surge_${VERSION}.txt
	tar cvf surge${VERSION}.tar surge${VERSION}/
