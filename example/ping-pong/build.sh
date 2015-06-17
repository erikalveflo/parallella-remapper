#!/bin/bash

set -e

ESDK=${EPIPHANY_HOME}
ELIBS=${ESDK}/tools/host/lib
EINCS=${ESDK}/tools/host/include
ELDF=${ESDK}/bsps/current/internal.ldf

DEVICE_PROGRAMS="e_main"

SRC="../../src/"

echo -e "==================================================================="

COMMAND="mkdir -p debug"
echo ${COMMAND}
eval ${COMMAND}

echo -e "\nBuild HOST application"
COMMAND="gcc main.c ${SRC}remapping.c ${SRC}pcg.c -o debug/main.elf -I ../../include/ -O2 -I ${EINCS} -L ${ELIBS} -lm -le-hal -le-loader -lpthread -Wno-unused-function -Wall"
echo ${COMMAND}
eval ${COMMAND}

echo -e "\nBuild DEVICE programs"
for PROGRAM in ${DEVICE_PROGRAMS} ; do
	IN=${PROGRAM}
	OUT=debug/${PROGRAM}
	echo "${PROGRAM}"

	echo -e "\nCompiling"
	COMMAND="e-gcc -c ${IN}.c -o ${OUT}.o -O3 -std=c99"
	echo ${COMMAND}
	eval ${COMMAND}

	echo -e "\nLinking"
	COMMAND="e-gcc ${OUT}.o -o ${OUT}.elf -O3 -std=c99 -T ${ELDF} -le-lib -Wall -Werror"
	echo ${COMMAND}
	eval ${COMMAND}

	echo -e "\nConvert ELF to SREC file"
	COMMAND="e-objcopy --srec-forceS3 --output-target srec ${OUT}.elf ${OUT}.srec"
	echo ${COMMAND}
	eval ${COMMAND}
done

echo -e "\nBuild complete\n"
