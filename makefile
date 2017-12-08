#以下是目录结构
INC_DIR		:= include#	这里放头文件.h
SRC_DIR		:= src#		这里放源文件.cpp
BUILD_DIR	:= build#	这里放编译过程中的链接文件
BIN_DIR		:= bin#		这里是编译后的文件

#以下是编译最后那个的参数
PETSC_INCLUDE := -I/${PETSC_DIR}/include -I/${PETSC_DIR}/${PETSC_ARCH}/include#	PETSC库所在地
INCLUDE		:= ${PETSC_INCLUDE} -I./${INC_DIR}#									增加头文件的引用目录
include ${PETSC_DIR}/lib/petsc/conf/variables#									PETSC需要
include ${PETSC_DIR}/lib/petsc/conf/rules#										PETSC需要
include ${PETSC_DIR}/lib/petsc/conf/test#										PETSC需要
PETSC_ALL_LIB := ${PETSC_MAT_LIB} ${PETSC_VEC_LIB} ${PETSC_KSP_LIB}#			PETSC需要的其他库也要加在这里

all : ${BIN_DIR}/main
	mpirun -n 8 ./${BIN_DIR}/main
#################
#更改运行参数改这里#
#################

${BIN_DIR}/main : ${BUILD_DIR}/main.o ${BUILD_DIR}/FEI.o
	@mkdir -p ${BIN_DIR}
	@-${CLINKER} $^ -o $@  ${PETSC_ALL_LIB}

${BUILD_DIR}/%.o : ${BUILD_DIR}/%.cpp
	@mkdir -p ${BUILD_DIR}
	${CXX} -c ${INCLUDE} ${CXXFLAGS} $^ -o $@

${BUILD_DIR}/%.cpp : ${SRC_DIR}/%.cpp
	@mkdir -p ${BUILD_DIR}
	python3 ./src/FEI.py $^ $@
	@mkdir -p ${SRC_DIR}/${BUILD_DIR}
	cp -f $@ ${SRC_DIR}/$@