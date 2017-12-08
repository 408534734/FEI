/*
	如果出现了不能写文件的操作，请在命令行执行以下命令：
	echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
*/

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <cmath>
#include <FEI.h>
#define pi M_PI
#define int PetscInt
#define float PetscScalar

;char help[] = "解一个偏微分方程\n\
-△u	= cos(3x)sin(πy)	(x,y) ∈ G=(0,π)x(0,1)\n\
u(x,0) = u(x,1) = 0, 0 <= x <= 1\n\
u(0,y) =  sin(πy)/(9+π²)\n\
u(π,y) = -sin(πy)/(9+π²)\n\n";

float x_start = 0, x_end = pi;
//x的取值范围
float y_start = 0, y_end = 1;
//y的取值范围
float h = 0.01;//这里设置为0.002需要8核运行1分钟左右，但设置为0.001就可能会在运行完之前就被操作系统击杀[手动捂脸]
//网格的边长为0.01
int x_size = (x_end-x_start)/h+1;
//x轴方向上的网格数
int y_size = (y_end-y_start)/h+1;
//y轴方向上的网格数

float f(float x, float y) {
	return cos((float)3*x)*sin(pi*y);
	float f = cos((float)3*x)*sin(pi*y);
}

enum matType{ AIJ, DENSE };
class EMat {//封装的矩阵
	private:
	PetscErrorCode ierr;
	PetscErrorCode init(matType type, int rows, int columns) {
		this -> globalRows = rows;
		this -> globalColumns = columns;
		ierr = MatCreate(PETSC_COMM_WORLD, &self);CHKERRQ(ierr);
		ierr = MatSetSizes(self, PETSC_DECIDE, PETSC_DECIDE, rows, columns);CHKERRQ(ierr);
		switch (type) {
			case AIJ:
				ierr = MatSetType(self, MATMPIAIJ);CHKERRQ(ierr);
				ierr = MatMPIAIJSetPreallocation(self, 5, NULL, 4, NULL);CHKERRQ(ierr);
				break;
			case DENSE:
				ierr = MatSetType(self, MATMPIDENSE);CHKERRQ(ierr);
				ierr = MatMPIDenseSetPreallocation(self, NULL);CHKERRQ(ierr);
				break;
		}
		ierr = MatGetOwnershipRange(self, &rowStart, &rowEnd);CHKERRQ(ierr);
	}
	public:
	Mat self;
	int globalRows, globalColumns;
	int rowStart, rowEnd;
	EMat(matType type, int rows, int columns) {
		init(type, rows, columns);
	}
	PetscErrorCode setValue(int row, int column, float value) {
		ierr = MatSetValue(self, row, column, value, INSERT_VALUES);CHKERRQ(ierr);
	}
	PetscErrorCode setValues(int row, int n, int *columns, float *values) {
		ierr = MatSetValues(self, 1, &row, n, columns, values, INSERT_VALUES);CHKERRQ(ierr);
	}
	void setValues(FEI &newFEI) {
		setValues(newFEI.row, newFEI.n, newFEI.columns, newFEI.values);
	}
	PetscErrorCode startAsm() {
		ierr = MatAssemblyBegin(self, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	PetscErrorCode endAsm() {
		ierr = MatAssemblyEnd(self, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	PetscErrorCode print() {
		ierr = MatView(self, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	PetscErrorCode setName(char *name) {
		PetscObjectSetName((PetscObject)self, name);
	}
};

class EVec {//封装的向量
	private:
	PetscErrorCode ierr;
	float *array;
	PetscErrorCode init(int size) {
		this->size = size;
		ierr = VecCreate(PETSC_COMM_WORLD, &self);CHKERRQ(ierr);
		ierr = VecSetType(self, VECMPI);CHKERRQ(ierr);
		ierr = VecSetSizes(self, PETSC_DECIDE, size);CHKERRQ(ierr);
		ierr = VecGetOwnershipRange(self, &start, &end);CHKERRQ(ierr);
		ierr = VecSet(self, 0);
	}
	public:
	Vec self;
	int size, start, end;
	EVec(int size) {
		init(size);
	}
	PetscErrorCode setValue(int position, float value) {
		ierr = VecSetValue(self, position, value, INSERT_VALUES);CHKERRQ(ierr);
	}
	void setValue(FEI &newFEI) {
		setValue(newFEI.row, newFEI.rightItem);
	}
	PetscErrorCode startAsm() {
		ierr = VecAssemblyBegin(self);CHKERRQ(ierr);
	}
	PetscErrorCode endAsm() {
		ierr = VecAssemblyEnd(self);CHKERRQ(ierr);
	}
	PetscErrorCode print() {
		ierr = VecView(self, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	PetscErrorCode setName(char *name) {
		PetscObjectSetName((PetscObject)self, name);
	}
	PetscErrorCode startGetValues() {
		ierr = VecGetArray(self, &array);CHKERRQ(ierr);
	}
	PetscErrorCode endGetValues() {
		ierr = VecRestoreArray(self, &array);CHKERRQ(ierr);
	}
	float& operator[](int i) {
		return array[i-start];
	}
};

class EKsp {//封装的KSP解法器
	private:
	PetscErrorCode ierr;
	public:
	EKsp(){}
	int calculate(EMat &A, EVec &x, EVec &b) {
		KSP ksp;
		ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp, A.self, A.self);CHKERRQ(ierr);
		ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT,PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
		ierr = KSPSolve(ksp, b.self, x.self);CHKERRQ(ierr);
		int iterationTimes;
		ierr = KSPGetIterationNumber(ksp, &iterationTimes);CHKERRQ(ierr);
		return iterationTimes;
	}
};

class MatlabPrinter {//封装的用于输出到MATLAB（ASCII形式）的打印机
	private:
	PetscViewer printer;
	void init(PetscObject obj) {
		const char *obj_name;
		char *full_name;
		PetscObjectGetName(obj, &obj_name);
		full_name = (char *)malloc(sizeof(char)*(strlen(obj_name)+10));
		strcpy(full_name, obj_name);
		strcpy(full_name+strlen(obj_name), "_output.m");
		//printf("%s\n", full_name);
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, full_name, &printer);
		PetscViewerPushFormat(printer, PETSC_VIEWER_ASCII_MATLAB);
	}
	public:
	MatlabPrinter(EMat &A) {
		init((PetscObject)A.self);
		MatView(A.self, printer);
	}
	MatlabPrinter(EVec &b) {
		init((PetscObject)b.self);
		VecView(b.self, printer);
	}
};
 
void row2position(int i, int &x, int &y) {//向量映射到矩阵
	x = i%x_size;
	y = i/x_size;
}

int position2row(int x, int y) {//矩阵映射到向量
	return y*x_size+x;
}

int uMap(int x, int y) {//矩阵映射到向量
	return y*x_size+x;
}

bool isBoundary(int x, int y) {//矩阵边界判别
	if (x == 0 || x == x_size-1)
		return true;
	if (y == 0 || y == y_size-1)
		return true;
	return false;
}

int main(int argc, char* argv[]) {
	//初始化MPI环境
	PetscInitialize(&argc, &argv, (char*)0, help);
	int size, rank;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	PetscErrorCode ierr;

	//系数矩阵和右端项
	EMat A(AIJ, x_size*y_size, x_size*y_size);
	EVec b(x_size*y_size);

	float *x = new float[x_size];
	for (int i = 0; i < x_size; i++)
		x[i] = h*i;
	float *y = new float[y_size];
	for (int i = 0; i < y_size; i++)
		y[i] = h*i;

	//设置系数矩阵
	if (rank == 0)
	for (int i = 0; i < x_size; i++)
		for (int j = 0; j < y_size; j++) {
			if (isBoundary(i, j)) {
				if (i == 0) {
					
	FEI newFEI(1, position2row(i,j));
	newFEI.set(uMap(i,j), 1);
	newFEI.rightItem = (sin(pi*y[j]))/((float)9+pow(pi, (float)2));

					A.setValues(newFEI);
					b.setValue(newFEI);
				}
				else if (i == x_size-1) {
					
	FEI newFEI(1, position2row(i,j));
	newFEI.set(uMap(i,j), 1);
	newFEI.rightItem = -(sin(pi*y[j]))/((float)9+pow(pi, (float)2));

					A.setValues(newFEI);
					b.setValue(newFEI);
				}
				else if (j == 0) {
					
	FEI newFEI(1, position2row(i,j));
	newFEI.set(uMap(i,j), 1);
	newFEI.rightItem = (float)0;

					A.setValues(newFEI);
					b.setValue(newFEI);
				}
				else if (j == y_size-1) {
					
	FEI newFEI(1, position2row(i,j));
	newFEI.set(uMap(i,j), 1);
	newFEI.rightItem = (float)0;

					A.setValues(newFEI);
					b.setValue(newFEI);
				}
			}
			else {
				
	FEI newFEI(5, position2row(i,j));
	newFEI.set(uMap(i,j), 1);
	newFEI.set(uMap(i-1,j), -((float)1)/((float)4));
	newFEI.set(uMap(i+1,j), -((float)1)/((float)4));
	newFEI.set(uMap(i,j-1), -((float)1)/((float)4));
	newFEI.set(uMap(i,j+1), -((float)1)/((float)4));
	newFEI.rightItem = (pow(h, (float)2))/((float)4)*f(x[i],y[j]);

				A.setValues(newFEI);
				b.setValue(newFEI);
			}
		}
	A.startAsm();
	b.startAsm();
	b.endAsm();
	A.endAsm();

	//解方程
	EVec xx(x_size*y_size);
	(new EKsp())->calculate(A, xx, b);

	//向量映射到矩阵
	xx.startGetValues();
	EMat ans(DENSE, y_size, x_size);
	for (int i = xx.start; i < xx.end; i++) {
		int xi, yi;
		row2position(i, xi, yi);
		ans.setValue(yi, xi, xx[i]);
	}
	xx.endGetValues();
	ans.startAsm();
	ans.endAsm();

	//打印矩阵
	ans.setName("ans");
	new MatlabPrinter(ans);

	return 0;
}
