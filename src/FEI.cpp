#include <FEI.h>
#include <petscsys.h>

#define int PetscInt
#define float PetscScalar

FEI::FEI(int n, int row) : n(n), row(row) {
	columns = new int[n];
	values = new float[n];
	top = 0;
}

void FEI::set(int column, float value) {
	this->columns[top] = column;
	this->values[top] = value;
	top++;
}
FEI::~FEI() {
	delete[] columns;
	delete[] values;
}