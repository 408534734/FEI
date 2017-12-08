#include <petscsys.h>
#define int PetscInt
#define float PetscScalar

class FEI {
public:
	int n, row, *columns;
	float *values, rightItem;
	int top;
	FEI(int n, int row);
	void set(int column, float value);
	~FEI();
};