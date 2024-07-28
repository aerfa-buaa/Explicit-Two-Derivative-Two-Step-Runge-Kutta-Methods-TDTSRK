/* Explicit Two-Derivative Two-Step Runge-Kutta Methods
 * Copyright (C) 2024
 * Author QinXueyu
 */

#include "xBurgers.hpp"

using namespace std;

int main()
{

	cout << "----------------------------------------" << endl;
	cout << "------------Start  Order----------------" << endl;

	CSolver *realsolver = NULL;

	realsolver = new CADSolver();

	realsolver->Driver();

	if (realsolver != NULL)
	{
		delete realsolver;
	}

	cout << "----------------------------------------" << endl;
	cout << "----------------end-------------------" << endl;
 
	return 0;
}