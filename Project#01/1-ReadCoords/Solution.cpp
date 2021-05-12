#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>

using namespace std;

int main()
{
	Molecule mol("D:/input/acetaldehyde.dat", 0);

	cout << "Number of atoms: " << mol.natom << endl;
	cout << "Input Cartesian coordinates:\n";
	mol.print_geom();

	return 0;
}