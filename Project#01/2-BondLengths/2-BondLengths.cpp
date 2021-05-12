#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;

int main()
{
    Molecule mol("D:/input/acetaldehyde.dat", 0);

    cout << "Number of atoms: " << mol.natom << endl;
    cout << "Input Cartesian coordinates:\n";
    mol.print_geom();

    cout << "Interatomic distances (bohr):\n";
    for (int i = 0; i < mol.natom; i++)
        for (int j = 0; j < i; j++)
            printf("%d %d %8.5f\n", i, j, mol.bond(i, j));

    return 0;
}