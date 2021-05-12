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

    cout << "\nBond angles:\n";
    for (int i = 0; i < mol.natom; i++) {
        for (int j = 0; j < i; j++) {
            for (int k = 0; k < j; k++) {
                if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0)
                    printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i, j, k) * (180.0 / acos(-1.0)));
            }
        }
    }

    return 0;
}