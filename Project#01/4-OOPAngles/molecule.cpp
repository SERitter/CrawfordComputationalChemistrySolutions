#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>

void Molecule::print_geom()
{
    for (int i = 0; i < natom; i++)
        printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

void Molecule::translate(double x, double y, double z)
{
    for (int i = 0; i < natom; i++) {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

double Molecule::bond(int atom1, int atom2)
{
    return sqrt((geom[atom1][0] - geom[atom2][0]) * (geom[atom1][0] - geom[atom2][0])
        + (geom[atom1][1] - geom[atom2][1]) * (geom[atom1][1] - geom[atom2][1])
        + (geom[atom1][2] - geom[atom2][2]) * (geom[atom1][2] - geom[atom2][2]));
}

double Molecule::angle(int atom1, int atom2, int atom3)
{
    return acos(unit(0, atom2, atom1) * unit(0, atom2, atom3) + unit(1, atom2, atom1) * unit(1, atom2, atom3) + unit(2, atom2, atom1) * unit(2, atom2, atom3));
}

double Molecule::oop(int atom1, int atom2, int atom3, int atom4)
{
    double ebcd_x = (unit(1, atom3, atom2) * unit(2, atom3, atom4) - unit(2, atom3, atom2) * unit(1, atom3, atom4));
    double ebcd_y = (unit(2, atom3, atom2) * unit(0, atom3, atom4) - unit(0, atom3, atom2) * unit(2, atom3, atom4));
    double ebcd_z = (unit(0, atom3, atom2) * unit(1, atom3, atom4) - unit(1, atom3, atom2) * unit(0, atom3, atom4));

    double exx = ebcd_x * unit(0, atom3, atom1);
    double eyy = ebcd_y * unit(1, atom3, atom1);
    double ezz = ebcd_z * unit(2, atom3, atom1);

    double theta = (exx + eyy + ezz) / sin(angle(atom2, atom3, atom4));

    if (theta < -1.0) theta = asin(-1.0);
    else if (theta > 1.0) theta = asin(1.0);
    else theta = asin(theta);

    return theta;
}

double Molecule::unit(int cart, int atom1, int atom2)
{
    return -(geom[atom1][cart] - geom[atom2][cart]) / bond(atom1, atom2);
}

Molecule::Molecule(const char* filename, int q)
{
    charge = q;

    // open filename
    std::ifstream is(filename);
    assert(is.good());

    // read the number of atoms from filename
    is >> natom;

    // allocate space
    zvals = new int[natom];
    geom = new double* [natom];
    for (int i = 0; i < natom; i++)
        geom[i] = new double[3];

    for (unsigned int i = 0; i < natom; i++)
        is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

    is.close();
}

Molecule::~Molecule()
{
    delete[] zvals;
    for (int i = 0; i < natom; i++)
    {
        delete[] geom[i];
    }
    delete[] geom;
}