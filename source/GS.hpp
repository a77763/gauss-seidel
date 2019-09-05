#ifndef _GS_
#define _GS_
#include <iostream>
#include <vector>
#include <mpi.h>
#include <cmath>
#include <algorithm>
using std::cout;
using std::endl;
using std::vector;
using std::begin;
using std::end;
#define LOWER 0
#define LEFT 1
#define RIGHT 2
#define ABOVE 3
void printvector(vector<double>m, int row, int col);
void init_vector(vector<double>& mat, double* boundary_temp, double interior_temp, int size);
int steady_state(vector<double>& u,double tolerance, int row, int col);
int mpi_steady_state(vector<double>& u,  double tolerance, int row, int col, int id, int procs);
#endif