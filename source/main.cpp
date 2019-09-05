#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <chrono>
#include "GS.hpp"
using std::cout;
using std::endl;
using std::vector;


void printvector2(vector<double> &m, int row, int col)
{
    for(int i=0;i<row ;++i)
    {
        cout<<"| ";
        for (int j=0;j < col;++j)
        {
            cout<<m[row*i+j]<<" | ";
        }
        cout<<endl;
    }
}


int main(int argc, char* argv[]){
    int print = 0;
    if(argc < 2){
        return -1;
    }
    if(argc == 3){
        print = atoi(argv[2]);
    }
    double boundary_temp[4] = {100,0,0,0};
    double interior_temp = 50;
    int size = atoi(argv[1]);
    double tolerance = 1/(powf(size,2));
    vector<double> u(size*size);
    vector<double> w(size*size);
    int id;
    int ierr;
    int procs;
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    ierr = MPI_Init ( &argc, &argv );
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &procs );
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );
    if(id == 0){
        init_vector(u, boundary_temp, interior_temp, size);
/*        int k = steady_state(u,tolerance,size,size);
        cout << "seq:"<<k<<endl;
        printvector(u,size,size);
        init_vector(u, boundary_temp, interior_temp, size);*/
        start = std::chrono::high_resolution_clock::now();
    }
    int iteration = mpi_steady_state(u,tolerance,size,size,id,procs);
    if(id==0){
        end = std::chrono::high_resolution_clock::now();
        auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        
        if(print){
            printvector(u,size,size);
        }else{
            cout << "iteration:"<<iteration<<" millisec:"<< millis<<endl;
        }
    }
    MPI_Finalize ();
    return 0;
}
