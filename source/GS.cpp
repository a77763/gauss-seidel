#include "GS.hpp"
void init_vector(vector<double>& res, double* boundary_temp, double interior_temp, int size)
{
    vector<double> mat(size*size);
    for(int i = 0; i < size; ++i){
        mat[i] = boundary_temp[LOWER];
        mat[i*size] = i < size-1 ? boundary_temp[LEFT]: boundary_temp[ABOVE];
        mat[i*size+size-1] = boundary_temp[RIGHT]; 
        mat[(size-1)*size+i] = boundary_temp[ABOVE]; 
    }
    for(int i = 1; i < size-1; ++i){
        for(int j = 1; j < size-1; ++j){
            mat[i*size+j] = interior_temp;
        }
    }
    res.swap(mat);
}

void printvector(vector<double> m, int row, int col)
{
    for(int i=0; i<row ;++i)
    {
        for (int j=0; j < col;++j)
        {
            cout<<m[col*i+j]<<";";
        }
        cout<<endl;
    }
}



void red_state(vector<double>& u, int row, int col, double& max)
{
    bool flag = true;
    int r = 0;
    int b = 0;
    double tmp = 0.0;
    for (int i=1; i < row-1; ++i){
        for (int j=1; j < col-1; ++j){
            if(flag){
                tmp = u[i*col+j];
                u[i*col+j]=  
                            (
                                u[(i-1)*col+j]  +
                                u[i*col+j-1]    +
                                u[i*col+j+1]    +
                                u[(i+1)*col+j]
                            ) *0.25;
                tmp = fabs(u[i*col+j] - tmp);;
                max = max < tmp ? tmp : max; 
                r++; 
            }
            else{
                b++;
            }
            flag = !flag;
        }
        if(r==b){
            flag = !flag;
        }
        r = b = 0;
    }
}

void black_state(vector<double>& u, int row, int col,  double& max)
{
    bool flag = false;
    int r = 0;
    int b = 0;
    double tmp = 0.0;
    for (int i=1; i < row-1; ++i){
        for (int j=1; j < col-1; ++j){
            if(flag){
                tmp = u[i*col+j];
                u[i*col+j]=  
                            (
                                u[(i-1)*col+j]  +
                                u[i*col+j-1]    +
                                u[i*col+j+1]    +
                                u[(i+1)*col+j]
                            ) *0.25;
                tmp = fabs(u[i*col+j] - tmp);
                max = max < tmp ? tmp : max; 
                b++; 
            }
            else{
                r++;
            }
            flag = !flag;
        }
        if(r==b){
            flag = !flag;
        }
        r = b = 0;
    }
}

int steady_state(vector<double>& u,double tolerance, int row, int col)
{
    int iteration = 0;
    double diff = 0;
    vector<double> temp(row*col);
    do{
        diff = 0;
        red_state(u,row,col,diff);
        black_state(u,row,col,diff);
        ++iteration;
    }while(diff>tolerance);

    return iteration;

}




void mpi_ss_comm(
    int id,
    vector<double>& u_local, 
    int row, 
    int col)
{
    int procs;
    vector<double> u_cpy_snd(col);   
    vector<double> u_cpy_fst(col);
    
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    if (id != 0){
        MPI_Request w_req_fst;
        std::copy(
            u_local.data() + col,
            u_local.data() + 2*col,
            u_cpy_fst.data());
        MPI_Isend(u_cpy_fst.data(), col, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD, &w_req_fst);
    }
    if (id + 1 != procs) {
        MPI_Request w_req_snd;
        std::copy(
            u_local.data() + u_local.size() - col*2,
            u_local.data() + u_local.size() - col,
            u_cpy_snd.data());
        MPI_Isend(u_cpy_snd.data(), col, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD, &w_req_snd);
    }

    if(id != 0)
    {
        MPI_Status recv_stat;
        MPI_Recv( 
            u_local.data(), col, MPI_DOUBLE,
            id-1, 0, MPI_COMM_WORLD, &recv_stat);
    }
    if(id+1 != procs)
    {
        MPI_Status recv_stat;
        MPI_Recv( 
            u_local.data() + u_local.size() - col, col, MPI_DOUBLE,
            id+1, 0, MPI_COMM_WORLD, &recv_stat);
    }
}


void mpi_ss_error(double local_diff, double & diff)
{
    MPI_Allreduce(&local_diff, &diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void mpi_ss_iteration(bool pattern, vector<double>& u_local, int row, int col,  double& max)
{
    if (pattern) {
        red_state(u_local,row,col, max);    
    }
    else{
        black_state(u_local,row,col,max);
    }
}

void mpi_ss_init(
    vector<int> dis,
    vector<int> count,
    vector<double>& u, 
    vector<double>& u_local, 
    int size){
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatterv(
        u.data(), count.data(), dis.data(), MPI_DOUBLE, 
        u_local.data(), size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }





void mpi_gather_result(int id, int local_row, int col, vector<int> dis, vector<int> count, vector<double> u_local, vector<double>& u){
    if(id == 0)
    {
        std::transform(dis.begin(),   dis.end(),   dis.begin(),   [col](int& d) {return d + col;});
        std::transform(count.begin(), count.end(), count.begin(), [col](int& c) {return c - 2*col;}); 
    }
    MPI_Gatherv(
        u_local.data() + col, (local_row-2)*col, MPI_DOUBLE, 
        u.data(), count.data(), dis.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
}



int mpi_steady_state(vector<double>& u, double tolerance, int row, int col, int id, int procs)
{
    int work = round((row-2)/procs);
    vector<int> dis;
    vector<int> count;
    int iteration=0;
    double diff = 0;
    double local_diff = 0;
    bool terminate = false;
    bool first = procs*work % 2 ? true : false;
    bool second = !first;
    int local_row;
    for (int p = 0; p < procs; ++p)
    {
        dis.push_back(p * work * col);
        count.push_back((work + 2)*col);
        local_row = (work + 2);
    }
    
    if ((count[procs-1] + dis[procs-1]) != row){
        count[procs-1] = (row*col - dis[procs-1]);
        if(id == procs-1)
            local_row = (row*col - dis[procs-1])/col; 
    }
    vector<double> u_local(count[id]);
    vector<double> w_local(count[id]);




    mpi_ss_init(dis, count, u, u_local, count[id]);
    do{

        local_diff = 0.0;
        mpi_ss_iteration(first, u_local, local_row, col, local_diff);
        mpi_ss_comm(id, u_local,local_row, col);
        mpi_ss_iteration(second, u_local, local_row, col, local_diff); 
        mpi_ss_error(local_diff, diff);
        if(diff > tolerance){
            mpi_ss_comm(id, u_local, local_row, col);
            ++iteration;
        }
        else
        {
            break;
        }
    }while(true);

    mpi_gather_result(id, local_row, col, dis, count, u_local, u);
    
    return iteration;
}

