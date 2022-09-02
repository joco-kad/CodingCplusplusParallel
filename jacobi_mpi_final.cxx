/**
 * Program for solving a system of equations by Jacobi iteration.
 * MPI edition.
 * authors: Jörg and Jan.
 * 
 * call as follows: ./NAME INPUTFILENAME OUTPUTFILENAME
 * whereas FILENAME is the name of a file storing a matrix M of size
 * m x (m+1). 
 * It solves for x According to: Ax = b
 * with M being the augmented coefficient matrix (A|b).
 */

//include awsome stuff
#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <mpi.h>
#include <time.h> 
#include <sstream>


using namespace std;

void JacobiStep(double * x,const double* M,const int mDim, int start, int steps);
void read(const char *name,double *&M, int &mDim);
int pos(const int i, const int j, const int mDim);
double maxnorm(const double* a,const double* b,const int mDim);
void writeToFile(const string filename,const double *x,const int mDim);

//begin main here
int main(int argc, char *argv[]){
//--------------------------------------------
//----------------MPI-------------------------
//--------------------------------------------
  //Begin with MPI statements
  MPI_Status status;
  int tasks, id;

  // MPI initialize 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &tasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  
  // memory stuff
  double precision = 1E-9;
  int mDim = 0;
  double *M;
  
  
  int finished = 0;
  double *x;
  int * recvcount = new int[tasks];
  int * offset = new int[tasks];
  double* c;
  int stepcounter=0;
  clock_t t;
 
  // Divide Bcast if mDim*mDim > n*n to avoid make sure data fits in network buffer
  int n = 5000;

  
  // Master task 0 reads Matrix from file then distributes dimension and matrix to all other tasks
  if (id==0){
    // read Matrix from file
    read(argv[1],M,mDim);
    
    // create guess for solution
    x = new double[mDim];
    for(int i = 0; i< mDim;i++)
      x[i] = 1.0;
    
    // prepare recvcount array which stores the number of steps every task does later
    // offset array which stores the offsets in receiving array
    for(int i = 0; i < tasks-1;i++){
      recvcount[i] = mDim/tasks;
      offset[i] = i*recvcount[i];//(mDim/tasks);
    }
    recvcount[tasks-1] = mDim/tasks + (mDim%tasks);
    offset[tasks-1] = (tasks-1)*(mDim/tasks);
    
    //clock time
    t = clock();
  }
 
  // send dim, x vector and matrix to all
  MPI_Bcast(&mDim,1,MPI_INT,0,MPI_COMM_WORLD);
  // all others need to create their arrays after they know the matrix dimension
  if(id>0){
   x = new double[mDim];
   M = new double[mDim*(mDim+1)];
  }  
  

  MPI_Bcast(x,mDim,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(recvcount,tasks,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(offset,tasks,MPI_INT,0,MPI_COMM_WORLD);
  // everybody creates an comparison vector (everbody should finish at the same time)
   c = new double[mDim];
      for(int i = 0; i < mDim; i++){
	c[i] = 0;
      }
      
  // if matrix is big divide broadcast in small packets    
  if((mDim*(mDim+1)) > (n*n)){
    int div = (mDim*(mDim+1))/(n*n);
    for(int i = 0; i<div; i++)
      MPI_Bcast(&M[i*((mDim*(mDim+1))/div)],(mDim*(mDim+1))/div,MPI_DOUBLE,0,MPI_COMM_WORLD);
    // remainder
    MPI_Bcast(&M[(div-1)*((mDim*(mDim+1))/div)],(mDim*(mDim+1))%(n*n),MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  else  
    MPI_Bcast(M,mDim*(mDim+1),MPI_DOUBLE,0,MPI_COMM_WORLD);
  


//--------------------------------------------
//----------------MAIN------------------------
//--------------------------------------------

// everybody works until we are done
  while(maxnorm(x,c,mDim)>precision){
    //update comparison vector
    for(int i = 0; i < mDim; i++){
	c[i] = x[i];
      }
    // everybody makes  jacobi steps for rows depending on his id
    JacobiStep(x,M,mDim,1+id*(mDim/tasks),mDim/tasks);// rest?
    // if there is a rest the task with the highest id has to do it
    if(id ==(tasks-1) && mDim%tasks!=0)
      JacobiStep(x,M,mDim,1+mDim-(mDim%tasks),mDim%tasks);

    // distribute solution to all tasks
    MPI_Allgatherv(&x[offset[id]], recvcount[id], MPI_DOUBLE, 
                   x,recvcount,offset, 
                   MPI_DOUBLE, MPI_