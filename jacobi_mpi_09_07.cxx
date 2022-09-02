/**
 * Program for solving a system of equations by Jacobi iteration.
 * Serial edition.
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

using namespace std;

void JacobiStep(double * &x,const double* M,const int mDim, int start, int steps);
void read(const char *name,double *&M, int &mDim);
int pos(const int i, const int j, const int mDim);
bool cmpr(const double* &a,const double* &b,const int mDim,const double precision);
void writeToFile(const string filename,const double *&x,const int mDim);

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
  //tags
  int mat = 1;
  int dim = 2;
  int xvec = 3;
  
  // memory stuff
  double precision = 0.0000001;
  int mDim = 0;
  double *M;
  bool finished = false;
  double *x;
  int * recvcount = new int[tasks];
  int * offset = new int[tasks];
  double* c;

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
      offset[i] = i*(mDim/tasks);
    }
    recvcount[tasks-1] = mDim/tasks + (mDim%tasks);
    offset[tasks-1] = (tasks-1)*(mDim/tasks);
    
    // id 0 needs comparison vector to check wether we are done.
    c = new double[mDim];
      for(int i = 0; i < mDim; i++){
	c[i] = x[i];
      }
//     cout << mDim << endl;
//     //M = new double[mDim*(mDim+1)];
//     for(int m = 1; m<=mDim; m++){
//       for(int n = 1; n<= mDim+1; n++){
// 	cout << M[pos(m,n,mDim)]<<"\t";
//     }
//     cout << endl;
//   }

   //send Matrix to all tasks
//    for(int dest = 1;dest<tasks;dest++){
//     // MPI_Send(&M,M.getnDim()*(M.getnDim()+1),MPI_DOUBLE,dest,mat,MPI_COMM_WORLD);
//    }
  }
  
  // send dim, x vector and matrix to all
  MPI_Bcast(&mDim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&M,mDim*(mDim+1),MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&x,mDim,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&recvcount,tasks,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&offset,tasks,MPI_INT,0,MPI_COMM_WORLD);
  
//--------------------------------------------
//----------------MAIN------------------------
//--------------------------------------------

// everybody works until we are done
  while(!finished){
    
    // everybody makes  jacobi steps for rows depending on his id
    JacobiStep(x,M,mDim,1+id*mDim/tasks,mDim/tasks);// rest?
    // if there is a rest the task with the highest id has to do it
    if(id ==(tasks-1) && mDim%tasks!=0)
      JacobiStep(x,M,1+mDim-(mDim%tasks),mDim%tasks);
    
    MPI_Allgatherv(x[offset[i]], recvcounts[i], MPI_DOUBLE, 
                   x,recvcounts,offset, 
                   MPI_DOUBLE, MPI_COMM_WORLD);
    
    // id 0 compares x to x from last step if to a certain precision equal stop iterating
    if((id == 0) && (cmpr(x,c,precision) == false)){
      //update comparison vector
      for(int i = 0; i < mDim; i++){
	c[i] = x[i];
      }
     else if(id == 0){
	 // print to file
	writeToFile(argv[2],x,mDim);
	 // we are done
	finished = true;
      }
    }
    // let everyone know we are done
    MPI_Bcast(&finished,1,MPI_BOOL,0,MPI_COMM_WORL