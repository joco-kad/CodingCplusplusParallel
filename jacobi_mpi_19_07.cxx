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
#include <stdlib.h>
	
#include <sstream>


using namespace std;

void JacobiStep(double * x,const double* M,const int mDim, int start, int steps);
void read(const char *name,double *&M, int &mDim);
int pos(const int i, const int j, const int mDim);
// bool cmpr(const double* a,const double* b,const int mDim,const double precision);
double cmpr(const double* a,const double* b,const int mDim);

void writeToFile(const string filename,const double *x,const int mDim);
double maxi(const double* a,const double* b,const int mDim);

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
  double precision = 1E-8;//= 0.0000001;
  int mDim = 0;
  double *M;
  
//   double *M2, *M3 , *M4;
  
  int finished = 0;
  double *x;
  int * recvcount = new int[tasks];
  int * offset = new int[tasks];
  double* c;
  int testcounter=0;
  double norm=0.0;
  double norm_old=1.0;
  double norm_pre=0.0;
  clock_t t;
 

  
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
    
    // id 0 needs comparison vector to check wether we are done.
    
        t = clock();
  }
 
  // send dim, x vector and matrix to all
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Bcast(&mDim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier( MPI_COMM_WORLD );
  if(id>0){
   x = new double[mDim];
   M = new double[mDim*(mDim+1)];
  }  
 
//   M2 = &M[mDim*(mDim+1)/4];
//   M3 = &M[(mDim*(mDim+1)/4)*2];
//   M4 = &M[(mDim*(mDim+1)/4)*3];
//   
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Bcast(x,mDim,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Bcast(recvcount,tasks,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Bcast(offset,tasks,MPI_INT,0,MPI_COMM_WORLD);
   c = new double[mDim];
      for(int i = 0; i < mDim; i++){
	c[i] = 0;
      } 
  //divide Matrix for sending
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Bcast(M,mDim*(mDim+1),MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Barrier( MPI_COMM_WORLD );
//   MPI_Bcast(M2,mDim*(mDim+1)/4,MPI_DOUBLE,0,MPI_COMM_WORLD);
//   MPI_Barrier( MPI_COMM_WORLD );
//   MPI_Bcast(M3,mDim*(mDim+1)/4,MPI_DOUBLE,0,MPI_COMM_WORLD);
//   MPI_Barrier( MPI_COMM_WORLD );
//   MPI_Bcast(M4,mDim*(mDim+1)/4,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Barrier( MPI_COMM_WORLD );


//--------------------------------------------
//----------------MAIN------------------------
//--------------------------------------------

// everybody works until we are done
  while(cmpr(x,c,mDim)>precision){
    for(int i = 0; i < mDim; i++){
	c[i] = x[i];
      }
    // everybody makes  jacobi steps for rows depending on his id
    JacobiStep(x,M,mDim,1+id*(mDim/tasks),mDim/tasks);// rest?
    // if there is a rest the task with the highest id has to do it
    if(id ==(tasks-1) && mDim%tasks!=0)
      JacobiStep(x,M,mDim,1+mDim-(mDim%tasks),mDim%tasks);

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allgatherv(&x[offset[id]], recvcount[id], MPI_DOUBLE, 
                   x,recvcount,offset, 
                   MPI_DOUBLE, MPI_COMM_WORLD);

    // id 0 compares x to x from last step if to a certain precision equal stop iterating
    double arschloch = cmpr(x,c,mDim);
    cout << " Rank : " << id << ", max Fehler : " << arschloch << endl;
//     if((arschloch <precision)){ //(id == 0) && 
//       Check for convergence
//       norm_old= norm;
//       norm=maxi(x,c,mDim);
//       norm_pre = norm/norm_old;
//       if(norm_old != 0.0)
//       cout<<"Vorfaktor (const.))="<<"\t"<<norm_pre<<endl;
//       update comparison vector
//     }
    
       if(id == 0 && (testcounter%1000)==0){
	 // print to file
	
	string Result;
	stringstream convert;
	convert.str(" ");
	convert << testcounter; 
	Result = convert.str();
	writeToFile(Result,x,mDim);
	 // we are done
	finished = 1;
	
      }

      // let everyone know we are done
     MPI_Barrier( MPI_COMM_WORLD );
     MPI_Bcast(&finished,1,MPI_INT,0,MPI_COMM_WORLD);
     testcounter++;
   }
   if(id == 0 ){
	 // print to file
	writeToFile(argv[2],x,mDim);
	 // we are done
	finished = 1;
	
      }
   delete[] c;
   cout << id << "\t" << testcounter << endl;
  if(id==0){ 
    t = clock() - t;
    fstream out;
    out.open ("benchgauss", std::fstream::out | std::fstream::app);
    out << tasks<< "\t" << ((double)t)/CLOCKS_PER_SEC << endl;
    out.close();
  }
  // free memory and finalize
  delete[] M;
  delete[] x;
  delete[] recvcount;
  delete[] offset;
  MPI_Finalize();
  return 0;
}


/**
* Calculates one Step of a Jacobi Iteration 	
* with x being the approx. solution, M the Matrix of the System and 
* b the right hand side.
* x = 1/aii (bi - Sum(aij * xj)) 
*/
void JacobiStep(double * x,const double* M,const int mDim, int start, int steps){

  // double h to store intermediate results
  // We need x for calculation and cannot overwrite it until
  // we have completed the whole step
  double h;
  
  //go through columns; initialize auxiliary variable
  for(int i =start; i< (start+steps); i++){
    h=0.0;
    // sum over  row
    for(int j = 1; j <= mDim; j++){
     if(i != j) h += M[pos(i,j,mDim)]*x[j-1];
    }
    h = (M[pos(i,mDim+1,mDim)] -h)/M[pos(i,i,mDim)];
    x[i-1]=h;
  }
}

/**
 * reads Data from file.Format file as follows
 * Add b vector to right hand side of Matrix A 
 * divided by tabs.
 */
void read(const char *name,double* &M, int &mDim){
  ifstream in(name);
  int i = 0;
  char c;
  // count all tabs in first line add 1 to get dimension
  while((c= in.get()) && (c != '\n')){
    if(c=='\t') 
      i++;
  }
  in.close();
  mDim = i;
  // read file and save into matrix
  in.open(name);
  M = new double[i*(i+1)];
  double h;
  for(int m = 1; m<=i; m++){
    for(int n = 1; n<= i+1; n++){
      in>>h;
      M[pos(m,n,i)] = h;
    }
  }
}

/**
 * Function calculates position of Matrix element in array.
 */
int pos(const int i, const int j, const int mDim){
	return ((i-1)*(mDim+1) + j-1);
}


/**
 * Functions to vectors returns true if identical to precision.
 */
// bool cmpr(const double* a,const double* b,const int mDim,const double precision){
//   double max = 0;
//   bool bla = false;
//           for(int i = 0; i<mDim;i++){
// 	    double max2 = abs(a[i]) - abs(b[i]);
// 	    if ( max2 > max)
// 	      max = max2;
// 	    if( abs(a[i] - b[i]) <  precision){
// 	      bla = true;
// 	    }
// 	   }
// 	   cout << "maximale abweichung : " << max << endl;
//         return bla;
// 	
// }
double cmpr(const double* a,const double* b,const int mDim){
  double max = 0;
          for(int i = 0; i<mDim;i++){
	    double max2 = abs(a[i]) - abs(b[i]);
	    if ( max2 > max)
	      max = max2;
	    }
	   cout << "maximale abweichung : " << max << endl;
        return max;
	
}


double maxi(const double* a,const double* b,const int mDim){
	double absolut=0.0;
	  for(int i = 0; i<mDim;i++){
	    if(absolut<a