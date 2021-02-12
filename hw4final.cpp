#include <mpi.h>                     //MPI Library
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <omp.h>                     //Open MP library
#define MTX_SIZE 10000

using namespace std;
//initialization of MPI Coding
int main(int argc, char* argv[])
{
	int rank, size, rc, tag, offset, rows;
	double start1= MPI_Wtime();
    
    srand(time(0));

// contigeous memory allocation

    double**a=new double*[MTX_SIZE];
	 double**b=new double*[MTX_SIZE];
	 double**results=new double*[MTX_SIZE];
	 double**sum=new double*[MTX_SIZE];
	 double**Average=new double*[MTX_SIZE];
	 int matrix_size =MTX_SIZE*MTX_SIZE;
     a[0]=new double[matrix_size];
	 b[0]=new double[matrix_size];
	 results[0]=new double[matrix_size];
	 sum[0]=new double[matrix_size];
	 Average[0]=new double[matrix_size];
	 
	for(int i = 1; i < MTX_SIZE; i++)
		{                    
		a[i] =&a[0][i*MTX_SIZE];
		b[i]=&b[0][i*MTX_SIZE];
		results[i]=&results[0][i*MTX_SIZE];
		sum[i]=&sum[0][i*MTX_SIZE];
		Average[i]=&Average[0][i*MTX_SIZE];
	}
	for (int i = 0; i < MTX_SIZE; i++) {
	    for (int j = 0; j < MTX_SIZE; j++) {
	        sum[i][j] = 0.0;
	    }
	}
// Limit for random generation

    double min = -199.999;
    double max = 199.999;
    double duration_sec_sum=0;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2) {
        cout << "c u later" << endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }
    //  THE MASTER CODE

	for(int i=1;i<=1000;i++){
	    if (rank == 0) {
	      
		        for (int i = 0; i < MTX_SIZE; i++) {
		            for (int j = 0; j < MTX_SIZE; j++) {
		                float v1 = (max - min) * ((((float)rand()) / (float)RAND_MAX)) + min;
		                float v2 = (max - min) * ((((float)rand()) / (float)RAND_MAX)) + min;
		                a[i][j] = v1;
		                b[i][j] = v2;
		            }
		        }
	
		        // Send matrix data to the workers
		        rows = MTX_SIZE / (size - 1);
		        offset = 0;
		        tag = 1; // From MASTER
		
 // MPI_Send (buffer, count, type, dest, tag, comm)
		
		        for (int dest = 1; dest <= size - 1; dest++) {
		            MPI_Send(&offset, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
		            MPI_Send(&rows, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
		            MPI_Send(&a[offset][0], rows * MTX_SIZE, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
		            MPI_Send(&b[0][0], MTX_SIZE * MTX_SIZE, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
		            offset += rows;
		        }
		
		    tag = 2; // From WORKER
		    int source;
		    // MPI_Recv(buffer, count, type, source, tag, comm, status)
		    for (int i = 1; i <= size - 1; i++) {
		        source = i;
		        MPI_Recv(&offset, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		        MPI_Recv(&rows, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		        MPI_Recv(&results[offset][0], rows * MTX_SIZE, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
		    }

	  }
   
	  
	    // END of MASTER CODE
	 
	    // START of WORKER CODE
	   
	
	    /*      MPI_Send (buffer, count, type, dest, tag, comm)
	                MPI_Recv (buffer, count, type, source, tag, comm, status)       */
	    if (rank > 0) {
	        tag = 1; //From MASTER
	        MPI_Recv(&offset, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
	        MPI_Recv(&rows, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
	        MPI_Recv(&a[0][0], rows * MTX_SIZE, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	        MPI_Recv(&b[0][0], MTX_SIZE * MTX_SIZE, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	
	        double start2= MPI_Wtime();
//Matrix Multiplication with open MP
	        #pragma omp parallel for
	        for (int k = 0; k < MTX_SIZE; k++) {
	            for (int i = 0; i < rows; i++) {
	                results[i][k] = 0.0;
	                for (int j = 0; j < MTX_SIZE; j++) {
	                    results[i][k] += a[i][j] * b[j][k];
	                }
	            }
	        }
	        double end2= MPI_Wtime();
		    double duration_sec = ( end2 - start2);
		    duration_sec_sum+= duration_sec;
	        
	        tag = 2; //From WORKER
	        MPI_Send(&offset, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
	        MPI_Send(&rows, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
	        MPI_Send(&results[0][0], rows * MTX_SIZE, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	    }
	    if (rank==0){
	    	for(int k = 0; k < MTX_SIZE; k++){
					for(int i = 0; i < rows; i++){
					    sum[i][k]+=results[i][k];
					}
				}	
	    	
		}
	    
	}
	
	// Printing of matrix
    if (rank==0){
    	cout << "final Matrix:"<< endl;
	   
	    for (int i = 0; i < MTX_SIZE; i++) {
	        for (int j = 0; j < MTX_SIZE; j++) {
	            cout << sum[i][j]/1000 << " ";
	        }
	        cout << endl;
	    }
	    cout<<endl;
	    double average_time=duration_sec_sum/1000;
	    cout<< "Average Time: "<< average_time<<" "<<endl;
	    double end1= MPI_Wtime();
	    double duration_sec_final = (end1 - start1);
        cout<<" Total Time:"<< duration_sec_final<<" ";
	}
	
    
    // END of WORKER Code
    
    
    MPI_Finalize(); 
	delete[]a;
	delete[]b;
	delete[]results;
	delete[]sum;
	delete[]Average;  
return 0;
}