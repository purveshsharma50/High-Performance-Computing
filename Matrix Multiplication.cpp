// Name- Purvesh Sharma

#include<iostream>
#include <cstdlib>
#include <ctime>
using namespace std;
#define MTX_SIZE 10000

int main(){
    double sumtime=0.0;
    clock_t begin = clock();                  // To find Total time of program
    float ** sum= new float*[MTX_SIZE];
    int count= 0;       int pulse=0;                      // to find count of program for finding average
    srand(time(0));
    float** a=new float*[MTX_SIZE];             // Dynamic memory
    float**b=new float*[MTX_SIZE];
    float **Average=new float*[MTX_SIZE];
    float** results=new float*[MTX_SIZE];
for(int i=0;i<MTX_SIZE;i++){
sum[i]=new float[MTX_SIZE];
a[i]=new float[MTX_SIZE];
b[i]=new float[MTX_SIZE];
Average[i]= new float[MTX_SIZE];
results[i]= new float[MTX_SIZE];
}
    float min=-199.999;
    float max=199.999;
    double averagetime=0.0;
    for(int i=0;i<1000;i++){                      // to execute program 1000 times
        count++;
    for(int i=0;i<MTX_SIZE;i++){                                //Matrix size
        for(int j=0;j<MTX_SIZE;j++){
            float r=(float)rand()/(float)RAND_MAX;               // to keep the program with in limit
            float v;
            v=min+r*(max-min);
            a[i][j]=v;
    }

}

    for(int i=0;i<MTX_SIZE;i++){
        for(int j=0;j<MTX_SIZE;j++){                           //Matrix Size
            float r2=(float)rand()/(float)RAND_MAX;           // To keep matrix with in limit
            float v2;
            v2=min+r2*(max-min);
            b[i][j]=v2;
    }

}

 // To store Result of product of Matrix 1 and Matrix 2

 clock_t begin = clock();                                 //To get the matrix multiplication time
                                         // To keep count for finding average time
for (int i = 0; i < MTX_SIZE; i++){                                // Loop to do product(rows x Column)
    for (int j = 0; j < MTX_SIZE; j++){
        for(int k=0;k<MTX_SIZE;k++){
        results[i][j] += a[i][k] * b[k][j];                    //Matrix Multiplication


        }
         sum[i][j]+=results[i][j];

}
}
pulse++;
 clock_t end = clock();
  double elapsed_secmat = double(end - begin) / CLOCKS_PER_SEC;      //To calculate multiplication time
 sumtime+=elapsed_secmat;
  averagetime=sumtime/pulse;


for (int i = 0; i < MTX_SIZE; i++){                         // Loop to print output Matrix
    for (int j = 0; j < MTX_SIZE; j++){
        Average[i][j]=sum[i][j]/count;
    }
}


}

        cout << "Average:"<<endl;
for (int i = 0; i < MTX_SIZE; i++){                         // Loop to print output Matrix
    for (int j = 0; j < MTX_SIZE; j++){


        cout<<Average[i][j]<<" ";
    }
    cout << endl;
}
 clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout<< " Calculation Time: "<<elapsed_secs<<endl;                    // Output Total Time
  cout<<"Average Time: "<<averagetime<<endl;                           //Output Average Multiplication Time

delete[]sum;
delete[]a;
delete[]b;
delete[]Average;
delete[]results;
return 0;
}
