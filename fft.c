

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>
#include <string.h>
#include <time.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
//21 + 1
#define numOfChars 22
#define Pi 3.14159265358979323846264338

// Global variables
int N;                            // Number of data points of the padded data (allways pow of 2)
double complex omega;             // omega in the fft-algorithm
int startIndex;                   // the global index where this processors elements starts
int endIndex;                     // the global index where this processors elements ends
int nofElementsOnTheProcessors;   // Number of elements on this processor (N/P)
int p;                            // This processor's index
int P;                            // Nof processors
int L;
int R;


int lenOfFile(char fileName[]) {

    /*
    Input: File name, to a file with doubles saved as strings seperated with \n.
    Outup: Number of elements in the file.
    */

    FILE* pf; //pointer file
    pf = fopen(fileName, "r");
    fseek(pf, 0, SEEK_END); // seek to end of file
    int size = ftell(pf); // get current file pointer
    fseek(pf, 0, SEEK_SET); // seek back to beginning of file
    int len = size/numOfChars;
    return len;
}

bool isPowOfTwo(int n){

    /*
    Input:  Integer n, we wish to check if power of 2.
    Output: True if power of two, false if 0 or not power of 2.
    */

    if (n == 0)
        return false;

    return (ceil(log2(n)) == floor(log2(n)));
}

int* int_to_bin(int k , int lenOfBink){

  /*
  Input: Integer k, which we want to know the binary represenation of
  Output: The binary represenation of k as an array.
  */

  int* bink = calloc(lenOfBink,sizeof(int));

  if(k == 0){
    return bink;
  }

  //printf("k = %i\n", k);
  int ind = lenOfBink - 1;

  int rem;

  while (ind>-1){
    rem = k%2;
    bink[ind] = rem;
    k = k / 2;
    ind--;
  }
  return bink;

}

int nextPowofTwo(int n){

  /*
  Input:  Integer n, we wish to know the next power of 2 if n not allready a
  power of 2.
  Output: True if power of two, false if 0 or not power of 2
  */

  return pow(2, ceil(log2(n)));
}


int* int_to_arr(long num, int lengthOfNum){

  /*
  Input: Integer num; Ie 112
  Output: Vector of num; ie [1,1,2]
  */

  int * returnVec = malloc(lengthOfNum*sizeof(int));
  //printf("lenOfNum: %i\n", lengthOfNum);
  while(lengthOfNum--){
    returnVec[lengthOfNum] = num%10;
    num = floor(num/10);
  }

  return returnVec;
}


int* binRepFlip(int num, int b_m, int m, int r){

  /*
  Input: Integer num we want the binary representation of with length r and the m-th bit set to b_m
  Output: Array with the binary representation of num with the m-th bit set to b_m
  Ie: num = 7, b_m=0, m=2 r = 4 ---> [0,1,0,1]
  */

  int *binaryRep = calloc(r,sizeof(int));
  int* binOfNum = malloc(ceil(log2(num)*sizeof(int)));

  int lenOfBin;
  if(num == 0||num == 1){
    lenOfBin = 1;
  } else if(isPowOfTwo(num)){
    lenOfBin = ceil(log2(num)) +1;
  } else {
    lenOfBin = ceil(log2(num));
  }

  binOfNum = int_to_bin(num, lenOfBin);


  int rr = r-lenOfBin;
  binaryRep[m] = b_m;
  for(int ind = 0; ind < lenOfBin; ind++){

    if(rr == m){;} // pass
    else{
      binaryRep[rr] = binOfNum[ind];
    }
    rr++;
  }

  return binaryRep;
}

int* binRep(int num, int r){

  /*
  Input: Integer num we want the binary representation of with length r
  Output: Array with the binary representation of num
  Ie: num = 7, r = 4 ---> [0,1,1,1]
  */

  int *binaryRep = calloc(r,sizeof(int));

  int lenOfBin;
  if(num == 0){
    return binaryRep;
  } else if(num == 1){
    lenOfBin = 1;
  } else if(isPowOfTwo(num)){
    lenOfBin = ceil(log2(num)) +1;
  } else {
    lenOfBin = ceil(log2(num));
  }

  int* binOfNum = malloc(ceil(log2(num)*sizeof(int)));
  binOfNum = int_to_bin(num, lenOfBin);

  int rr = r-lenOfBin;
  for(int ind = 0; ind < lenOfBin; ind++){
    binaryRep[rr] = binOfNum[ind];
    rr++;
  }

  return binaryRep;
}

int binRepToInt(int* binRep,int lengthOfBinRep){

  /*
  Input: pointer to array-representaion of a int:
  Output: The represenation as an int of the base 10 represenation of the input
  */

  int returnVal = 0;
  int pow = 1;
  for(int ii = 0; ii<lengthOfBinRep; ii++){
    returnVal = returnVal + binRep[lengthOfBinRep-ii-1]*pow;
    pow = pow*2;
  }

  return returnVal;
}

void reverseArr(int* arr,int n){

  /*
  Reverses the array
  */

  for(int left = 0, right = n-1; left<right; left++, right--){
    int tt = arr[left];
    arr[left] = arr[right];
    arr[right] = tt;
  }
}

int createOmegaPow(int* bini, int m, int r){

  /*
  Input: binary representation in an array of indeex i
  Output: the power omega should be raised to
  */

  int* binRepOfPow = calloc(r,sizeof(int));

  for(int ii = 0; ii<m+1; ii++){
    binRepOfPow[ii] = bini[m-ii];
  }

  int decRepOfPow = binRepToInt(binRepOfPow,r);
  return decRepOfPow;
}

int* bitFlip(int* flip, int m, int r){
  /*
  Input: Flip the m-th bit in flip
  Output: the fliped version of Flip
  */

  int* returnVec = flip;
  if(flip[m] == 0){returnVec[m] = 1;}
  else{returnVec[m] = 0;}
  return returnVec;

}

double complex* fftp_seq(double complex * X, int n){

  /*
  Input: vector x, we wish to compute the fourier transform of.
         int n, number of elements on this processor.
  Output: the fourier transform of x.
  */

  double complex* returnVec = malloc(nofElementsOnTheProcessors * sizeof(double complex));

  int r = log2(N);  // num of bits in index (N is global)
  int d = log2(P);  // P = 2^d
  double complex * R =  malloc(nofElementsOnTheProcessors * sizeof(double complex));
  double complex * S =  malloc(nofElementsOnTheProcessors * sizeof(double complex));

  for(int ii = 0; ii < nofElementsOnTheProcessors; ii++) {R[ii] = X[ii];}

  for(int m = 0; m<r; m++){                                     // Outer loop

    for(int ii = 0; ii<nofElementsOnTheProcessors; ii++) {S[ii] = R[ii];}
    int localInd = 0;

    double complex* recivedS = malloc(nofElementsOnTheProcessors * sizeof(double complex));

    int commWith;
    for(int i = startIndex; i<endIndex; i++){                   // inner loop

      // Work for the fft step
      int* j_b = calloc(r,sizeof(int)); int* k_b = calloc(r,sizeof(int));  // j and k in binary
      int j; int k;        // j and k in base 10
      j_b = binRepFlip(i, 0, m, r);
      k_b = binRepFlip(i, 1, m, r);
      j = binRepToInt(j_b, r);
      k = binRepToInt(k_b, r);
      int* bini;            // i in binary
      bini = binRep(i,r);
      int omegaPow = createOmegaPow(bini, m, r);
      //printf("i = %i, j = %i, k = %i, pow = %i\n", i, j, k, omegaPow);
      if(localInd == 0){    // Communications only needed in the first inner iteration and before the d-th outer iteration
        if(m < d){

          int* commWithRowBin = malloc(sizeof(int)*r);
          int commWithRow;

          commWithRowBin = bitFlip(bini, m, r);
          commWithRow = binRepToInt(commWithRowBin, r);
          commWith =  floor(commWithRow/L);

          MPI_Status status;
          MPI_Sendrecv(S, nofElementsOnTheProcessors,
                      MPI_C_DOUBLE_COMPLEX, commWith,
                      123, recivedS,
                      nofElementsOnTheProcessors,
                      MPI_C_DOUBLE_COMPLEX,commWith,
                      123, MPI_COMM_WORLD,
                      &status);
        }
      }

      // fft computations
      //printf("i = %i \t j = %i\t k = %i\t pow = %i\n",i,j,k,omegaPow);
      if(m >= d){                                // No comm needed
        //printf("i = %i \t j = %i\t k = %i\t pow = %i\n",i,j,k,omegaPow);
        int local_j = j - p*L;
        int local_k = k - p*L;

        R[localInd] = S[local_j] + S[local_k]*cpow(omega, omegaPow);
        //printf("\t R[%i] = %.15f + %.15f*I\n----\n",i,creal(R[localInd]), cimag(R[localInd]));
      }
      else{                                      // Communications needed

        if(commWith > p){ // We got S[k]:s
          int local_j = j - p*L;
          int local_k = k - p*L;
          //printf("IF: pp = %i, loc_j = %i, loc_k = %i\n",p, local_j, local_k);
          R[localInd] = S[local_j] + recivedS[localInd]*cpow(omega, omegaPow);
        }
        else{           // we got S[j]:s
          int local_k = k - p*L;
          //printf("ELSE: pp = %i, loc_j = %i, loc_k = %i\n",p, local_j, local_k);
          R[localInd] = recivedS[localInd] + S[local_k]*cpow(omega, omegaPow);
        }
      }
      localInd++;
    }                                             // End inner loop
  }                                               // End inner loop
  for(int ii = 0; ii<nofElementsOnTheProcessors; ii++){returnVec[ii] = R[ii];}
  free(R);
  free(S);
  return returnVec;
}

void main(int argc, char** argv) {

    //char fileName[] = "mp3_dataCH1.txt";
    char fileName[] = "wav_dataCH1.txt";

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes, P is global
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    //printf("P=%i\n", P);

    // Get the rank of the process, p is global
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

    // Check if P is pow of 2.
    if(p == 0){
      if(!isPowOfTwo(P)){
        printf("\nNumber of processors MUST be a power of 2\n");
        return;
      }
    }

    // number of data points in the file
    int N_data = lenOfFile(fileName);

    //  make sure N is pow of 2, else make it the next pow of 2, N is global
    if(!isPowOfTwo(N_data)) {
      N = nextPowofTwo(N_data);   // Data did not cointain pow of 2 data points
    }
    else {
      N = N_data;                 // N_data was allready a pow of 2
    }

    // Omega in the algorithm
    omega = cexp((2*I*Pi)/N); // omega is global

    L = N / P;
    R = N % P;  // R = 0
    nofElementsOnTheProcessors = (N + P - p - 1) / P;       // nof elements on THIS procsessor
    startIndex = p * L + MIN(p, R);                         // Is global
    endIndex = startIndex + nofElementsOnTheProcessors;     // Is global

    // File reading and save data in malloc arrays
    int startByte = startIndex * numOfChars; // Which byte in the file the procsessor should start reading from
    double complex* data = malloc(nofElementsOnTheProcessors * sizeof(double complex));

    FILE* dataFile;
    dataFile = fopen(fileName, "r");
    if(!dataFile){
      printf("proc:%i - Error opening file\n",p);
      return;
    };

    fseek(dataFile, startByte, SEEK_SET);
    double num;
    for (int i = 0; i<nofElementsOnTheProcessors; i++){
        if (i+startIndex < N_data){
          fscanf(dataFile, "%lf",&num);
          data[i] = num + 0.0*I;
        }
        else {
          data[i] = 0.0 + 0.0*I;
        }
    }

    fclose(dataFile);
    /*
    for(int ii = 0; ii < nofElementsOnTheProcessors; ii++){
      printf("proc = %i\t %e\n",p, creal(data[ii]));
    }
    */
    printf("proc: %i - Reading OK\n",p);

    // Start fft and timing
    double complex* fftOfData = malloc(sizeof(double)*nofElementsOnTheProcessors);

    clock_t start, end;
    double cpu_time_used;
    int TIMES = 1;
    start = clock(); // START TIMER
    for(int ii = 0; ii<TIMES; ii++){
      fftOfData = fftp_seq(data, nofElementsOnTheProcessors); // FFT
    }
    end = clock(); // END TIMER

    free(data);

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cpu_time_used = cpu_time_used/TIMES;
    FILE* output_time;
    output_time = fopen("clocking.dat", "a");
    fprintf(output_time, "%i \t %f \t %i \t %i \n", p, cpu_time_used, P, N);
    fclose(output_time);

    // Checking output
    /*
    for(int ii = 0; ii < nofElementsOnTheProcessors; ii++){
      printf("proc = %i\t %e + %e*1i\n",p, creal(fftOfData[ii]),cimag(fftOfData[ii]));
    }
    */

    // File writing
    int toLeft = p - 1;
    int toRight = p + 1;

    // We have to fix the boundarys
    if (toLeft == -1){
      toLeft = MPI_PROC_NULL;
    }
    if (toRight == P){
      toRight = MPI_PROC_NULL;
    }

    int NextGO = 1;
    int timeToGo;
    FILE *output;
    if(p == 0){
      output = fopen("outputFFT.dat", "w");
      for(int i = 0; i<nofElementsOnTheProcessors; i++){
        fprintf(output, "%e + %e*1i\n",creal(fftOfData[i]),cimag(fftOfData[i]));
      }
      fclose(output);
      MPI_Send(&NextGO, 1, MPI_INT, toRight, 100, MPI_COMM_WORLD);
    }
    else{
      MPI_Recv(&timeToGo, 1, MPI_INT, toLeft, 100, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      output = fopen("outputFFT.dat", "a");
      for(int i = 0; i<nofElementsOnTheProcessors; i++){
        fprintf(output, "%e + %e*1i\n",creal(fftOfData[i]),cimag(fftOfData[i]));
        //printf(output, "%.13d + %.13d*1i\n",creal(data[i]),cimag(data[i]));
      }
      fclose(output);
      MPI_Send(&NextGO, 1, MPI_INT, toRight, 100, MPI_COMM_WORLD);
    }
    free(fftOfData);
    MPI_Finalize();

}


/*
// Not uesed - reqursive fft
double complex* fftp_rec(double complex * x,int N){


  //Input: vector x, we wish to compute the fourier transform of.
  //Output: the fourier transform of x.


  double complex omega_N = cexp((-2*I*Pi)/N);

  // N == 1 means end of recursion
  double complex* returnVec = malloc(N * sizeof(double complex));

  if(N==1){
    returnVec[0] = x[0] + 0*I; // make it complex
  }
  else{
    double complex w[N/2];
    double complex evenIndex_x[N/2];
    double complex oddIndex_x[N/2];

    for(int ind = 0; ind < N/2; ind++){
      complex double k = ind + 0*I;
      w[ind] = cpow(omega_N, k);
      oddIndex_x[ind] = x[2*ind + 1];
      evenIndex_x[ind] = x[2*ind];
    }

    double complex* u = fftp_rec(evenIndex_x, N/2);
    double complex* v_ = fftp_rec(oddIndex_x, N/2); // will be dot multipled with w

    // Dot multiply v_ with w and create return vec
    double complex v;
    for(int ind = 0; ind < N/2; ind++){
      v = w[ind]*v_[ind];
      returnVec[ind] = u[ind] + v;
      returnVec[N/2 + ind] = u[ind] - v;
    }
  }

  return returnVec;

}

*/
