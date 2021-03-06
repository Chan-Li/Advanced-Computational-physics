#include <stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#define N 20000
#define N 20000
int main(int argc, char* argv[])
{
    int myid,numprocs;
    int namelen,numsend,numrec;
    int i,j,l,k;
    int outcome;
    double start, stop;
    double test;
    double temp=0;
    double mode=0;
    int row_num;
    double *vec = NULL;
    double *matrix_=NULL;
    double *result=NULL;
    double *buffer = NULL;
    long double *scale = NULL;
    
    
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    row_num = N/(numprocs-1);
    vec=malloc(N*sizeof(double));
    matrix_=malloc(N*N*sizeof(double));
    result=malloc(N*sizeof(double));
    buffer = malloc(N*sizeof(double));
    test=0;
    long double *result_scaling=NULL;
    result_scaling=malloc(N*sizeof(long double));

    double *ans = NULL;
    ans=malloc(row_num*sizeof(double));
    scale=malloc(row_num*sizeof(long double));
    for(i=0;i<row_num;i++)
        ans[i]=0;
    for(i=0;i<N;i++)
        result[i]=0;

    if (numprocs<=1)
    {
        MPI_Finalize();
        return 0;
    }

    if (myid == 0)
    {
        start = MPI_Wtime();
        srand((unsigned)time(NULL)+123 * (0 + 2324817813));
        for(j=0;j<N;j++)
        {
            vec[j]=rand()/(RAND_MAX+1.0);
        }
        srand((unsigned)time(NULL)+123 * (1 + 2324817813));
        for(j=0;j<(N);j++)
        {
            for(i=0;i<N;i++)
            {
                matrix_[j*N+i]=rand()/(RAND_MAX+1.0);
            }
        }
    }
    MPI_Bcast(vec,N,MPI_INT,0,MPI_COMM_WORLD);
    if (myid == 0)
       {
           numsend=0;
           while (numsend<N)
           {
               for (l=1; l<numprocs; l++)
               {
                   for(j=0;j<(N);j++)
                   {
                       buffer[j]=matrix_[(numsend)*N+j];//rand()%10;
                   }
                    MPI_Send(buffer,(N),MPI_DOUBLE_PRECISION,l,l,MPI_COMM_WORLD);
                    numsend+=1;
               }
           }
           if (numsend>=N)
           {
               for (l=1; l<numprocs; l++)
               {
                   MPI_Send(NULL,0,MPI_INT,l,0,MPI_COMM_WORLD);
               }
               
           }
           outcome=0;
           while(outcome<N)
           {
               for(k=1;k<numprocs;k++)
               {
                   MPI_Recv(&test,1, MPI_DOUBLE_PRECISION, k, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                   result[outcome] = test;
                   outcome+=1;
               }
           }
           mode=0;
           for(k=0;k<N;k++)
           {
               mode+=(result[k]*result[k]);
           }
           if (mode>1e10)
           {
               mode = sqrt(mode/1e10)*(1e5);
           }
           else
           {
               mode = sqrt(mode);
           }
           //mode = (sqrt(mode));
           //printf("Mode%f\n",mode);
           
       }
    //
    
     else
       {
           numrec=0;
           while (numrec<N)
           {
               MPI_Recv(buffer,(N), MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
               if (status.MPI_TAG==0)
               {
                   break;
               }
               else
               {
                   numrec+=1;
                   test = 0;
                   for(i=0;i<N;i++)
                   {
                       test = test+vec[i]*buffer[i];
                   }
                   ans[numrec-1]=test;
                   MPI_Send(&test,1,MPI_DOUBLE_PRECISION,0,99,MPI_COMM_WORLD);
               }
               
           }
           
                                        
       }
    MPI_Bcast(&mode,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);
    if (myid>0)
    {
        for (i=0; i<row_num; i++)
        {
            scale[i] = ans[i]/mode;
            temp=scale[i];
            MPI_Send(&temp,1,MPI_DOUBLE_PRECISION,0,99,MPI_COMM_WORLD);
        }
        //printf("The mode %f\n",mode);

    }
    if (myid==0)
    {
        outcome=0;
        while(outcome<N)
        {
            for(k=1;k<numprocs;k++)
            {
                MPI_Recv(&temp,1, MPI_DOUBLE_PRECISION, k, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                result_scaling[outcome] = temp;
                outcome+=1;
            }
        }
        for(k=0;k<10;k++)
        {
            printf("The scaled vector %Lf\n",result_scaling[k]);
        }
        stop = MPI_Wtime();

        printf("%lf\n",stop-start);
        free(result);
        free(vec);
        free(matrix_);
        free(buffer);
        free(result_scaling);
        free(scale);
        free(ans);

    }

    MPI_Finalize();
    return 0;


}
