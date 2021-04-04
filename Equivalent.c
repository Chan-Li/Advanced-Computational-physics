#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#define N 20000
#define M 20000
int main()
{
    double *vec=NULL;
    double *mat=NULL;
    double *result_extra=NULL;
    double *mat2=NULL;
    double mode;
    double *result2=NULL;
    double *result_scaling=NULL;
    int myid;
    double start, stop;
    int numprocs;
    int row_num,row_num2;
    int i,j,k;
    double *result=NULL;
    double *all_together=NULL;
    vec=malloc(M*sizeof(double));
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    row_num=(int)N/numprocs;
    row_num2 = N%numprocs;
    result_extra=malloc(row_num2*sizeof(double));
    mat2=malloc(M*row_num2*sizeof(double));
    srand((unsigned)time(NULL)+123 * (5+12324817813));
    for(j=0;j<M;j++)
        vec[j]=rand()/(RAND_MAX+1.0);
    if(myid==0)
    {
        start = MPI_Wtime();
    }

    
    mat=malloc(row_num*M*sizeof(double));
    result=malloc(row_num*sizeof(double));
    result2=malloc(row_num*sizeof(double));
    

    
    for(i=0;i<row_num;i++)
        result[i]=0;
    for(i=0;i<row_num2;i++)
        result_extra[i]=0;
    
    
    srand((unsigned)time(NULL)+123 * (myid + 2324817813));
    for(i=0;i<row_num;i++)
        for(j=0;j<M;j++)
            mat[i*M+j]=rand()/(RAND_MAX+1.0);
    srand((unsigned)time(NULL)+123 * (myid + 2324817812));
    for(i=0;i<row_num2;i++)
        for(j=0;j<M;j++)
            mat2[i*M+j]=rand()/(RAND_MAX+1.0);
    MPI_Barrier(MPI_COMM_WORLD);

    for(i=0;i<row_num;i++)
    {
        for(j=0;j<M;j++)
        {
            
            result[i]+=mat[i*M+j]*vec[j];
        }
    }




    MPI_Barrier(MPI_COMM_WORLD);
        
    if(myid==0)
    {
        all_together=malloc(N*sizeof(double));
        for(k=0;k<N;k++)
        {
            all_together[k]=0;
        }

        MPI_Gather
            (
            result,
            row_num,
            MPI_DOUBLE,
            all_together,
            row_num,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
            );
    }
    else
    {

        MPI_Gather
            (
            result,
            row_num,
            MPI_DOUBLE,
            all_together,
            row_num,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
            );
    }

   
    MPI_Barrier(MPI_COMM_WORLD);

   
    if(myid==0)
    {
        if (row_num2>0)
        {
            for(i=0;i<(row_num2);i++)
            {
                for(j=0;j<(M);j++)
                {
                    all_together[row_num*numprocs+i]+= vec[j]*mat2[i*M+j];
                }
                
            }
        }
//        for(k=0;k<10;k++)
//                {
//                    printf("This is the result%f\n",all_together[k]);
//                }


        
        for(k=0;k<N;k++)
        {
            //printf("Result %f\n",all_together[k]);
            mode+=all_together[k]*all_together[k];
        }
        mode = sqrt(mode);
        //printf("The mode is %f\n",mode);
     }
    MPI_Bcast(&mode,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);

//
    for(i=0;i<row_num;i++)
    {
        result2[i] = result[i]/mode;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid==0)
    {
        
        result_scaling=malloc(N*sizeof(double));

        
        MPI_Gather
            (
            result2,
            row_num,
            MPI_DOUBLE,
            result_scaling,
            row_num,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
            );
    }
    else
    {
        MPI_Gather
            (
            result2,
            row_num,
            MPI_DOUBLE,
            result_scaling,
            row_num,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
            );
    }


   
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0)
    {
        if (row_num2>0)
        {
            for(i=0;i<row_num2;i++)
            {
                result_scaling[row_num*(numprocs)+i]=all_together[row_num*(numprocs)+i]/mode;
            }
        }
        for(k=0;k<10;k++)
        {
            printf("This is the scaling result%f\n",result_scaling[k]);
        }

        stop = MPI_Wtime();
        printf("%lf\n",stop-start);
    }



    

    MPI_Finalize();
    /**********************/

    free(vec);
    free(mat);
    free(result);
    free(all_together);
    free(result2);
    free(result_scaling);

    return 0;
}
