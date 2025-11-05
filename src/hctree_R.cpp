#include <Rcpp.h>
#include <iostream>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ordering.h"

#ifdef __cplusplus
extern "C"
{
#endif

class htreeNode
{
public:
    int left;
    int right;
    double distance;
};

int* ivector1(int length)
{
  int *out;
  out=new int[length+1];
  //out--;
  return out;
}

void free_ivector1(int* vec)
{
  //vec++;
  delete[] vec;
}

double* dvector1(int length)
{
  double *out;
  out=new double[length+1];
  //out--;
  return out;
}

void free_dvector1(double* vec)
{
  //vec++;
  delete[] vec;
}

int** imatrix1(int nrow,int ncol)
{
  int **out,i=0;
  out=new int*[nrow+1];
  for(i=0;i<=nrow;i++){
    out[i]=new int[ncol+1];
    //out[i]--;
  }
  //out--;
  return out;
}

void free_imatrix1(int** mat,int nrow)
{
  int i=0;
  //mat++;
  for(i=0;i<=nrow;i++){
    //mat[i]++;
    delete[] mat[i];
  }
  delete[] mat;
}

double** dmatrix1(int nrow,int ncol)
{
  double **out;
  int i=0;
  out=new double*[nrow+1];
  for(i=0;i<=nrow;i++){
    out[i]=new double[ncol+1];
    //out[i]--;
  }
  //out--;
  return out;
}

void free_dmatrix1(double** mat,int nrow)
{
  int i=0;
  //mat++;
  for(i=0;i<=nrow;i++){
    //mat[i]++;
    delete[] mat[i];
  }
  delete[] mat;
}

    double** fun1dto2dArray(double *input_ptr, int n, int p)
    {
        /*
        double** output_ptr = new double* [p];
        for(int i=0; i<p; i++)
        {
            output_ptr[i] = new double[n];   
        }
        */
        double** output_ptr = dmatrix1(n, p);

        for(int i=0; i<n*p; i++)
        {
            int x = i/n;
            int y = i%n;
            output_ptr[x+1][y+1]=input_ptr[i];
        }

        return output_ptr;
    }


static int nodecompare(const void* a, const void* b)
/* Helper function for qsort. */
{
    const htreeNode* node1 = (const htreeNode*)a;
    const htreeNode* node2 = (const htreeNode*)b;
    const double term1 = node1->distance;
    const double term2 = node2->distance;
    if (term1 < term2) return -1;
    if (term1 > term2) return +1;
    return 0;
}

static double find_closest_pair(int n, double** distmatrix, int* ip, int* jp)
/*
This function searches the distance matrix to find the pair with the shortest
distance between them. The indices of the pair are returned in ip and jp; the
distance itself is returned by the function.

n          (input) int
The number of elements in the distance matrix.

distmatrix (input) double**
A ragged array containing the distance matrix. The number of columns in each
row is one less than the row index.

ip         (output) int*
A pointer to the integer that is to receive the first index of the pair with
the shortest distance.

jp         (output) int*
A pointer to the integer that is to receive the second index of the pair with
the shortest distance.
*/
{
    int i=0, j=0;
    double temp=0;
    double distance=distmatrix[1][0];
    *ip = 1;
    *jp = 0;
    for(i=1; i< n; i++)
    {
        for(j=0; j<i; j++)
        {
            temp=distmatrix[i][j];
            if(temp<distance)
            {
                distance=temp;
                *ip=i;
                *jp=j;
            }
        }
    }
    return distance;
}

double dmin(double value1,double value2)
{
    if(value1 >value2)return value2;
    else return value1;
}

double dmax(double value1,double value2)
{
    if(value1 >value2)return value1;
    else return value2;
}



void htree_single_d(double **_data,int n,int **merge,double *hgt)
{
/*
  int i,j,n2,iter,*minwith,*active,ptmin,pm1,pm2;
  double **data,min,max,*minval;
  n2=2*n-1;
  //allocate
  data=dmatrix1(n2,n2);
  active=ivector1(n2);
  minwith=ivector1(n2);
  minval=dvector1(n2);
  //
  max=_data[1][1];
  //_data[1][7]=_data[1][7];
  //_data[7][1]=_data[7][1];
  for(i=1;i<n;i++){
    for(j=(i+1);j<=n;j++){
      if(_data[i][j]>max){
        max=_data[i][j];
      }
    }
  }
  max=max+1.0;
  for(i=1;i<=n;i++){
    merge[i][1]=merge[i][2]=0;
    hgt[i]=_data[i][i];
    active[i]=1;
    minval[i]=max;
    for(j=1;j<=n;j++){
      data[i][j]=_data[i][j];
      if((i != j) && (minval[i] > data[i][j])){
        minval[i]=data[i][j];
        minwith[i]=j;
      }
    }
  }
  for(j=(n+1);j<=n2;j++){
    active[j]=0;
  }
  iter=n+1;
  while(iter<=n2){
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(minval[i] < min){
          ptmin=i;
          min=minval[i];
        }
      }
    }
    hgt[iter]=min;
    if(ptmin < minwith[ptmin]){
      pm1=merge[iter][1]=ptmin;
      pm2=merge[iter][2]=minwith[ptmin];
    }else{
      pm1=merge[iter][1]=minwith[ptmin];
      pm2=merge[iter][2]=ptmin;
    }
    active[pm1]=0;
    active[pm2]=0;
    //
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(data[i][pm1]<data[i][pm2]){
          data[i][iter]=data[iter][i]=data[i][pm1];
        }else{
          data[i][iter]=data[iter][i]=data[i][pm2];
        }
        //
        if((minwith[i] == pm1) || (minwith[i] == pm2)){
          minwith[i]=iter;
          minval[i]=data[i][iter];
        }
        //
        if(min > data[iter][i]){
          min=data[iter][i];
          ptmin=i;
        }
      }
    }
    minval[iter]=min;
    minwith[iter]=ptmin;
    active[iter]=1;
    iter++;
  }
  free_dmatrix1(data,n2);
  free_dvector1(minval);
  free_ivector1(minwith);
  free_ivector1(active);
*/

    int i=0,j=0,k=0;
    int nelements = n;
    int nnodes = n-1;
    int* vector;
    double* temp;
    int* index;
    htreeNode* result;

    temp = new double [nnodes];
    index = new int [nelements];
    vector = new int [nnodes];
    result = new htreeNode [nnodes];

    for(i=1;i<=n;i++)
    {
        merge[i][1]=merge[i][2]=0;
        hgt[i]=_data[i][i];
    }

    for(i=0; i<nnodes; i++)
    {
        vector[i] = i;
        result[i].distance = DBL_MAX;
    }


    for(i=0; i<n; i++)
    {
        for(j=0; j<i; j++)
        {
            temp[j] = _data[i+1][j+1];
        }
        for(j=0; j<i; j++)
        {
            k=vector[j];
            if(result[j].distance >= temp[j])
            {
                if(result[j].distance < temp[k])
                {
                    temp[k] = result[j].distance;
                }
                result[j].distance = temp[j];
                vector[j] = i;
            }
            else if(temp[j] < temp[k])
            {
                temp[k] = temp[j];
            }
        }
        for(j=0; j<i; j++)
        {
            if(result[j].distance >= result[vector[j]].distance)
            {
                vector[j] = i;
            }
        }
    }

    delete[] temp;

    for(i=0; i<nnodes; i++)
    {
        result[i].left = i;
    }
    qsort(result, nnodes, sizeof(htreeNode), nodecompare);

    for(i=0; i<nelements; i++)
    {
        index[i] = i;
    }
    for(i=0; i<nnodes; i++)
    {

        j = result[i].left;
        k = vector[j];
        result[i].left = index[j];
        result[i].right = index[k];
        index[k] = -i-1;


    }

    for(i=1; i<=nnodes; i++)
    {
        merge[n+i][1]=result[i-1].left;
        if(merge[n+i][1]<0)
            merge[n+i][1]=n-merge[n+i][1];
        else
            merge[n+i][1]=merge[n+i][1]+1;
        merge[n+i][2]=result[i-1].right;
        if(merge[n+i][2]<0)
            merge[n+i][2]=n-merge[n+i][2];
        else
            merge[n+i][2]=merge[n+i][2]+1;
        hgt[n+i]=result[i-1].distance;
    }

    delete[] vector;
    delete[] index;
    delete[] result;

}//End of Function

//

void htree_single_s(double **_data,int n,int **merge,double *hgt)
{
  int i=0,j=0,n2=0,iter=0,*minwith,*active,ptmin=0,pm1=0,pm2=0;
  double **data,min=0,max=0,*minval;
  n2=2*n-1;
  //allocate
  data=dmatrix1(n2,n2);
  active=ivector1(n2);
  minwith=ivector1(n2);
  minval=dvector1(n2);
  //
  max=_data[1][1];
  for(i=1;i<n;i++){
    for(j=(i+1);j<=n;j++){
      if(_data[i][j]<max){
        max=_data[i][j];
      }
    }
  }
  max=max-1.0;
  for(i=1;i<=n;i++){
    merge[i][1]=merge[i][2]=0;
    hgt[i]=_data[i][i];
    active[i]=1;
    minval[i]=max;
    for(j=1;j<=n;j++){
      data[i][j]=_data[i][j];
      if((i != j) && (minval[i] < data[i][j])){
        minval[i]=data[i][j];
        minwith[i]=j;
      }
    }
  }
  for(j=(n+1);j<=n2;j++){
    active[j]=0;
  }
  iter=n+1;
  while(iter<=n2){
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(minval[i] > min){
          ptmin=i;
          min=minval[i];
        }
      }
    }
    hgt[iter]=min;
    if(ptmin < minwith[ptmin]){
      pm1=merge[iter][1]=ptmin;
      pm2=merge[iter][2]=minwith[ptmin];
    }else{
      pm1=merge[iter][1]=minwith[ptmin];
      pm2=merge[iter][2]=ptmin;
    }
    active[pm1]=0;
    active[pm2]=0;
    //
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(data[i][pm1]>data[i][pm2]){
          data[i][iter]=data[iter][i]=data[i][pm1];
        }else{
          data[i][iter]=data[iter][i]=data[i][pm2];
        }
        //
        if((minwith[i] == pm1) || (minwith[i] == pm2)){
          minwith[i]=iter;
          minval[i]=data[i][iter];
        }
        //
        if(min < data[iter][i]){
          min=data[iter][i];
          ptmin=i;
        }
      }
    }
    minval[iter]=min;
    minwith[iter]=ptmin;
    active[iter]=1;
    iter++;
  }
  free_dmatrix1(data,n2);
  free_dvector1(minval);
  free_ivector1(minwith);
  free_ivector1(active);
}

//
// Hierachical Complete Linkage Tree
//

void htree_complete_d(double **_data,int n,int **merge,double *hgt)
{
/*
  int i,j,n2,iter,*minwith,*active,ptmin,pm1,pm2;
  double **data,min,max,*minval;
  n2=2*n-1;
  //allocate
  data=dmatrix1(n2,n2);
  active=ivector1(n2);
  minwith=ivector1(n2);
  minval=dvector1(n2);
  //
  max=_data[1][1];
  for(i=1;i<n;i++){
    for(j=(i+1);j<=n;j++){
      if(_data[i][j]>max){
        max=_data[i][j];
      }
    }
  }
  max=max+1.0;
  for(i=1;i<=n;i++){
    merge[i][1]=merge[i][2]=0;
    hgt[i]=_data[i][i];
    active[i]=1;
    minval[i]=max;
    for(j=1;j<=n;j++){
      data[i][j]=_data[i][j];
      if((i != j) && (minval[i] > data[i][j])){
        minval[i]=data[i][j];
        minwith[i]=j;
      }
    }
  }
  for(j=(n+1);j<=n2;j++){
    active[j]=0;
  }
  iter=n+1;
  while(iter<=n2){
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(minval[i] < min){
          ptmin=i;
          min=minval[i];
        }
      }
    }
    hgt[iter]=min;
    if(ptmin < minwith[ptmin]){
      pm1=merge[iter][1]=ptmin;
      pm2=merge[iter][2]=minwith[ptmin];
    }else{
      pm1=merge[iter][1]=minwith[ptmin];
      pm2=merge[iter][2]=ptmin;
    }
    active[pm1]=0;
    active[pm2]=0;
    //
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(data[i][pm1]>data[i][pm2]){
          data[i][iter]=data[iter][i]=data[i][pm1];
        }else{
          data[i][iter]=data[iter][i]=data[i][pm2];
        }
        //
        if((minwith[i] == pm1) || (minwith[i] == pm2)){
          minval[i]=max;
          for(j=1;j<iter;j++){
            if(active[j]){
              if(j != i){
                if(minval[i] > data[i][j]){
                  minval[i]=data[i][j];
                  minwith[i]=j;
                }
              }
            }
          }
        }
        //
        if(min > data[iter][i]){
          min=data[iter][i];
          ptmin=i;
        }
      }
    }
    minval[iter]=min;
    minwith[iter]=ptmin;
    active[iter]=1;
    iter++;
  }
  free_dmatrix1(data,n2);
  free_dvector1(minval);
  free_ivector1(minwith);
  free_ivector1(active);
*/
    int i,j;
    int m;
    int nelements=n;
    int* clusterid;
    htreeNode* result;
    double **data = new double *[nelements];
    for(i=0;i<nelements;i++)
    {
        data[i] = new double [nelements];
    }

    clusterid = new int [nelements];
    result = new htreeNode [nelements-1];

    for(i=1;i<=nelements;i++)
    {
        merge[i][1]=merge[i][2]=0;
        hgt[i]=_data[i][i];
    }

    for(i=0;i<nelements;i++)
    {
        for(j=0;j<nelements;j++)
        {
            data[i][j]=_data[i+1][j+1];
        }
    }
    /* Setup a list specifying to which cluster a gene belongs */
    for(j=0; j<nelements; j++)
    {
        clusterid[j]=j;
    }

    for(m=nelements; m>1; m--)
    {
        int is=1;
        int js=0;
        result[nelements-m].distance=find_closest_pair(m, data, &is, &js);

        /* Fix the distances */
        for(j=0; j<js; j++)
        {
            data[js][j]=dmax(data[is][j],data[js][j]);
        }
        for(j=js+1; j<is; j++)
        {
            data[j][js]=dmax(data[is][j],data[j][js]);
        }
        for(j=is+1; j<m; j++)
        {
            data[j][js]=dmax(data[j][is],data[j][js]);
        }

        for(j=0; j<is; j++)
        {
            data[is][j]=data[m-1][j];
        }
        for(j=is+1; j<m-1; j++)
        {
            data[j][is]=data[m-1][j];
        }

        /* Update clusterids */
        result[nelements-m].left=clusterid[is];
        result[nelements-m].right=clusterid[js];
        clusterid[js]=m-nelements-1;
        clusterid[is]=clusterid[m-1];
    }

    for(i=1; i<=nelements-1; i++)
    {
        merge[n+i][1]=result[i-1].left;
        if(merge[n+i][1]<0)
            merge[n+i][1]=n-merge[n+i][1];
        else
            merge[n+i][1]=merge[n+i][1]+1;
        merge[n+i][2]=result[i-1].right;
        if(merge[n+i][2]<0)
            merge[n+i][2]=n-merge[n+i][2];
        else
            merge[n+i][2]=merge[n+i][2]+1;
        hgt[n+i]=result[i-1].distance;
    }

    for(i=0;i<nelements;i++)
    {
        delete [] data[i];
    }
    delete [] data;
    delete[] clusterid;
    delete[] result;

}//End of Function

//

void htree_complete_s(double **_data,int n,int **merge,double *hgt)
{
  int i,j,n2,iter,*minwith,*active,ptmin,pm1,pm2;
  double **data,min,max,*minval;
  n2=2*n-1;
  //allocate
  data=dmatrix1(n2,n2);
  active=ivector1(n2);
  minwith=ivector1(n2);
  minval=dvector1(n2);
  //
  max=_data[1][1];
  for(i=1;i<n;i++){
    for(j=(i+1);j<=n;j++){
      if(_data[i][j]<max){
        max=_data[i][j];
      }
    }
  }
  max=max-1.0;
  for(i=1;i<=n;i++){
    merge[i][1]=merge[i][2]=0;
    hgt[i]=_data[i][i];
    active[i]=1;
    minval[i]=max;
    for(j=1;j<=n;j++){
      data[i][j]=_data[i][j];
      if((i != j) && (minval[i] < data[i][j])){
        minval[i]=data[i][j];
        minwith[i]=j;
      }
    }
  }
  for(j=(n+1);j<=n2;j++){
    active[j]=0;
  }
  iter=n+1;
  while(iter<=n2){
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(minval[i] > min){
          ptmin=i;
          min=minval[i];
        }
      }
    }
    hgt[iter]=min;
    if(ptmin < minwith[ptmin]){
      pm1=merge[iter][1]=ptmin;
      pm2=merge[iter][2]=minwith[ptmin];
    }else{
      pm1=merge[iter][1]=minwith[ptmin];
      pm2=merge[iter][2]=ptmin;
    }
    active[pm1]=0;
    active[pm2]=0;
    //
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(data[i][pm1]<data[i][pm2]){
          data[i][iter]=data[iter][i]=data[i][pm1];
        }else{
          data[i][iter]=data[iter][i]=data[i][pm2];
        }
        //
        if((minwith[i] == pm1) || (minwith[i] == pm2)){
          minval[i]=max;
          for(j=1;j<iter;j++){
            if(active[j]){
              if(j != i){
                if(minval[i] < data[i][j]){
                  minval[i]=data[i][j];
                  minwith[i]=j;
                }
              }
            }
          }
        }
        //
        if(min < data[iter][i]){
          min=data[iter][i];
          ptmin=i;
        }
      }
    }
    minval[iter]=min;
    minwith[iter]=ptmin;
    active[iter]=1;
    iter++;
  }
  free_dmatrix1(data,n2);
  free_dvector1(minval);
  free_ivector1(minwith);
  free_ivector1(active);
}

//
// Hierachical WPGMA Linkage Tree
//

void htree_average_d(double **_data,int n,int **merge,double *hgt)
{
/*
  int i,j,n2,iter,*minwith,*active,ptmin,pm1,pm2;
  double **data,min,max,*minval;
  n2=2*n-1;
  //allocate
  data=dmatrix1(n2,n2);
  active=ivector1(n2);
  minwith=ivector1(n2);
  minval=dvector1(n2);
  //
  max=_data[1][1];
  for(i=1;i<n;i++){
    for(j=(i+1);j<=n;j++){
      if(_data[i][j]>max){
        max=_data[i][j];
      }
    }
  }
  max=max+1.0;
  for(i=1;i<=n;i++){
    merge[i][1]=merge[i][2]=0;
    hgt[i]=_data[i][i];
    active[i]=1;
    minval[i]=max;
    for(j=1;j<=n;j++){
      data[i][j]=_data[i][j];
      if((i != j) && (minval[i] > data[i][j])){
        minval[i]=data[i][j];
        minwith[i]=j;
      }
    }
  }
  for(j=(n+1);j<=n2;j++){
    active[j]=0;
  }
  iter=n+1;
  while(iter<=n2){
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(minval[i] < min){
          ptmin=i;
          min=minval[i];
        }
      }
    }
    hgt[iter]=min;
    if(ptmin < minwith[ptmin]){
      pm1=merge[iter][1]=ptmin;
      pm2=merge[iter][2]=minwith[ptmin];
    }else{
      pm1=merge[iter][1]=minwith[ptmin];
      pm2=merge[iter][2]=ptmin;
    }
    active[pm1]=0;
    active[pm2]=0;
    //
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        data[i][iter]=data[iter][i]=(data[i][pm1]+data[i][pm2])/2.0;
        //
        if((minwith[i] == pm1) || (minwith[i] == pm2)){
          minval[i]=max;
          for(j=1;j<iter;j++){
            if(active[j]){
              if(j != i){
                if(minval[i] > data[i][j]){
                  minval[i]=data[i][j];
                  minwith[i]=j;
                }
              }
            }
          }
        }
        //
        if(min > data[iter][i]){
          min=data[iter][i];
          ptmin=i;
        }
      }
    }
    minval[iter]=min;
    minwith[iter]=ptmin;
    active[iter]=1;
    iter++;
  }
  free_dmatrix1(data,n2);
  free_dvector1(minval);
  free_ivector1(minwith);
  free_ivector1(active);
*/
    int i,j;
    int m;
    int nelements=n;
    int* clusterid;
    int* number;
    htreeNode* result;
    double **data = new double *[nelements];
    for(i=0;i<nelements;i++)
    {
        data[i] = new double [nelements];
    }

    clusterid = new int [nelements];
    number = new int [nelements];
    result = new htreeNode [nelements-1];

    for(i=1;i<=nelements;i++)
    {
        merge[i][1]=merge[i][2]=0;
        hgt[i]=_data[i][i];
    }

    for(i=0;i<nelements;i++)
    {
        for(j=0;j<nelements;j++)
        {
            data[i][j]=_data[i+1][j+1];
        }
    }

  /* Setup a list specifying to which cluster a gene belongs, and keep track
   * of the number of elements in each cluster (needed to calculate the
   * average). */
    for(j=0; j<nelements; j++)
    {
        number[j]=1;
        clusterid[j]=j;
    }

    for(m=nelements; m>1; m--)
    {
        int sum;
        int is=1;
        int js=0;
        result[nelements-m].distance=find_closest_pair(m, data, &is, &js);

        /* Save result */
        result[nelements-m].left=clusterid[is];
        result[nelements-m].right=clusterid[js];

        /* Fix the distances */
        sum = number[is] + number[js];
        for(j=0; j<js; j++)
        {
            data[js][j]=data[is][j]*number[is]+data[js][j]*number[js];
            data[js][j] /= sum;
        }
        for(j=js+1; j<is; j++)
        {
            data[j][js]=data[is][j]*number[is]+data[j][js]*number[js];
            data[j][js] /= sum;
        }
        for(j=is+1; j<m; j++)
        {
            data[j][js]=data[j][is]*number[is]+data[j][js]*number[js];
            data[j][js] /= sum;
        }

        for(j=0; j<is; j++)
        {
            data[is][j]=data[m-1][j];
        }
        for(j=is+1; j<m-1; j++)
        {
            data[j][is]=data[m-1][j];
        }

        /* Update number of elements in the clusters */
        number[js]=sum;
        number[is]=number[m-1];

        /* Update clusterids */
        clusterid[js]=m-nelements-1;
        clusterid[is]=clusterid[m-1];
    }

    for(i=1; i<=nelements-1; i++)
    {
        merge[n+i][1]=result[i-1].left;
        if(merge[n+i][1]<0)
            merge[n+i][1]=n-merge[n+i][1];
        else
            merge[n+i][1]=merge[n+i][1]+1;
        merge[n+i][2]=result[i-1].right;
        if(merge[n+i][2]<0)
            merge[n+i][2]=n-merge[n+i][2];
        else
            merge[n+i][2]=merge[n+i][2]+1;
        hgt[n+i]=result[i-1].distance;
    }

    for(i=0;i<nelements;i++)
    {
        delete [] data[i];
    }
    delete [] data;
    delete[] clusterid;
    delete[] result;
    delete[] number;

    // // print merge 和 hgt
    // std::cout << "hgt (heights):" << std::endl;
    // for (int i = 1; i <= 2 * n - 1; ++i) {
    //     std::cout << "hgt[" << i << "]: " << hgt[i] << std::endl;
    // }

    // std::cout << "merge (clusters):" << std::endl;
    // for (int i = 1; i <= 2 * n - 1; ++i) {
    //     std::cout << "merge[" << i << "]: " << merge[i][1] << ", " << merge[i][2] << std::endl;
    // }

/*
    FILE *fp;
    int n2=2*n-1;
    AnsiString tempstring;
    fp=fopen("after.txt","w");
    for(i=1;i<=n2;i++)
    {
        tempstring=hgt[i];
        fputs(tempstring.c_str(),fp);
        fputs(" ",fp);
    }
    fputs("\n",fp);
    for(i=1;i<=n2;i++)
    {
        tempstring=merge[i][1];
        fputs(tempstring.c_str(),fp);
        fputs(" ",fp);
    }
    fputs("\n",fp);
    for(i=1;i<=n2;i++)
    {
        tempstring=merge[i][2];
        fputs(tempstring.c_str(),fp);
        fputs(" ",fp);
    }
    fputs("\n",fp);

    fclose(fp);
*/
}//End of Function

//

void htree_average_s(double **_data,int n,int **merge,double *hgt)
{
  int i,j,n2,iter,*minwith,*active,ptmin,pm1,pm2;
  double **data,min,max,*minval;
  n2=2*n-1;
  //allocate
  data=dmatrix1(n2,n2);
  active=ivector1(n2);
  minwith=ivector1(n2);
  minval=dvector1(n2);
  //
  max=_data[1][1];
  for(i=1;i<n;i++){
    for(j=(i+1);j<=n;j++){
      if(_data[i][j]<max){
        max=_data[i][j];
      }
    }
  }
  max=max-1.0;
  for(i=1;i<=n;i++){
    merge[i][1]=merge[i][2]=0;
    hgt[i]=_data[i][i];
    active[i]=1;
    minval[i]=max;
    for(j=1;j<=n;j++){
      data[i][j]=_data[i][j];
      if((i != j) && (minval[i] < data[i][j])){
        minval[i]=data[i][j];
        minwith[i]=j;
      }
    }
  }
  for(j=(n+1);j<=n2;j++){
    active[j]=0;
  }
  iter=n+1;
  while(iter<=n2){
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(minval[i] > min){
          ptmin=i;
          min=minval[i];
        }
      }
    }
    hgt[iter]=min;
    if(ptmin < minwith[ptmin]){
      pm1=merge[iter][1]=ptmin;
      pm2=merge[iter][2]=minwith[ptmin];
    }else{
      pm1=merge[iter][1]=minwith[ptmin];
      pm2=merge[iter][2]=ptmin;
    }
    active[pm1]=0;
    active[pm2]=0;
    //
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        data[i][iter]=data[iter][i]=(data[i][pm1]+data[i][pm2])/2.0;
        //
        if((minwith[i] == pm1) || (minwith[i] == pm2)){
          minval[i]=max;
          for(j=1;j<iter;j++){
            if(active[j]){
              if(j != i){
                if(minval[i] < data[i][j]){
                  minval[i]=data[i][j];
                  minwith[i]=j;
                }
              }
            }
          }
        }
        //
        if(min < data[iter][i]){
          min=data[iter][i];
          ptmin=i;
        }
      }
    }
    minval[iter]=min;
    minwith[iter]=ptmin;
    active[iter]=1;
    iter++;
  }
  free_dmatrix1(data,n2);
  free_dvector1(minval);
  free_ivector1(minwith);
  free_ivector1(active);
}

//
// Hierachical UPGMA Linkage Tree
//

void htree_upgma_d(double **_data,int n,int **merge,double *hgt)
{
  int i,j,n2,iter,*minwith,*active,ptmin,pm1,pm2,*elements;
  double **data,min,max,*minval,c1,c2,cs;
  n2=2*n-1;
  //allocate
  data=dmatrix1(n2,n2);
  active=ivector1(n2);
  minwith=ivector1(n2);
  minval=dvector1(n2);
  elements=ivector1(n2);
  //
  max=_data[1][1];
  for(i=1;i<n;i++){
    for(j=(i+1);j<=n;j++){
      if(_data[i][j]>max){
        max=_data[i][j];
      }
    }
  }
  max=max+1.0;
  for(i=1;i<=n;i++){
    merge[i][1]=merge[i][2]=0;
    hgt[i]=_data[i][i];
    active[i]=1;
    elements[i]=1;
    minval[i]=max;
    for(j=1;j<=n;j++){
      data[i][j]=_data[i][j];
      if((i != j) && (minval[i] > data[i][j])){
        minval[i]=data[i][j];
        minwith[i]=j;
      }
    }
  }
  for(j=(n+1);j<=n2;j++){
    active[j]=0;
  }
  iter=n+1;
  while(iter<=n2){
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(minval[i] < min){
          ptmin=i;
          min=minval[i];
        }
      }
    }
    hgt[iter]=min;
    if(ptmin < minwith[ptmin]){
      pm1=merge[iter][1]=ptmin;
      pm2=merge[iter][2]=minwith[ptmin];
    }else{
      pm1=merge[iter][1]=minwith[ptmin];
      pm2=merge[iter][2]=ptmin;
    }
    active[pm1]=0;
    active[pm2]=0;
    //
    min=max;
    ptmin=0;
    c1=double(elements[pm1]);
    c2=double(elements[pm2]);
    cs=c1+c2;
    for(i=1;i<iter;i++){
      if(active[i]){

        data[i][iter]=data[iter][i]=(c1*data[i][pm1]+c2*data[i][pm2])/cs;
        //
        if((minwith[i] == pm1) || (minwith[i] == pm2)){
          minval[i]=max;
          for(j=1;j<iter;j++){
            if(active[j]){
              if(j != i){
                if(minval[i] > data[i][j]){
                  minval[i]=data[i][j];
                  minwith[i]=j;
                }
              }
            }
          }
        }
        //
        if(min > data[iter][i]){
          min=data[iter][i];
          ptmin=i;
        }
      }
    }
    minval[iter]=min;
    minwith[iter]=ptmin;
    active[iter]=1;
    elements[iter]=int(cs);
    iter++;
  }
  free_dmatrix1(data,n2);
  free_dvector1(minval);
  free_ivector1(minwith);
  free_ivector1(active);
  free_ivector1(elements);
}

//

void htree_upgma_s(double **_data,int n,int **merge,double *hgt)
{
  int i,j,n2,iter,*minwith,*active,ptmin,pm1,pm2,*elements;
  double **data,min,max,*minval,c1,c2,cs;
  n2=2*n-1;
  //allocate
  data=dmatrix1(n2,n2);
  active=ivector1(n2);
  minwith=ivector1(n2);
  minval=dvector1(n2);
  elements=ivector1(n2);
  //
  max=_data[1][1];
  for(i=1;i<n;i++){
    for(j=(i+1);j<=n;j++){
      if(_data[i][j]<max){
        max=_data[i][j];
      }
    }
  }
  max=max-1.0;
  for(i=1;i<=n;i++){
    merge[i][1]=merge[i][2]=0;
    hgt[i]=_data[i][i];
    active[i]=1;
    elements[i]=1;
    minval[i]=max;
    for(j=1;j<=n;j++){
      data[i][j]=_data[i][j];
      if((i != j) && (minval[i] < data[i][j])){
        minval[i]=data[i][j];
        minwith[i]=j;
      }
    }
  }
  for(j=(n+1);j<=n2;j++){
    active[j]=0;
  }
  iter=n+1;
  while(iter<=n2){
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        if(minval[i] > min){
          ptmin=i;
          min=minval[i];
        }
      }
    }
    hgt[iter]=min;
    if(ptmin < minwith[ptmin]){
      pm1=merge[iter][1]=ptmin;
      pm2=merge[iter][2]=minwith[ptmin];
    }else{
      pm1=merge[iter][1]=minwith[ptmin];
      pm2=merge[iter][2]=ptmin;
    }
    active[pm1]=0;
    active[pm2]=0;
    //
    c1=double(elements[pm1]);
    c2=double(elements[pm2]);
    cs=c1+c2;
    min=max;
    ptmin=0;
    for(i=1;i<iter;i++){
      if(active[i]){
        data[i][iter]=data[iter][i]=(c1*data[i][pm1]+c2*data[i][pm2])/cs;
        //
        if((minwith[i] == pm1) || (minwith[i] == pm2)){
          minval[i]=max;
          for(j=1;j<iter;j++){
            if(active[j]){
              if(j != i){
                if(minval[i] < data[i][j]){
                  minval[i]=data[i][j];
                  minwith[i]=j;
                }
              }
            }
          }
        }
        //
        if(min < data[iter][i]){
          min=data[iter][i];
          ptmin=i;
        }
      }
    }
    minval[iter]=min;
    minwith[iter]=ptmin;
    elements[iter]=int(cs);
    active[iter]=1;
    iter++;
  }
  free_dmatrix1(data,n2);
  free_dvector1(minval);
  free_ivector1(minwith);
  free_ivector1(active);
  free_ivector1(elements);
}



//void hctree_sort(double *alldataPtr, int *externalOrder, int *output_left_ptr, int *output_right_ptr, double *output_hgt, int *output_order, int nrow, int ncol, int orderType, int flipType)
void hctree_sort(const Rcpp::NumericMatrix& distance_matrix, int *externalOrder, int *output_left_ptr, int *output_right_ptr, double *output_hgt, int *output_order, int nrow, int ncol, int orderType, int flipType)
{
    
    //double **alldata = fun1dto2dArray(alldataPtr,nrow,ncol);
    // 建立 row-major 的 2D vector
    //std::vector<std::vector<double>> alldata(nrow+1, std::vector<double>(ncol+1));
    double** alldata = dmatrix1(nrow, ncol);
    // Copy ── 仍用 (i, j) 索引
    for (int i = 0; i < nrow; ++i)
      for (int j = 0; j < ncol; ++j)
        alldata[i+1][j+1] = distance_matrix(i, j);


    int **merge;
    double *hgt;
    merge=imatrix1(2*nrow-1,2);
    hgt=dvector1(2*nrow-1);

    if(orderType==0)
    {
        htree_single_d(alldata, nrow, merge, hgt);
    }
    else if(orderType==1)
    {
        htree_complete_d(alldata, nrow, merge, hgt);
    }
    else if(orderType==2)
    {
        htree_average_d(alldata, nrow, merge, hgt);
    }

    /*
    int* output_array = ellipse_sort_d(alldata, nrow, ncol, corr_row, recursive_n);

    for(int i=0; i<nrow; i++)
    {
        output_ptr[i] =  output_array[i];   
    }
    
    printf("result:\n");
    for(int i=1; i<=2*nrow-1; i++)
        printf("%d, %d, %d, %f\n",i,merge[i][1],merge[i][2],hgt[i]);
    printf("end\n");
    */
    List *ol;
    int *torder;
    torder=ivector1(nrow);
    if(flipType==0)
    {
      for(int i=1;i<=nrow;i++)
        torder[i]=i;
      ordering_external(torder,merge,nrow);
    }
    else if(flipType==1)
    {
      for(int i=1;i<=nrow;i++)
        torder[i]=externalOrder[i-1]+1;
      ordering_external(torder,merge,nrow);
    }
    else if(flipType==2)
    {
      ordering_Uncle_d(alldata,nrow,merge);
    }
    else if(flipType==3)
    {
      ordering_GrandPa_d(alldata,nrow,merge);
    }
    else
    {
      for(int i=1;i<=nrow;i++)
        torder[i]=i;
      ordering_external(torder,merge,nrow);      
    }
    
    
    ol=getOrderList(merge,2*nrow-1);
    int* obj_ordering=list2vector(ol);

    for(int i=0; i<ol->length(); i++)
      output_order[i] = obj_ordering[i+1]-1;

    for(int i=0; i<nrow-1; i++)
    {
      output_hgt[i] = hgt[i+1+nrow];
      output_left_ptr[i] = merge[i+1+nrow][1]-1;
      output_right_ptr[i] = merge[i+1+nrow][2]-1;
    }

    //for(int i=1;i<=nrow;i++)
        //printf("r2e: %d\n",torder[i]); 
    //for(int i=0; i<ol->length(); i++)
        //printf("%d\n",output_order[i]);   

    free_dmatrix1(alldata, nrow);
    free_imatrix1(merge, 2*nrow-1);
    free_dvector1(hgt);
    //free_ivector1(obj_ordering);
    //obj_ordering++;
    delete[] obj_ordering;
    free_ivector1(torder);
    delete ol;
        
}

#ifdef __cplusplus
}
#endif

// [[Rcpp::export]]
Rcpp::List hctree_sort_R(Rcpp::NumericMatrix distance_matrix, Rcpp::IntegerVector externalOrder, 
                         int orderType, int flipType) {
    int nrow = distance_matrix.nrow();
    int ncol = distance_matrix.ncol();

    // 準備輸出
    Rcpp::IntegerVector output_left_ptr(nrow - 1);
    Rcpp::IntegerVector output_right_ptr(nrow - 1);
    Rcpp::NumericVector output_hgt(nrow - 1);
    Rcpp::IntegerVector output_order(nrow);

    // 調用 C++ 的階層分群函數
    /*hctree_sort(&distance_matrix[0], externalOrder.begin(),
                output_left_ptr.begin(), output_right_ptr.begin(),
                output_hgt.begin(), output_order.begin(),
                nrow, ncol, orderType, flipType);
*/
    hctree_sort(distance_matrix, externalOrder.begin(),
                output_left_ptr.begin(), output_right_ptr.begin(),
                output_hgt.begin(), output_order.begin(),
                nrow, ncol, orderType, flipType);
    // 返回結果
    return Rcpp::List::create(
        Rcpp::Named("left") = output_left_ptr,
        Rcpp::Named("right") = output_right_ptr,
        Rcpp::Named("height") = output_hgt,
        Rcpp::Named("order") = output_order
    );
}
