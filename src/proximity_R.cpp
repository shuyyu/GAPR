#include <Rcpp.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <cstdio>

#ifdef __cplusplus
extern "C"
{
#endif
static double missing_value = -1.0E+10;

double** fun1dto2dArray_Proximity(double *input_ptr, int n, int p)
{
        double** output_ptr = new double* [n];
        for(int i=0; i<n; i++)
        {
            output_ptr[i] = new double[p];   
        }

        for(int i=0; i<n*p; i++)
        {
            int x = i/p;
            int y = i%p;
            output_ptr[x][y]=input_ptr[i];
        }

        //delete [] input_ptr;

        return output_ptr;
}

double* fun2dto1dArray_Proximity(double **input_ptr, int n, int p)
{
        double* output_ptr = new double[n*p];

        for(int i=0; i<n; i++)
        {
          for(int j=0; j<p; j++)
          {
            output_ptr[i*p+j]=input_ptr[i][j];  
          }
        }
/*

        for(int i=0; i<n; i++)
        {
            delete [] input_ptr[i];
        }
        delete [] input_ptr;
*/
        return output_ptr;
}

double* L2_distance(double** data, int no_row, int no_col, bool isDataContainMissingValue)
{
  int i, j, k;

  double** dis;
  dis = new double*[no_row];
  for(i=0;i<no_row;i++)
  {
    dis[i]=new double[no_row];
  }

  double s = 0;
  int missing_count = 0;

  if (isDataContainMissingValue)
  {
    for(i = 0; i<no_row; i++)
    {
      for(k = i+1; k<no_row; k++)
      {
        s = 0;
        missing_count = 0;
        for(j = 0; j<no_col; j++)
        {
          if (data[i][j]!=missing_value&&data[k][j]!=missing_value)
          {
            s += (data[i][j]-data[k][j])*(data[i][j]-data[k][j]);
            missing_count++;
          }
        }

        if (missing_count==0)
        {
          dis[i][k] = missing_value;
        }
        else
        {
          dis[i][k] = sqrt(s*no_row/(double) missing_count);
        }

        dis[k][i] = dis[i][k];
      }
      dis[i][i] = 0.0;
    }
  }
  else
  {
    for(i = 0; i<no_row; i++)
    {
      for(k = i+1; k<no_row; k++)
      {
        s = 0;
        for(j = 0; j<no_col; j++)
        {
          s += (data[i][j]-data[k][j])*(data[i][j]-data[k][j]);
        }
        dis[i][k] = sqrt(s);
        dis[k][i] = dis[i][k];
        //printf("dis%d %d: %f\n",i,k,dis[i][k]); 
      }
      dis[i][i] = 0.0;
    }
  }
  //printf("dis: %f\n",dis[1][0]); 
  return fun2dto1dArray_Proximity(dis, no_row, no_row);
}

double* L2_distanceT(double** data, int no_row, int no_col, bool isDataContainMissingValue){
  int i, j, k;

  double** dis;
  dis = new double*[no_col];
  for(i=0;i<no_col;i++)
  {
    dis[i]=new double[no_col];
  }
        double s = 0;
        int missing_count = 0;

        if (isDataContainMissingValue){

            for(i = 0; i<no_col; i++){
                for(k = i; k<no_col; k++){
                    s = 0;
                    missing_count = 0;
                    for(j = 0; j<no_row; j++){
                        if (data[j][i]!=missing_value&&data[j][k]!=missing_value){
                            s += (data[j][i]-data[j][k])*(data[j][i]-data[j][k]);
                            missing_count++;
                        }
                    }
                    if (missing_count==0){
                        dis[i][k] = missing_value;
                    }
                    else{
                        dis[i][k] = sqrt(s*no_row/(double) missing_count);
                    }

                    dis[k][i] = dis[i][k];
                }
                dis[i][i] = 0.0;
            }
        }
        else{
            for(i = 0; i<no_col; i++){
                for(k = i; k<no_col; k++){
                    s = 0;
                    for(j = 0; j<no_row; j++){
                        s += (data[j][i]-data[j][k])*(data[j][i]-data[j][k]);
                    }
                    dis[i][k] = sqrt(s);
                    dis[k][i] = dis[i][k];
                }
                dis[i][i] = 0.0;
            }
        }

        return fun2dto1dArray_Proximity(dis, no_col, no_col);

}

//// corr for row
double* corr_dataT(double** data, int data_obs, int data_var, bool isabs)
{
  ////for missing value
        int i, j, k;
        double miss_mean_i = 0;
        double miss_mean_j = 0;
        int new_count;
        //double sum;

        double** cov;
        cov = new double*[data_obs];
        for(i=0; i<data_obs; i++)
        {
          cov[i]=new double[data_obs];
        }
        double var_i, var_j, mycov;

        int pairwise;

        for(i = 0; i<data_obs; i++){
            for(j = i+1; j<data_obs; j++){

                new_count = 0;
                miss_mean_i = 0;
                miss_mean_j = 0;

                for(k = 0; k<data_var; k++){

                    if ((data[i][k]!=missing_value)&&(data[j][k]!=missing_value)){
                        new_count++;
                        miss_mean_i += data[i][k];
                        miss_mean_j += data[j][k];
                    }
                }

                miss_mean_i = miss_mean_i/new_count;
                miss_mean_j = miss_mean_j/new_count;

                mycov = 0.0;
                var_i = 0.0;
                var_j = 0.0;
                pairwise = 0;
                for(k = 0; k<data_var; k++){
                    if ((data[i][k]!=missing_value)&&(data[j][k]!=missing_value)){
                        mycov += (data[i][k]-miss_mean_i)*(data[j][k]-miss_mean_j);
                        var_i += (data[i][k]-miss_mean_i)*(data[i][k]-miss_mean_i);
                        var_j += (data[j][k]-miss_mean_j)*(data[j][k]-miss_mean_j);
                    }
                    else{
                        pairwise++;
                    }
                }

                if ((-0.00000001<var_i&&var_i<0.00000001)||(-0.00000001<var_j&&var_j<0.00000001)){
                    //cov[i][j]=missing_value;
                    cov[i][j] = 0;
                }
                else{
                    //corr
                    if (isabs){
                        cov[i][j] = fabs(mycov/sqrt(var_i*var_j));
                    }
                    else{
                        cov[i][j] = mycov/sqrt(var_i*var_j);
                    }

                }

                if (pairwise==data_var){
                    cov[i][j] = missing_value;
                }

                cov[j][i] = cov[i][j];

            }
            cov[i][i] = 1.0;
        }

        return fun2dto1dArray_Proximity(cov, data_obs, data_obs);
}

double* corr_data(double** data, int no_row, int no_col, bool isabs)
{
  ////missing value: pairwise
        int i, j, k;
        double miss_mean_i = 0;
        double miss_mean_j = 0;
        int new_count = 0;
        double var_i = 0.0;
        double var_j = 0.0;
        //double sum;
        double cov;
        double** cor;
        cor = new double*[no_col];
        for(i=0; i<no_col; i++)
        {
          cor[i]=new double[no_col];
        }
        int pairwise; //use when two rows/columns are pairwise missing

        for(i = 0; i<no_col; i++){
            for(j = i+1; j<no_col; j++){
                new_count = 0;
                miss_mean_i = 0.;
                miss_mean_j = 0.;
                cov = 0.;
                var_i = 0.;
                var_j = 0.;
                pairwise = 0;
                for(k = 0; k<no_row; k++){
                    if ((data[k][i]!=missing_value)&&(data[k][j]!=missing_value)){
                        new_count++;
                        miss_mean_i += data[k][i];
                        miss_mean_j += data[k][j];
                        cov += (data[k][i])*(data[k][j]);
                        var_i += (data[k][i])*(data[k][i]);
                        var_j += (data[k][j])*(data[k][j]);
                    }
                    else{
                        pairwise++;
                    }

                }
                miss_mean_i = miss_mean_i/(double) new_count;
                miss_mean_j = miss_mean_j/(double) new_count;
                cov = cov-new_count*miss_mean_i*miss_mean_j;
                var_i = var_i-new_count*miss_mean_i*miss_mean_i;
                var_j = var_j-new_count*miss_mean_j*miss_mean_j;

                if ((-0.00000001<var_i&&var_i<0.00000001)||(-0.00000001<var_j&&var_j<0.00000001)){
                    cor[i][j] = 0;
                    //cor[i][j]=missing_value;
                    //cor=cov(x,y)/V(x)V(y)
                    //
                }
                else{
                    //corr
                    if (isabs){
                        cor[i][j] = fabs(cov/sqrt(var_i*var_j));
                    }
                    else{
                        cor[i][j] = cov/sqrt(var_i*var_j);
                    }

                }
                if (pairwise==no_row){
                    cor[i][j] = missing_value;
                }

                cor[j][i] = cor[i][j];

            }
            cor[i][i] = 1.0;

        }

        return fun2dto1dArray_Proximity(cor, no_col, no_col);
}

double* kendallsTau_Row(double** data, int no_row, int no_col)
{
        int m, n, i, j;
        int con = 0;
        int dis = 0;
        int exx = 0;
        int exy = 0;

        double denomx;
        double denomy;
        double** tau;
        tau = new double*[no_row];
        for(i=0; i<no_row; i++)
        {
          tau[i]=new double[no_row];
        }
        
        for(m = 0; m<no_row; m++){
            for(n = (m+1); n<no_row; n++){

                con = 0;
                dis = 0;
                exx = 0;
                exy = 0;

                for(i = 0; i<no_col; i++){
                    for(j = 0; j<i; j++){
                        double x1 = data[m][i];
                        double x2 = data[m][j];
                        double y1 = data[n][i];
                        double y2 = data[n][j];

                        if (x1<x2&&y1<y2){
                            con++;
                        }
                        if (x1>x2&&y1>y2){
                            con++;
                        }
                        if (x1<x2&&y1>y2){
                            dis++;
                        }
                        if (x1>x2&&y1<y2){
                            dis++;
                        }
                        if (x1==x2&&y1!=y2){
                            exx++;
                        }
                        if (x1!=x2&&y1==y2){
                            exy++;
                        }
                    }
                }
                denomx = con+dis+exx;
                denomy = con+dis+exy;

                if (denomx==0){
                    tau[m][n] = 1.0;
                }
                else if (denomy==0){
                    tau[m][n] = 1.0;
                }
                else{
                    tau[m][n] = (double) (con-dis)/sqrt(denomx*denomy);
                }
                tau[n][m] = tau[m][n];
                tau[m][m] = 1;

            }
        }

        return fun2dto1dArray_Proximity(tau, no_row, no_row);
}

double* kendallsTau_Col(double** data, int no_row, int no_col)
{
        int m, n, i, j;
        int con = 0;
        int dis = 0;
        int exx = 0;
        int exy = 0;

        double denomx;
        double denomy;
        double **tau;
        tau = new double*[no_col];
        for(i=0; i<no_col; i++)
        {
          tau[i]=new double[no_col];
        }

        for(m = 0; m<no_col; m++){
            for(n = (m+1); n<no_col; n++){

                con = 0;
                dis = 0;
                exx = 0;
                exy = 0;

                for(i = 0; i<no_row; i++){
                    for(j = 0; j<i; j++){
                        double x1 = data[i][m];
                        double x2 = data[j][m];
                        double y1 = data[i][n];
                        double y2 = data[j][n];

                        if (x1<x2&&y1<y2){
                            con++;
                        }
                        if (x1>x2&&y1>y2){
                            con++;
                        }
                        if (x1<x2&&y1>y2){
                            dis++;
                        }
                        if (x1>x2&&y1<y2){
                            dis++;
                        }
                        if (x1==x2&&y1!=y2){
                            exx++;
                        }
                        if (x1!=x2&&y1==y2){
                            exy++;
                        }
                    }
                }
                denomx = con+dis+exx;
                denomy = con+dis+exy;

                if (denomx==0){
                    tau[m][n] = 1.0;
                }
                else if (denomy==0){
                    tau[m][n] = 1.0;
                }
                else{
                    tau[m][n] = (double) (con-dis)/sqrt(denomx*denomy);
                }
                tau[n][m] = tau[m][n];
                tau[m][m] = 1;

            }
        }

        return fun2dto1dArray_Proximity(tau, no_col, no_col);
}

    /** Swap two positions.
     @param a    Array.
     @param i    Line number.
     @param j    Line number.
     */

void swap(double* a, int* order, int i, int j){
        double T;

        T = a[i];
        a[i] = a[j];
        a[j] = T;
        int t;

        t = order[i];
        order[i] = order[j];
        order[j] = t;
}

void QuickSort(double* a, int* order, int lo0, int hi0){

        int lo = lo0;
        int hi = hi0;
        double mid;

        if (hi0>lo0){
            //Arbitrarily establishing partition element as the midpoint of the array.
            mid = a[(lo0+hi0)/2];

            //loop through the array until indices cross
            while (lo<=hi){
                //find the first element that is greater than or equal to the partition element starting from the left Index.
                while ((lo<hi0)&&(a[lo]<mid)){
                    ++lo;
                }
                //find an element that is smaller than or equal to the partition element starting from the right Index.
                while ((hi>lo0)&&(a[hi]>mid)){
                    --hi;
                }
                //if the indexes have not crossed, swap
                if (lo<=hi){
                    swap(a, order, lo, hi);
                    ++lo;
                    --hi;
                }
            }

            //If the right index has not reached the left side of array must now sort the left partition.
            if (lo0<hi){
                QuickSort(a, order, lo0, hi);
            }

            //If the left index has not reached the right side of array must now sort the right partition.
            if (lo<hi0){
                QuickSort(a, order, lo, hi0);
            }

        }
}



void sort(double* a, int len, int* order){
        QuickSort(a, order, 0, len-1);
}

double* getrank(double* array, int no_obs, bool double_id){

        int i, j;
        double* A = new double[no_obs];
        for(i = 0; i<no_obs; i++){
          A[i] = array[i];
        }
        int* order = new int[no_obs];

        for (i = 0; i<no_obs; i++){
            order[i] = i;
        }
        sort(A, no_obs, order);

        int* rank = new int[no_obs];
        double* rank_tie = new double[no_obs];

        double sum = 1;
        double t = A[0];
        int count = 1;

        rank[order[0]] = 0;
        rank_tie[order[0]] = (double) (1);

        for (i = 1; i<no_obs; i++){

            rank[order[i]] = i;

            if (t==A[i]){

                sum = sum+(i+1);
                count++;

                //last element
                if (i==(no_obs-1)){
                    for (j = 0; j<count; j++){
                        rank_tie[order[i-j]] = (double) sum/(double) count;
                    }
                    rank_tie[order[i]] = (double) sum/(double) count;
                }

            }
            else{ //until different

                t = A[i];
                if (count!=1){
                    for (j = 0; j<count; j++){
                        rank_tie[order[i-1-j]] = (double) sum/(double) count;
                    }
                }
                ////this element 
                rank_tie[order[i]] = (double) (i+1);

                sum = i+1;
                count = 1;
            }
        }
        return rank_tie;
}

double** getRankRow_tie(double** data, int no_row, int no_col){

        int i, j;
        double* vector = new double[no_col];
        double** rank;
        rank = new double*[no_row];
        for(i=0; i<no_row; i++)
        {
          rank[i]=new double[no_col];
        }
        double* rank_tie;

        //mySorting mysort;
        bool double_id = true;
        //mysort = new mySorting();

        for(i = 0; i<no_row; i++){

            for(j = 0; j<no_col; j++){
                vector[j] = data[i][j];
            }

            rank_tie = getrank(vector, no_col, double_id); //double[] array, boolean copyArray)
            for(j = 0; j<no_col; j++){
                rank[i][j] = rank_tie[j];
            }
        }

        return rank;

}

double** getRankCol_tie(double** data, int no_row, int no_col){

        int i, j;
        double* vector = new double[no_row];
        //double[][] rank = new double[no_row][no_col];
        double** rank;
        rank = new double*[no_row];
        for(i=0; i<no_row; i++)
        {
          rank[i]=new double[no_col];
        }
        //mySorting mysort;
        bool double_id = true;
        double* rank_tie;
        //mysort = new mySorting();

        for(j = 0; j<no_col; j++){

            for(i = 0; i<no_row; i++){
                vector[i] = data[i][j];
            }

           // mysort = new mySorting();
            rank_tie = getrank(vector, no_row, double_id); //double[] array, boolean copyArray)
            for(i = 0; i<no_row; i++){
                rank[i][j] = rank_tie[i];
            }
            //System.out.println("last row rank: " + rank[no_row-1][j]);
        }

        return rank;

}

double* atan_data(double** data, int no_row, int no_col){

        ////missing value: pairwise
        int i, j, k;
        double miss_mean_i = 0;
        double miss_mean_j = 0;
        int new_count = 0;
        double var_i = 0.0;
        double var_j = 0.0;
        //double sum;
        double cov;
        double cor;
        double answer;
        double** corM;
        corM = new double*[no_col];
        for(i=0; i<no_col; i++)
        {
          corM[i]=new double[no_col];
        }
        int pairwise; //use when two rows/columns are pairwise missing

        for(i = 0; i<no_col; i++){
            for(j = i+1; j<no_col; j++){
                new_count = 0;
                miss_mean_i = 0.;
                miss_mean_j = 0.;
                cov = 0.;
                var_i = 0.;
                var_j = 0.;
                pairwise = 0;
                for(k = 0; k<no_row; k++){
                    if ((data[k][i]!=missing_value)&&(data[k][j]!=missing_value)){
                        new_count++;
                        miss_mean_i += data[k][i];
                        miss_mean_j += data[k][j];
                        cov += (data[k][i])*(data[k][j]);
                        var_i += (data[k][i])*(data[k][i]);
                        var_j += (data[k][j])*(data[k][j]);
                    }
                    else{
                        pairwise++;
                    }

                }
                miss_mean_i = miss_mean_i/(double) new_count;
                miss_mean_j = miss_mean_j/(double) new_count;
                cov = cov-new_count*miss_mean_i*miss_mean_j;
                var_i = var_i-new_count*miss_mean_i*miss_mean_i;
                var_j = var_j-new_count*miss_mean_j*miss_mean_j;

                if ((-0.00000001<var_i&&var_i<0.00000001)||(-0.00000001<var_j&&var_j<0.00000001)){
                    corM[i][j] = 0;
                    //cor[i][j]=missing_value;
                    //cor=cov(x,y)/V(x)V(y)
                    //
                }
                else{
                    //corr
                    answer = 2*cov/fabs(var_i-var_j);
                    cor = cov/sqrt(var_i*var_j);
                    corM[i][j] = 2*fabs(cor)*atan(answer)/M_PI;

                }
                if (pairwise==no_row){
                    corM[i][j] = missing_value;
                }

                corM[j][i] = corM[i][j];

            }
            corM[i][i] = 1.0;

        }

        return fun2dto1dArray_Proximity(corM, no_col, no_col);
}

//// corr for row
double* atan_dataT(double** data, int data_obs, int data_var){
        ////for missing value

        int i, j, k;
        double miss_mean_i = 0;
        double miss_mean_j = 0;
        int new_count;
        //double sum;
        double** covM;
        covM = new double*[data_obs];
        for(i=0; i<data_obs; i++)
        {
          covM[i]=new double[data_obs];
        }
        double var_i, var_j, mycov;
        double cor;
        double answer;

        int pairwise;

        for(i = 0; i<data_obs; i++){
            for(j = i+1; j<data_obs; j++){

                new_count = 0;
                miss_mean_i = 0;
                miss_mean_j = 0;

                for(k = 0; k<data_var; k++){

                    if ((data[i][k]!=missing_value)&&(data[j][k]!=missing_value)){
                        new_count++;
                        miss_mean_i += data[i][k];
                        miss_mean_j += data[j][k];
                    }
                }

                miss_mean_i = miss_mean_i/new_count;
                miss_mean_j = miss_mean_j/new_count;

                mycov = 0.0;
                var_i = 0.0;
                var_j = 0.0;
                pairwise = 0;
                for(k = 0; k<data_var; k++){
                    if ((data[i][k]!=missing_value)&&(data[j][k]!=missing_value)){
                        mycov += (data[i][k]-miss_mean_i)*(data[j][k]-miss_mean_j);
                        var_i += (data[i][k]-miss_mean_i)*(data[i][k]-miss_mean_i);
                        var_j += (data[j][k]-miss_mean_j)*(data[j][k]-miss_mean_j);
                    }
                    else{
                        pairwise++;
                    }
                }

                if ((-0.00000001<var_i&&var_i<0.00000001)||(-0.00000001<var_j&&var_j<0.00000001)){
                    //cov[i][j]=missing_value;
                    covM[i][j] = 0;
                }
                else{
                    //corr
                    answer = 2*mycov/fabs(var_i-var_j);
                    cor = mycov/sqrt(var_i*var_j);
                    covM[i][j] = 2*fabs(cor)*atan(answer)/M_PI;
                }

                if (pairwise==data_var){
                    covM[i][j] = missing_value;
                }

                covM[j][i] = covM[i][j];

            }
            covM[i][i] = 1.0;
        }

        return fun2dto1dArray_Proximity(covM, data_obs, data_obs);
}

double* CityBlock(double** data, int no_row, int no_col, bool isDataContainMissingValue){

        int i, j, k;

        double** dis;
        dis = new double*[no_row];
        for(i=0;i<no_row;i++)
        {
          dis[i]=new double[no_row];
        }
        double s = 0;
        int missing_count = 0;

        if (isDataContainMissingValue){

            for(i = 0; i<no_row; i++){
                for(k = i+1; k<no_row; k++){
                    s = 0;
                    missing_count = 0;
                    for(j = 0; j<no_col; j++){
                        if (data[i][j]!=missing_value&&data[k][j]!=missing_value){
                            s += fabs(data[i][j]-data[k][j]);
                            missing_count++;
                        }
                    }
                    //very slow by using pow
                    //dis[i][k]=Math.pow(s, 1/2.0);
                    //d * sqrt(all/common)
                    //dis[i][k]=Math.sqrt(s*no_col/(double)missing_count);
                    if (missing_count==0){
                        dis[i][k] = missing_value;
                    }
                    else{
                        dis[i][k] = sqrt(s*no_row/(double) missing_count);
                    }

                    dis[k][i] = dis[i][k];
                }
                dis[i][i] = 0.0;
            }

        }
        else{
            for(i = 0; i<no_row; i++){
                for(k = i+1; k<no_row; k++){
                    s = 0;
                    for(j = 0; j<no_col; j++){
                        s += fabs(data[i][j]-data[k][j]);
                    }
                    dis[i][k] = sqrt(s);
                    dis[k][i] = dis[i][k];
                }
                dis[i][i] = 0.0;
            }

        }
        return fun2dto1dArray_Proximity(dis, no_row, no_row);
}

double* CityBlockT(double** data, int no_row, int no_col, bool isDataContainMissingValue){
        int i, j, k;

        double** dis;
        dis = new double*[no_col];
        for(i=0;i<no_col;i++)
        {
          dis[i]=new double[no_col];
        }
        double s = 0;
        int missing_count = 0;

        if (isDataContainMissingValue){

            for(i = 0; i<no_col; i++){
                for(k = i; k<no_col; k++){
                    s = 0;
                    missing_count = 0;
                    for(j = 0; j<no_row; j++){
                        if (data[j][i]!=missing_value&&data[j][k]!=missing_value){
                            s += fabs(data[j][i]-data[j][k]);
                            missing_count++;
                        }
                    }
                    if (missing_count==0){
                        dis[i][k] = missing_value;
                    }
                    else{
                        dis[i][k] = sqrt(s*no_row/(double) missing_count);
                    }

                    dis[k][i] = dis[i][k];
                }
                dis[i][i] = 0.0;
            }
        }
        else{
            for(i = 0; i<no_col; i++){
                for(k = i; k<no_col; k++){
                    s = 0;
                    for(j = 0; j<no_row; j++){
                        s += fabs(data[j][i]-data[j][k]);
                    }
                    dis[i][k] = sqrt(s);
                    dis[k][i] = dis[i][k];
                }
                dis[i][i] = 0.0;
            }
        }
        return fun2dto1dArray_Proximity(dis, no_col, no_col);
}

double* Ucorr_data(double** data, int no_row, int no_col, bool isabs)
{
        ////missing value: pairwise
        int i, j, k;
        //int new_count = 0;
        double var_i = 0.0;
        double var_j = 0.0;
        //double sum;
        double cov;
        double** cor;
        cor = new double*[no_col];
        for(i=0; i<no_col; i++)
        {
          cor[i]=new double[no_col];
        }
        int pairwise; //use when two rows/columns are pairwise missing

        for(i = 0; i<no_col; i++){
            for(j = i+1; j<no_col; j++){
                //new_count = 0;
                cov = 0.;
                var_i = 0.;
                var_j = 0.;
                pairwise = 0;
                for(k = 0; k<no_row; k++){
                    if ((data[k][i]!=missing_value)&&(data[k][j]!=missing_value)){
                        //new_count++;
                        cov += (data[k][i])*(data[k][j]);
                        var_i += (data[k][i])*(data[k][i]);
                        var_j += (data[k][j])*(data[k][j]);
                    }
                    else{
                        pairwise++;
                    }

                }

                if ((-0.00000001<var_i&&var_i<0.00000001)||(-0.00000001<var_j&&var_j<0.00000001)){
                    cor[i][j] = 0;
                }
                else{
                    //corr
                    if (isabs){
                        cor[i][j] = fabs(cov/sqrt(var_i*var_j));
                    }
                    else{
                        cor[i][j] = cov/sqrt(var_i*var_j);
                    }

                }
                if (pairwise==no_row){
                    cor[i][j] = missing_value;
                }

                cor[j][i] = cor[i][j];

            }
            cor[i][i] = 1.0;

        }
        return fun2dto1dArray_Proximity(cor, no_col, no_col);
}

//// corr for row
double* Ucorr_dataT(double** data, int data_obs, int data_var, bool isabs)
{
        ////for missing value
        int i, j, k;
        //int new_count;
        //double sum;

        double** cov;
        cov = new double*[data_obs];
        for(i=0; i<data_obs; i++)
        {
          cov[i]=new double[data_obs];
        }
        double var_i, var_j, mycov;

        int pairwise;

        for(i = 0; i<data_obs; i++){
            for(j = i+1; j<data_obs; j++){

                //new_count = 0;
                mycov = 0.0;
                var_i = 0.0;
                var_j = 0.0;
                pairwise = 0;
                for(k = 0; k<data_var; k++){
                    if ((data[i][k]!=missing_value)&&(data[j][k]!=missing_value)){
                        mycov += (data[i][k]*data[j][k]);
                        var_i += (data[i][k]*data[i][k]);
                        var_j += (data[j][k]*data[j][k]);
                    }
                    else{
                        pairwise++;
                    }
                }

                if ((-0.00000001<var_i&&var_i<0.00000001)||(-0.00000001<var_j&&var_j<0.00000001)){
                    //cov[i][j]=missing_value;
                    cov[i][j] = 0;
                }
                else{
                    //corr
                    if (isabs){
                        cov[i][j] = fabs(mycov/sqrt(var_i*var_j));
                    }
                    else{
                        cov[i][j] = mycov/sqrt(var_i*var_j);
                    }

                }

                if (pairwise==data_var){
                    cov[i][j] = missing_value;
                }

                cov[j][i] = cov[i][j];

            }
            cov[i][i] = 1.0;
        }
        return fun2dto1dArray_Proximity(cov, data_obs, data_obs);
}

double* cov_data(double** data, int no_row, int no_col, bool isabs)
{
        ////missing value: pairwise
        int i, j, k;
        double miss_mean_i = 0;
        double miss_mean_j = 0;
        int new_count = 0;
        double var_i = 0.0;
        double var_j = 0.0;
        //double sum;
        double cov;
        double** cor;
        cor = new double*[no_col];
        for(i=0; i<no_col; i++)
        {
          cor[i]=new double[no_col];
        }
        int pairwise; //use when two rows/columns are pairwise missing

        for(i=0; i<no_col; i++){
            for(j=i; j<no_col; j++){
                new_count = 0;
                miss_mean_i = 0.;
                miss_mean_j = 0.;
                cov = 0.;
                var_i = 0.;
                var_j = 0.;
                pairwise = 0;
                for(k = 0; k<no_row; k++){
                    if ((data[k][i]!=missing_value)&&(data[k][j]!=missing_value)){
                        new_count++;
                        miss_mean_i += data[k][i];
                        miss_mean_j += data[k][j];
                        cov += (data[k][i])*(data[k][j]);
                        var_i += (data[k][i])*(data[k][i]);
                        var_j += (data[k][j])*(data[k][j]);
                    }
                    else{
                        pairwise++;
                    }

                }
                miss_mean_i = miss_mean_i/(double) new_count;
                miss_mean_j = miss_mean_j/(double) new_count;
                cov = cov-new_count*miss_mean_i*miss_mean_j;
                var_i = var_i-new_count*miss_mean_i*miss_mean_i;
                var_j = var_j-new_count*miss_mean_j*miss_mean_j;

                cor[i][j] = cov;

                if (pairwise==no_row){
                    cor[i][j] = missing_value;
                }
                else  //maokao fixed..(20120829)
                {
                  cor[i][j] = cov/(new_count-1);
                }

                cor[j][i] = cor[i][j];
            }

        }
        return fun2dto1dArray_Proximity(cor, no_col, no_col);
}

//// corr for row
double* cov_dataT(double** data, int data_obs, int data_var, bool isabs)
{
        ////for missing value
        int i, j, k;
        double miss_mean_i = 0;
        double miss_mean_j = 0;
        int new_count;
        //double sum;
        double** cov;
        cov = new double*[data_obs];
        for(i=0; i<data_obs; i++)
        {
          cov[i]=new double[data_obs];
        }
        //double var_i, var_j, mycov;
        double mycov;
        int pairwise;

        for(i = 0; i<data_obs; i++){
            for(j = i; j<data_obs; j++){

                new_count = 0;
                miss_mean_i = 0;
                miss_mean_j = 0;

                for(k = 0; k<data_var; k++){

                    if ((data[i][k]!=missing_value)&&(data[j][k]!=missing_value)){
                        new_count++;
                        miss_mean_i += data[i][k];
                        miss_mean_j += data[j][k];
                    }
                }

                miss_mean_i = miss_mean_i/new_count;
                miss_mean_j = miss_mean_j/new_count;

                mycov = 0.0;
                //var_i = 0.0;
                //var_j = 0.0;
                pairwise = 0;
                for(k = 0; k<data_var; k++){
                    if ((data[i][k]!=missing_value)&&(data[j][k]!=missing_value)){
                        mycov += (data[i][k]-miss_mean_i)*(data[j][k]-miss_mean_j);
                        //var_i += (data[i][k]-miss_mean_i)*(data[i][k]-miss_mean_i);
                        //var_j += (data[j][k]-miss_mean_j)*(data[j][k]-miss_mean_j);
                    }
                    else{
                        pairwise++;
                    }
                }

                cov[i][j] = mycov;


                if (pairwise==data_var){
                    cov[i][j] = missing_value;
                }
                else  //maokao fixed..(20120829)
                {
                  cov[i][j] = mycov/(new_count-1);
                }

                cov[j][i] = cov[i][j];

            }

        }
        return fun2dto1dArray_Proximity(cov, data_obs, data_obs);
}


    ///////////////////////////////////////////////////////////////////
    //                                                               //
    //              BinaryProximityMeasure for Col                   //
    //                                                               //
    ///////////////////////////////////////////////////////////////////

double* BinaryProximityMeasure_Col(double** data, int no_row, int no_col, int Measure_ids){

        int i, j, k;
        double** cor;
        cor = new double*[no_col];
        for(i=0; i<no_col; i++)
        {
          cor[i]=new double[no_col];
        }
        int table[2][2];

        double numerator = 0;
        double denominator = 0;

        ////missing value: pairwise
        double a = 0, b = 0, c = 0, d = 0;

        //    1  0
        //1 | a  b|
        //0 | c  d|

        int pairwise = 0; //use when two rows/columns are pairwise missing

        for(i = 0; i<no_col; i++){
            for(j = i+1; j<no_col; j++){

                table[0][0] = 0; //d
                table[0][1] = 0; //c
                table[1][0] = 0; //b
                table[1][1] = 0; //a
                pairwise = 0;
                for(k = 0; k<no_row; k++){
                    if ((data[k][i]!=missing_value)&&(data[k][j]!=missing_value)){
                        table[(int) data[k][i]][(int) data[k][j]]++;

                    }
                    else{
                        pairwise++;
                    }
                }
                d = (double) table[0][0];
                c = (double) table[0][1];
                b = (double) table[1][0];
                a = (double) table[1][1];


                if(Measure_ids==20){  //if(Measure_ids.startsWith("Hamman")){
                    numerator = a+d-b-c;
                    denominator = a+b+c+d;
                }
                else if(Measure_ids==21){  //else if(Measure_ids.startsWith("Jaccard")){
                    numerator = a;
                    denominator = a+b+c;
                }
                else if(Measure_ids==22){  //else if(Measure_ids.startsWith("Phi")){
                    numerator = a*d-b*c;
                    denominator = sqrt((a+b)*(a+c)*(d+b)*(d+c));
                }
                else if(Measure_ids==23){  //else if(Measure_ids.startsWith("Rao")){
                    numerator = a;
                    denominator = a+b+c+d;
                }
                else if(Measure_ids==24){  //else if(Measure_ids.startsWith("Rogers")){
                    numerator = a+d;
                    denominator = a+2*b+2*c+d;
                }
                else if(Measure_ids==25){  //else if(Measure_ids.startsWith("Simple")){
                    numerator = a+d;
                    denominator = a+b+c+d;
                }
                else if(Measure_ids==26){  //else if(Measure_ids.startsWith("Sneath")){
                    numerator = a;
                    denominator = a+2*b+2*c;
                }
                else if(Measure_ids==27){  //else if(Measure_ids.startsWith("Yule")){
                    numerator = a*d-b*c;
                    denominator = a*d+b*c;
                }


                if (-0.00000001<denominator&&denominator<0.00000001){
                    if (-0.00000001<numerator&&numerator<0.00000001){
                        cor[i][j] = 0;
                    }
                    else{
                        cor[i][j] = missing_value;
                    }
                    //cor[i][j]=0;
                    //cor=cov(x,y)/V(x)V(y)
                }
                else{
                    cor[i][j] = numerator/denominator;
                }

                if (pairwise==no_row){
                    cor[i][j] = missing_value;
                }

                cor[j][i] = cor[i][j];

            }
            cor[i][i] = 1.0;
        }

        return fun2dto1dArray_Proximity(cor, no_col, no_col);

}

    ///////////////////////////////////////////////////////////////////
    //                                                               //
    //              BinaryProximityMeasure for Row                   //
    //                                                               //
    ///////////////////////////////////////////////////////////////////
double* BinaryProximityMeasure_Row(double** data, int no_row, int no_col, int Measure_ids){

        int i, j, k;
        double** cor;
        cor = new double*[no_row];
        for(i=0; i<no_row; i++)
        {
          cor[i]=new double[no_row];
        }
        int table[2][2];

        double numerator = 0;
        double denominator = 0;

        ////missing value: pairwise
        
        double a = 0, b = 0, c = 0, d = 0;

        int pairwise = 0; //use when two rows/columns are pairwise missing

        for(i = 0; i<no_row; i++){
            for(j = i+1; j<no_row; j++){

                table[0][0] = 0; //d
                table[0][1] = 0; //c
                table[1][0] = 0; //b
                table[1][1] = 0; //a
                pairwise = 0;
                for(k = 0; k<no_col; k++){
                    if ((data[i][k]!=missing_value)&&(data[j][k]!=missing_value)){
                        table[(int) data[i][k]][(int) data[j][k]]++;
                    }
                    else{
                        pairwise++;
                    }
                }
                d = (double) table[0][0];
                c = (double) table[0][1];
                b = (double) table[1][0];
                a = (double) table[1][1];

                /*
                v_Measure.addElement("Hamman: a/(b+c)");
                v_Measure.addElement("Jaccard: a/(a+b+c)");
                v_Measure.addElement("Phi: (ad-bc)/sqrt{(a+b)(a+c)(b+d)(c+d)}");
                v_Measure.addElement("Rao: a/(a+b+c+d)");
                v_Measure.addElement("Rogers: (a+d)/(a+2b+2c+d)");
                v_Measure.addElement("Simple Match: (a+d)/(a+b+c+d)");
                v_Measure.addElement("Sneath: a/(a+2b+2c)");
                v_Measure.addElement("Yule: (ad-bc)/(ad+bc)");
                */


                if(Measure_ids==20){  //if(Measure_ids.startsWith("Hamman")){
                    numerator = a+d-b-c;
                    denominator = a+b+c+d;
                }
                else if(Measure_ids==21){  //else if(Measure_ids.startsWith("Jaccard")){
                    numerator = a;
                    denominator = a+b+c;
                }
                else if(Measure_ids==22){  //else if(Measure_ids.startsWith("Phi")){
                    numerator = a*d-b*c;
                    denominator = sqrt((a+b)*(a+c)*(d+b)*(d+c));
                }
                else if(Measure_ids==23){  //else if(Measure_ids.startsWith("Rao")){
                    numerator = a;
                    denominator = a+b+c+d;
                }
                else if(Measure_ids==24){  //else if(Measure_ids.startsWith("Rogers")){
                    numerator = a+d;
                    denominator = a+2*b+2*c+d;
                }
                else if(Measure_ids==25){  //else if(Measure_ids.startsWith("Simple")){
                    numerator = a+d;
                    denominator = a+b+c+d;
                }
                else if(Measure_ids==26){  //else if(Measure_ids.startsWith("Sneath")){
                    numerator = a;
                    denominator = a+2*b+2*c;
                }
                else if(Measure_ids==27){  //else if(Measure_ids.startsWith("Yule")){
                    numerator = a*d-b*c;
                    denominator = a*d+b*c;
                }

                if (-0.00000001<denominator&&denominator<0.00000001){
                    if (-0.00000001<numerator&&numerator<0.00000001){
                        cor[i][j] = 0;
                    }
                    else{

                        //cor[i][j]=0;
                        cor[i][j] = missing_value;
                        //cor=cov(x,y)/V(x)V(y)
                    }
                }
                else{
                    cor[i][j] = numerator/denominator;
                }

                if (pairwise==no_col){
                    cor[i][j] = missing_value;
                }

                cor[j][i] = cor[i][j];
            }
            cor[i][i] = 1.0;
            //printf("dis: %f\n",cor[i][j]); 
        }
        
        return fun2dto1dArray_Proximity(cor, no_row, no_row);

}


void computeProximity(double *alldataPtr, double *outputProx, int nrow, int ncol, int proxType, int side, int isContainMissingValue)
{
    // printf("start to compute proximity.\n"); 

    double **alldata = fun1dto2dArray_Proximity(alldataPtr,nrow,ncol);
    double *proxArray = nullptr;

    bool isDataContainMissingValue;
    if(isContainMissingValue==0)
      isDataContainMissingValue=false;
    else
      isDataContainMissingValue=true;

    if(proxType==0)
    {
      if(side==0) //row
      {
        //proxArray = new double[nrow*nrow];
        proxArray = L2_distance(alldata, nrow, ncol, isDataContainMissingValue);
      }     
      else
      {
        //proxArray = new double[ncol*ncol];
        proxArray = L2_distanceT(alldata, nrow, ncol, isDataContainMissingValue);
      }
    }
    else if(proxType==1)
    {
      if(side==0) //row
      {
        proxArray = corr_dataT(alldata, nrow, ncol, false);
      }     
      else
      {
        proxArray = corr_data(alldata, nrow, ncol, false);
      }
    }
    else if(proxType==2)
    {
      if(side==0) //row
      {
        proxArray = kendallsTau_Row(alldata, nrow, ncol);
      }     
      else
      {
        proxArray = kendallsTau_Col(alldata, nrow, ncol);
      }
    }
    else if(proxType==3)
    {
      if(side==0) //row
      {
        double** rankData = getRankRow_tie(alldata, nrow, ncol);
        proxArray = corr_dataT(rankData, nrow, ncol, false);
      }     
      else
      {
        double** rankData = getRankCol_tie(alldata, nrow, ncol);
        proxArray = corr_data(rankData, nrow, ncol, false);
      }
    }
    else if(proxType==4)
    {
      if(side==0) //row
      {
        proxArray = atan_dataT(alldata, nrow, ncol);
      }     
      else
      {
        proxArray = atan_data(alldata, nrow, ncol);
      }
    }
    else if(proxType==5)
    {
      if(side==0) //row
      {
        proxArray = CityBlock(alldata, nrow, ncol, isDataContainMissingValue);
      }     
      else
      {
        proxArray = CityBlockT(alldata, nrow, ncol, isDataContainMissingValue);
      }
    }
    else if(proxType==6)
    {
      if(side==0) //row
      {
        proxArray = corr_dataT(alldata, nrow, ncol, true);
      }     
      else
      {
        proxArray = corr_data(alldata, nrow, ncol, true);
      }
    }
    else if(proxType==7)
    {
      if(side==0) //row
      {
        proxArray = Ucorr_dataT(alldata, nrow, ncol, false);
      }     
      else
      {
        proxArray = Ucorr_data(alldata, nrow, ncol, false);
      }
    }
    else if(proxType==8)
    {
      if(side==0) //row
      {
        proxArray = Ucorr_dataT(alldata, nrow, ncol, true);
      }     
      else
      {
        proxArray = Ucorr_data(alldata, nrow, ncol, true);
      }
    }
    else if(proxType==9)
    {
      if(side==0) //row
      {
        proxArray = cov_dataT(alldata, nrow, ncol, false);
      }     
      else
      {
        proxArray = cov_data(alldata, nrow, ncol, false);
      }
    }
    else if(proxType>=20 && proxType<=27)
    {
      if(side==0) //row
      {
        proxArray = BinaryProximityMeasure_Row(alldata, nrow, ncol, proxType);
      }     
      else
      {
        proxArray = BinaryProximityMeasure_Col(alldata, nrow, ncol, proxType);
      }
    }

    if(proxArray){
        if(side==0) //row
        { 
          for(int i=0; i<(nrow*nrow); i++)
              outputProx[i] = proxArray[i];
        }
        else
        { 
          for(int i=0; i<(ncol*ncol); i++)
              outputProx[i] = proxArray[i];
        }
    }

    for(int i=0; i<nrow; i++){
      delete[] alldata[i];
    }
    delete[] alldata;

    if(proxArray)
        delete[] proxArray;


    //printf("outputProx: %f\n",outputProx[1]); 

}

#ifdef __cplusplus
}
#endif
 
// [[Rcpp::export]]
Rcpp::NumericMatrix computeProximity_R(Rcpp::NumericMatrix data, int proxType, int side, int isContainMissingValue) {
    int nrow = data.nrow();
    int ncol = data.ncol();

    // Convert R data to 1D array
    double* alldata = new double[nrow * ncol];
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            alldata[i * ncol + j] = data(i, j);
        }
    }

    // The output array size is determined by side
    int outputSize = (side == 0) ? (nrow * nrow) : (ncol * ncol);
    double* outputProx = new double[outputSize];

    // Calling C++ function
    computeProximity(alldata, outputProx, nrow, ncol, proxType, side, isContainMissingValue);

    // Convert the output to R matrix
    Rcpp::NumericMatrix result((side == 0) ? nrow : ncol, (side == 0) ? nrow : ncol);
    for (int i = 0; i < result.nrow(); ++i) {
        for (int j = 0; j < result.ncol(); ++j) {
            result(i, j) = outputProx[i * result.ncol() + j];
        }
    }

    // Clean up the memory
    delete[] alldata;
    delete[] outputProx;

    return result;
}

