#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#ifdef __cplusplus
extern "C"
{
#endif

double** fun1dto2dArray_r2e(double *input_ptr, int n, int p)
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

        return output_ptr;
    }

int* ivector1_r2e(int length)
{
  int *out;
  out=new int[length];
  //out--;
  return out;
}

void free_ivector1_r2e(int* vec)
{
  //vec++;
  delete[] vec;
}

double* dvector1_r2e(int length)
{
  double *out;
  out=new double[length];
  //out--;
  return out;
}

void free_dvector1_r2e(double* vec)
{
  //vec++;
  delete[] vec;
}

int** imatrix1_r2e(int nrow,int ncol)
{
  int **out,i=0;
  out=new int*[nrow];
  for(i=0;i<nrow;i++){
    out[i]=new int[ncol];
    //out[i]--;
  }
  //out--;
  return out;
}

void free_imatrix1_r2e(int** mat,int nrow)
{
  int i=0;
  //mat++;
  for(i=0;i<nrow;i++){
    //mat[i]++;
    delete[] mat[i];
  }
  delete[] mat;
}

double** dmatrix1_r2e(int nrow,int ncol)
{
  double **out;
  int i=0;
  out=new double*[nrow];
  for(i=0;i<nrow;i++){
    out[i]=new double[ncol];
    //out[i]--;
  }
  //out--;
  return out;
}

void free_dmatrix1_r2e(double** mat,int nrow)
{
  int i=0;
  //mat++;
  for(i=0;i<nrow;i++){
    //mat[i]++;
    delete[] mat[i];
  }
  delete[] mat;
}

//
// Proximility Matrix
//

void correlation_prox(double ** data, int nrow, int ncol, double **corr_row)
{
    //TODO: Add your source code here
  int i,j,k;
  double **work,**corr,*sv,mean,temp;
  work=dmatrix1_r2e(nrow,ncol);
  corr=dmatrix1_r2e(nrow,nrow);
  sv=dvector1_r2e(nrow);
  for(i=0;i<nrow;i++){
    mean=0.0;
    for(j=0;j<ncol;j++){
      mean+=data[i][j];
    }
    mean/=double(ncol);
    sv[i]=0.0;
    for(j=0;j<ncol;j++){
      work[i][j]=data[i][j]-mean;
      sv[i]+=work[i][j]*work[i][j];
    }
    sv[i]=sqrt(sv[i]);
  }
  for(i=0;i<nrow-1;i++){
    corr[i][i]=1.0;
    for(j=i;j<nrow;j++){
      if((sv[i] == 0.0) && (sv[j] == 0.0)){
        if(data[i][1] == data[j][1]){
          temp=1.0;
        }else{
          temp=0.0;
        }
      }else if((sv[i] == 0.0) || (sv[j] == 0.0)){
        temp=0.0;
      }else{
        temp=0.0;
        for(k=0;k<ncol;k++){
          temp+=work[i][k]*work[j][k];
        }
        temp/=(sv[i]*sv[j]);
      }
      corr[i][j]=corr[j][i]=temp;
    }
  }
  corr[nrow-1][nrow-1]=1.0;

  for(i=0;i<nrow;i++){
    for(j=0;j<ncol;j++){
      corr_row[i][j]=corr[i][j];
    }
  }
  free_dmatrix1_r2e(work,nrow);
  free_dvector1_r2e(sv);
  free_dmatrix1_r2e(corr,nrow);
}

int* ellipse_sort_d(double **alldata,int nrow,int ncol,double **corr_row,int recursive_n)
{
    //TODO: Add your source code here
    int i=0,j=0;
    double ellipse_temp=0;

    int min_i=0,min_j=0;
    double minrij=0;
    double r1=0;
    double *leng,*ang;
    double *new_x,*new_y;

    int new_x_max=0,new_x_min=0;
    double dc_max=0;
    double ellipse_a=0,ellipse_b=0,ellipse_c=0;
    int cy=0;
    //double zoomratio;

    leng = new double[nrow];
    ang = new double[nrow];
    new_x = new double[nrow];
    new_y = new double[nrow];

    if(recursive_n==0)
    {
        //correlation_row(alldata,nrow,ncol,corr_row);
        correlation_prox(alldata,nrow,nrow,corr_row);
    }
    else
    {
        correlation_prox(corr_row,nrow,nrow,corr_row);
    }

    //correlation_prox(corr_row,nrow,nrow,corr_row);
    recursive_n++;
    //find_correlation_min(nrow,nrow,&minrij,&min_i,&min_j,corr_row);
    //§ä¥Xcorrelation·í¤¤µ´¹ï­È³Ì¤p­Èªºi,j
    minrij=1;
    for(i=0;i<nrow;i++)
    {
        for(j=0;j<ncol;j++)
        {
            if(fabs(corr_row[i][j])<minrij&&corr_row[i][j]!=0)
            {
                minrij=fabs(corr_row[i][j]);
                min_i=i;
                min_j=j;
            }
        }
    }

    //Rij¦pªG¬O¥¿ªº,®Ú¾Ú(0,0)´î45«×, ¦pªG¬O­tªº«h¥[45«×
    if(corr_row[min_i][min_j]>0)
    {
        r1=-M_PI_4;
    }
    else
    {
        r1=M_PI_4;
    }

    //compute_leng_ang(corr_row, min_i, min_j, nrow, ang, leng);
    //­pºâ©Ò¦³ÂI¨ì(0,0)ªºªø«×»P¨¤«×
    for(i=0;i<nrow;i++)
    {
        leng[i]=sqrt(corr_row[min_i][i]*corr_row[min_i][i]+corr_row[min_j][i]*corr_row[min_j][i]);

        if (corr_row[min_i][i]!=0&&0>corr_row[min_j][i]&&corr_row[min_i][i]>0)    //²Ä¤@¶H­­
            ang[i]=-1*atan((corr_row[min_j][i]-0)/(corr_row[min_i][i]-0))*180/M_PI;
        else if (corr_row[min_i][i]!=0&&0>corr_row[min_j][i]&&0>corr_row[min_i][i]) //²Ä¤G¶H­­
            ang[i]=180-atan((corr_row[min_j][i]-0)/(corr_row[min_i][i]-0))*180/M_PI;
        else if (corr_row[min_i][i]!=0&&corr_row[min_j][i]>0&&0>corr_row[min_i][i]) //²Ä¤T¶H­­
            ang[i]=180-atan((corr_row[min_j][i]-0)/(corr_row[min_i][i]-0))*180/M_PI;
        else if (corr_row[min_i][i]!=0&&corr_row[min_j][i]>0&&corr_row[min_i][i]>0) //²Ä¥|¶H­­
            ang[i]=360-atan((corr_row[min_j][i]-0)/(corr_row[min_i][i]-0))*180/M_PI;
        else if (corr_row[min_i][i]==0&&0>corr_row[min_j][i])        //90«×
            ang[i]=90;
        else if (corr_row[min_i][i]==0&&corr_row[min_j][i]>0)        //270«×
            ang[i]=270;

        ang[i]=ang[i]*M_PI/180;
    }

    //compute_newxy(new_x, new_y, ang, leng, r1);
    //­pºâ·sªº®y¼Ð
    for(i=0;i<nrow;i++)
    {
        new_x[i]=leng[i]*cos(ang[i]+r1);
        new_y[i]=leng[i]*sin(ang[i]+r1);
    }


    //find_newx_maxmin(new_x, new_y, &new_x_max, &new_x_min);
    //§äx³Ì¤jªºÂI®y¼Ð©Mx³Ì¤pªºÂI®y¼Ð
    double tempmax,tempmin;
    tempmax=fabs(new_x[0]);
    new_x_max=0;
    tempmin=fabs(new_y[0]);
    new_x_min=0;
    for(i=0;i<nrow;i++)
    {
        if(fabs(new_x[i])>tempmax && new_x[i]!=0)
        {
            tempmax=fabs(new_x[i]);
            new_x_max=i;
        }

        if(fabs(new_y[i])>tempmin)
        {
            tempmin=fabs(new_y[i]);
            new_x_min=i;
        }
    }


    //®Ú¾Ú¾ò¶ê¤½¦¡±a¤J¤W­±¨â­Ó®y¼Ð¡A¨D¥Xa»Pb
    ellipse_a=(new_y[new_x_min]*new_y[new_x_min]-new_y[new_x_max]*new_y[new_x_max])/(new_x[new_x_max]*new_x[new_x_max]*new_y[new_x_min]*new_y[new_x_min]-new_x[new_x_min]*new_x[new_x_min]*new_y[new_x_max]*new_y[new_x_max]);
    ellipse_b=(new_x[new_x_min]*new_x[new_x_min]-new_x[new_x_max]*new_x[new_x_max])/(new_x[new_x_min]*new_x[new_x_min]*new_y[new_x_max]*new_y[new_x_max]-new_x[new_x_max]*new_x[new_x_max]*new_y[new_x_min]*new_y[new_x_min]);

    if(ellipse_a<ellipse_b)
    {
        cy=0;
        ellipse_a=sqrt(1/ellipse_a);
        ellipse_b=sqrt(1/ellipse_b);
    }
    else
    {
        cy=1;
        ellipse_temp=ellipse_b;
        ellipse_b=sqrt(1/ellipse_a);
        ellipse_a=sqrt(1/ellipse_temp);
    }

    ellipse_c=sqrt(fabs(ellipse_a*ellipse_a-ellipse_b*ellipse_b));

    //Memo1->Lines->Add("c: "+(AnsiString)ellipse_c);
    //find_Dc_max(&dc_max, new_x, new_y, ellipse_a, ellipse_b, ellipse_c, cy);
    //­pºâ¨CÂI¨ì¨âµJ¶Zªº¶ZÂ÷©M(dcij),¨Ã§ä¥Xdcij-2*aªº³Ì¤j­È
    double temp_max,dcij;
    //int tempmaxi;
    dc_max=0;

    if(cy==0)
    {
        for(i=0;i<nrow;i++)
        {
            dcij=sqrt((new_x[i]-ellipse_c)*(new_x[i]-ellipse_c)+(new_y[i]*new_y[i]))+sqrt((new_x[i]+ellipse_c)*(new_x[i]+ellipse_c)+(new_y[i]*new_y[i]));
            temp_max=fabs(dcij-2*ellipse_a);

            //temp_max=dcij-2*ellipse_a;
            if(temp_max>dc_max)
            {
                dc_max=temp_max;
                //tempmaxi=i;
            }
        }
    }
    else
    {
        for(i=0;i<nrow;i++)
        {
            dcij=sqrt((new_y[i]-ellipse_c)*(new_y[i]-ellipse_c)+(new_x[i]*new_x[i]))+sqrt((new_y[i]+ellipse_c)*(new_y[i]+ellipse_c)+(new_x[i]*new_x[i]));
            temp_max=fabs(dcij-2*ellipse_a);
            if(temp_max>dc_max)
            {
                dc_max=temp_max;
                //tempmaxi=i;
            }
        }
    }

    if(dc_max>0.001)
    {
        if(ang)
            delete[] ang;
        if(leng)
            delete[] leng;
        if(new_x)
            delete[] new_x;
        if(new_y)
            delete[] new_y;
        return ellipse_sort_d(alldata,nrow,ncol,corr_row,recursive_n);
    }
    else
    {
        int *sortang;
        int *sortang1;

        double tempang=0;
        double maxdist=0;
        double *ang1;
        int startang=0,tempi=0;
        sortang = new int[nrow];
        sortang1 = new int[nrow];
        ang1 = new double[nrow];
        for(i=0;i<nrow;i++)
        {
            sortang[i]=i;
            ang1[i]=ang[i];
        }
        for(i=0;i<nrow;i++)
        {
            for(j=0;j<nrow;j++)
            {
                if(ang1[i]<ang1[j])
                {
                    tempang=ang1[j];
                    ang1[j]=ang1[i];
                    ang1[i]=tempang;
                    tempi=sortang[i];
                    sortang[i]=sortang[j];
                    sortang[j]=tempi;
                }
            }
        }

        for(i=0;i<nrow-1;i++)
        {
            tempang=ang1[i+1]-ang1[i];
            if(tempang>maxdist)
            {
                maxdist=tempang;
                startang=i+1;
            }
        }

        j=0;
        for(i=startang;i<nrow;i++)
        {
            sortang1[j]=sortang[i];
            j++;
        }
        for(i=0;i<startang;i++)
        {
            sortang1[j]=sortang[i];
            j++;
        }

        if(ang)
            delete[] ang;
        if(leng)
            delete[] leng;
        if(new_x)
            delete[] new_x;
        if(new_y)
            delete[] new_y;
        delete[] sortang;

        delete[] ang1;

        return sortang1;
    }
}

void ellipse_sort(double *alldataPtr, int *output_ptr, int nrow, int ncol, int recursive_n)
{
    double **alldata = fun1dto2dArray_r2e(alldataPtr,nrow,ncol);
    double **corr_row;
    corr_row=dmatrix1_r2e(nrow,nrow);

    int* output_array = ellipse_sort_d(alldata, nrow, ncol, corr_row, recursive_n);

    for(int i=0; i<nrow; i++)
    {
        output_ptr[i] =  output_array[i];   
    }

    //free_dmatrix1_r2e(alldata, nrow);
    for(int i = 0; i < nrow; i++){
        delete[] alldata[i];
    }
    delete[] alldata;
    free_dmatrix1_r2e(corr_row, nrow);
    //free_ivector1_r2e(output_array);
    delete[] output_array;

}

#ifdef __cplusplus
}
#endif

// [[Rcpp::export]]
Rcpp::IntegerVector ellipse_sort_R(Rcpp::NumericMatrix data) {
    int nrow = data.nrow();
    int ncol = data.ncol();
    
    // Output array used to store sorted resultsv
    int* output = new int[nrow];

    // Calling C++ function
    ellipse_sort(&data[0], output, nrow, ncol, 0); // &data[0] passes the matrix in as a pointer
    
    // The result of the conversion is a vector of R
    Rcpp::IntegerVector result(nrow);
    for (int i = 0; i < nrow; i++) {
        result[i] = output[i] + 1; // R's indexes start at 1
    }
    
    // Clean up the memory
    delete[] output;
    
    return result;
}
