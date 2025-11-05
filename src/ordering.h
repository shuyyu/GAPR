#include "list.h"

#ifdef __cplusplus
extern "C"
{
#endif

List* getOrderList(int **merge,int root);
double select_average_list(double **data,List *rlist,List *clist);
double select_average(double **data,List *list);

void ordering_GrandPa_d(double **data,int n,int **merge);
void ordering_GrandPa_s(double **data,int n,int **merge);

void ordering_Uncle_d(double **data,int n,int **merge);
void ordering_Uncle_s(double **data,int n,int **merge);

void ordering_external(int *eo,int **merge,int n);

#ifdef __cplusplus
}
#endif 
