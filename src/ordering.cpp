//
// Ordering 
//

#include "ordering.h"


//
// declare
//

int isTerminal(int **merge,int pt);

int* ivector1_3(int length)
{
  int *out;
  out=new int[length+1];
  //out--;
  return out;
}

void free_ivector1_3(int* vec)
{
  //vec++;
  delete[] vec;
}

//
// utility 
//

List* getOrderList(int **merge,int root)
{
  int p,flag;
  gInt *obj;
  Stack s;
  List *lt;
  lt=new List();
  p=root;
  flag=1;
  while(flag){
    while(p>0){
      s.push(new gInt(p));
      p=merge[p][1];
    }
    lt->append(s.pop());
    if(!s.empty()){
      obj=(gInt*) s.pop();
      p=obj->data;
      delete obj;
      if(isTerminal(merge,merge[p][2])){
        s.push(new gInt(merge[p][2]));
        p=0;
      }else{
        p=merge[p][2];
      }
    }
    if((p == 0) && (s.empty())) flag=0;
  }
  return lt;
}

int isTerminal(int **merge,int pt)
{
  if((merge[pt][1] == 0) && (merge[pt][2] == 0)){
    return 1;
  }else{
    return 0;
  }
}

double select_average_list(double **data,List *rlist,List *clist)
{
  double out;
  int i,j,pi,pj,nr,nc;
  rlist->first();
  clist->first();
  nr=rlist->length();
  nc=clist->length();
  out=0.0;
  for(i=1;i<=nr;i++){
    pi=((gInt*) rlist->GetAt(i))->data;
    for(j=1;j<=nc;j++){
      pj=((gInt*) clist->GetAt(j))->data;
      out+=data[pi][pj];
    }
  }
  out/=double(nr*nc);
  return out;
}

double select_average(double **data,List *list)
{
  double out;
  int i,j,n,*row;
  n=list->length();
  row=list2vector(list);
  out=0.0;
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++){
      out+=data[row[i]][row[j]];
    }
  }
  out/=double(n*n);
  free_ivector1_3(row);
  return out;
}

//
// GrandPa - distance
//

void ordering_GrandPa_d_left(double **data,int **merge,int root,List *rlt);
void ordering_GrandPa_d_right(double **data,int **merge,int root,List *rlt);

void ordering_GrandPa_d(double **data,int n,int **merge)
{
  List *left,*right,*tlist;
  double pleft,pright;
  int root,temp;
  root=2*n-1;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average(data,left);
  pright=select_average(data,right);
  if(pleft<pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
    tlist=left;
    left=right;
    right=tlist;
  }
  ordering_GrandPa_d_left(data,merge,merge[root][1],right);
  ordering_GrandPa_d_right(data,merge,merge[root][2],left);
  delete left;
  delete right;
}

void ordering_GrandPa_d_left(double **data,int **merge,int root,List *rlt)
{
  int temp;
  double pleft,pright;
  List *left,*right;
  if(isTerminal(merge,root)) return;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average_list(data,left,rlt);
  pright=select_average_list(data,right,rlt);
  if(pleft<pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
  }
  delete left;
  delete right;
  ordering_GrandPa_d_left(data,merge,merge[root][1],rlt);
  ordering_GrandPa_d_left(data,merge,merge[root][2],rlt);
}

void ordering_GrandPa_d_right(double **data,int **merge,int root,List *rlt)
{
  int temp;
  double pleft,pright;
  List *left,*right;
  if(isTerminal(merge,root)) return;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average_list(data,left,rlt);
  pright=select_average_list(data,right,rlt);
  if(pleft>pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
  }
  delete left;
  delete right;
  ordering_GrandPa_d_right(data,merge,merge[root][1],rlt);
  ordering_GrandPa_d_right(data,merge,merge[root][2],rlt);
}

//
// GrandPa - simility
//

void ordering_GrandPa_s_left(double **data,int **merge,int root,List *rlt);
void ordering_GrandPa_s_right(double **data,int **merge,int root,List *rlt);

void ordering_GrandPa_s(double **data,int n,int **merge)
{
  List *left,*right,*tlist;
  double pleft,pright;
  int root,temp;
  root=2*n-1;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average(data,left);
  pright=select_average(data,right);
  if(pleft>pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
    tlist=left;
    left=right;
    right=tlist;
  }
  ordering_GrandPa_s_left(data,merge,merge[root][1],right);
  ordering_GrandPa_s_right(data,merge,merge[root][2],left);
  delete left;
  delete right;
}

void ordering_GrandPa_s_left(double **data,int **merge,int root,List *rlt)
{
  int temp;
  double pleft,pright;
  List *left,*right;
  if(isTerminal(merge,root)) return;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average_list(data,left,rlt);
  pright=select_average_list(data,right,rlt);
  if(pleft>pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
  }
  delete left;
  delete right;
  ordering_GrandPa_s_left(data,merge,merge[root][1],rlt);
  ordering_GrandPa_s_left(data,merge,merge[root][2],rlt);
}

void ordering_GrandPa_s_right(double **data,int **merge,int root,List *rlt)
{
  int temp;
  double pleft,pright;
  List *left,*right;
  if(isTerminal(merge,root)) return;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average_list(data,left,rlt);
  pright=select_average_list(data,right,rlt);
  if(pleft<pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
  }
  delete left;
  delete right;
  ordering_GrandPa_s_right(data,merge,merge[root][1],rlt);
  ordering_GrandPa_s_right(data,merge,merge[root][2],rlt);
}

//
// Uncle - distance
//

void ordering_Uncle_d_left(double **data,int **merge,int root,List *rlt);
void ordering_Uncle_d_right(double **data,int **merge,int root,List *rlt);

void ordering_Uncle_d(double **data,int n,int **merge)
{
  List *left,*right,*tlist;
  double pleft,pright;
  int root,temp;
  root=2*n-1;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average(data,left);
  pright=select_average(data,right);
  if(pleft<pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
    tlist=left;
    left=right;
    right=tlist;
  }
  ordering_Uncle_d_left(data,merge,merge[root][1],right);
  ordering_Uncle_d_right(data,merge,merge[root][2],left);
}

void ordering_Uncle_d_left(double **data,int **merge,int root,List *rlt)
{
  int temp;
  double pleft,pright;
  List *left,*right,*tlist;
  if(isTerminal(merge,root)) {
    delete rlt;
    return;
  }
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average_list(data,left,rlt);
  pright=select_average_list(data,right,rlt);
  delete rlt;
  if(pleft<pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
    tlist=left;
    left=right;
    right=tlist;
  }
  ordering_Uncle_d_left(data,merge,merge[root][1],right);
  ordering_Uncle_d_right(data,merge,merge[root][2],left);
}

void ordering_Uncle_d_right(double **data,int **merge,int root,List *rlt)
{
  int temp;
  double pleft,pright;
  List *left,*right,*tlist;
  if(isTerminal(merge,root)) {
    delete rlt;
    return;
  }
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average_list(data,left,rlt);
  pright=select_average_list(data,right,rlt);
  delete rlt;
  if(pleft>pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
    tlist=left;
    left=right;
    right=tlist;
  }
  ordering_Uncle_d_left(data,merge,merge[root][1],right);
  ordering_Uncle_d_right(data,merge,merge[root][2],left);
}

//
// Uncle - simility
//

void ordering_Uncle_s_left(double **data,int **merge,int root,List *rlt);
void ordering_Uncle_s_right(double **data,int **merge,int root,List *rlt);

void ordering_Uncle_s(double **data,int n,int **merge)
{
  List *left,*right,*tlist;
  double pleft,pright;
  int root,temp;
  root=2*n-1;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average(data,left);
  pright=select_average(data,right);
  if(pleft>pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
    tlist=left;
    left=right;
    right=tlist;
  }
  ordering_Uncle_s_left(data,merge,merge[root][1],right);
  ordering_Uncle_s_right(data,merge,merge[root][2],left);
}

void ordering_Uncle_s_left(double **data,int **merge,int root,List *rlt)
{
  int temp;
  double pleft,pright;
  List *left,*right,*tlist;
  if(isTerminal(merge,root)) {
    delete rlt;
    return;
  }
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average_list(data,left,rlt);
  pright=select_average_list(data,right,rlt);
  delete rlt;
  if(pleft>pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
    tlist=left;
    left=right;
    right=tlist;
  }
  ordering_Uncle_s_left(data,merge,merge[root][1],right);
  ordering_Uncle_s_right(data,merge,merge[root][2],left);
}

void ordering_Uncle_s_right(double **data,int **merge,int root,List *rlt)
{
  int temp;
  double pleft,pright;
  List *left,*right,*tlist;
  if(isTerminal(merge,root)) {
    delete rlt;
    return;
  }
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  pleft=select_average_list(data,left,rlt);
  pright=select_average_list(data,right,rlt);
  delete rlt;
  if(pleft<pright){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
    tlist=left;
    left=right;
    right=tlist;
  }
  ordering_Uncle_s_left(data,merge,merge[root][1],right);
  ordering_Uncle_s_right(data,merge,merge[root][2],left);
}

//
// External
//

void external(int *ieo,int **merge,int root);

void ordering_external(int *eo,int **merge,int n)
{
  int *ieo,i,root;
  ieo=ivector1_3(n);
  for(i=1;i<=n;i++){
    ieo[eo[i]]=i;
  }
  root=2*n-1;
  external(ieo,merge,root);
  free_ivector1_3(ieo);
}

void external(int *ieo,int **merge,int root)
{
  int c1,c2,i,j,pi,pj,nl,nr,temp;
  List *left,*right;
  if(isTerminal(merge,root)) return;
  left=getOrderList(merge,merge[root][1]);
  right=getOrderList(merge,merge[root][2]);
  nl=left->length();
  nr=right->length();
  c1=0;
  c2=0;
  left->first();
  for(i=1;i<=nl;i++){
    pi=((gInt*) left->val())->data;
    right->first();
    for(j=1;j<=nr;j++){
      pj=((gInt*) right->val())->data;
      if(ieo[pi] > ieo[pj]) c1++;
      else c2++;
      right->next();
    }
    left->next();
  }
  if(c1 > c2){
    temp=merge[root][1];
    merge[root][1]=merge[root][2];
    merge[root][2]=temp;
  }
  delete left;
  delete right;
  external(ieo,merge,merge[root][1]);
  external(ieo,merge,merge[root][2]);
}
