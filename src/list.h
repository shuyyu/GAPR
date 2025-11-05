#ifdef __cplusplus
extern "C"
{
#endif

#ifndef NULL
#define NULL 0
#endif

class gObject
{
  public:
  virtual ~gObject();
};

class gInt : public gObject
{
public:
  int data;
  gInt(int _data);
  virtual ~gInt();
};

class ListNode
{
protected:
  ListNode *_next;
  ListNode *_prev;
public:
  gObject *_val;
  ListNode(gObject *val);
  ListNode();
  virtual ~ListNode();
  ListNode *next();
  ListNode *prev();
  ListNode *insert(ListNode *b);
  ListNode *remove();
  void splice(ListNode *b);
friend class List;
};

class List : public gObject
{
private:
  ListNode *header;
  ListNode *win;
  int _length,_index;
public:
  List();
  ~List();
  gObject* insert(gObject* val);
  gObject* append(gObject* val);
  List* append(List* l);
  gObject* prepend(gObject* val);
  gObject* remove();
  gObject* val(gObject* val);
  gObject* val();
  gObject* GetAt(int _pt);
  gObject* next();
  gObject* prev();
  gObject* first();
  gObject* last();
  int length();
  int isFirst();
  int isLast();
  int isHead();
};

class Stack : public gObject
{
private:
  List s;
public:
  Stack();
  virtual ~Stack();
  void push(gObject *v);
  gObject* pop();
  int empty();
  int size();
  gObject* top();
  gObject* nextToTop();
  gObject* bottom();
};

//
// utility function
//

int* list2vector(List *list);


#ifdef __cplusplus
}
#endif 
