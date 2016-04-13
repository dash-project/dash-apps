#ifndef _DGRAPH
#define _DGRAPH

#define BLOCK_SIZE  128
#define SMALL_BLOCK_SIZE 32

//#include "npbparams.h"

//#define fielddim 4
//#define FEAT_MAX_LEN ((NUM_SAMPLES) * (fielddim) * (2))
//
typedef struct DGNode_s DGNode;
typedef struct DGArc_s DGArc;

struct DGArc_s{
  int id;
  DGNode *tail,*head;
  int length,width,attribute,maxWidth;
};

struct DGNode_s{
  int maxInDegree,maxOutDegree;
  int inDegree,outDegree;
  int id;
  char *name;
  DGArc **inArc,**outArc;
  int depth,height,width;
  int color,attribute,address,verified;
  void *feat;
};

typedef struct{
  int maxNodes,maxArcs;
  int id;
  char *name;
  int numNodes,numArcs;
  DGNode **node;
  DGArc **arc;
} DGraph;

DGArc *newArc(DGNode *tl,DGNode *hd);
void arcShow(DGArc *ar);
DGNode *newNode(char const * nm);
void nodeShow(DGNode* nd);

DGraph* newDGraph(char const * nm);
int AttachNode(DGraph *dg,DGNode *nd);
int AttachArc(DGraph *dg,DGArc* nar);
void graphShow(DGraph *dg,int DetailsLevel);

#endif
