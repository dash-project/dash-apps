#ifndef _DGRAPH
#define _DGRAPH

#define BLOCK_SIZE  128
#define SMALL_BLOCK_SIZE 32

#include "npbparams.h"

typedef struct DGArc_s DGArc;
typedef struct DGNode_s DGNode;
typedef struct DGNode_new_s DGNodeNew;
typedef struct DGraph_s DGraph;

#define fielddim 4
#define FEAT_MAX_LEN ((NUM_SAMPLES) * (fielddim) * (2))
typedef struct DGNode_Feat_s {
  int len;
  double val[FEAT_MAX_LEN];
} DGNode_Feat;

struct DGArc_s {
  int id;
  DGNode *tail,*head;
  int length,width,attribute,maxWidth;
};

struct DGNode_s {
  DGArc **inArc, **outArc;
  int maxInDegree,maxOutDegree;
  int inDegree,outDegree;
  int id;
  char name[BLOCK_SIZE];
  int in[4], out[4];
  int depth,height,width;
  int color,attribute,address,verified;
  DGNode_Feat feat;
};

struct DGNode_new_s {
  int maxInDegree,maxOutDegree;
  int inDegree,outDegree;
  int id;
  char name[BLOCK_SIZE];
  int in[4], out[4];
  int depth,height,width;
  int color,attribute,address,verified;
  DGNode_Feat feat;
};

struct DGraph_s{
  int maxNodes,maxArcs;
  int id;
  char *name;
  int numNodes,numArcs;
  DGNode ** node;
  DGArc ** arc;
};

DGArc *newArc(DGNode *tl,DGNode *hd);
void arcShow(DGArc *ar);
DGNode *newNode(const char nm[BLOCK_SIZE]);
void nodeShow(DGNode* nd);

DGraph* newDGraph(char const *nm);
int AttachNode(DGraph *dg,DGNode *nd);
int AttachArc(DGraph *dg,DGArc* nar);
void graphShow(DGraph *dg,int DetailsLevel);

#endif
