#ifndef RF_NODE_OPS_H
#define RF_NODE_OPS_H
Node *makeNode(unsigned int xSize);
void freeNodeGeneric(Node *parent);
void freeNodeNew(Node *parent);
void setParent(
  Node *daughter,
  Node *parent
);
void setLeftDaughter(
   Node *daughter,
   Node *parent
);
void setRightDaughter(
  Node *daughter,
  Node *parent
);
void stackMPSign(Node *node, unsigned int mpIndexSize);
void unstackMPSign(Node *node);
void stackFMPSign(Node *node, unsigned int fmpIndexSize);
void unstackFMPSign(Node *node);
void stackNodeLMPIndex(Node *node, unsigned int size);
void unstackNodeLMPIndex(Node *node);
void stackNodeFLMPIndex(Node *node, unsigned int size);
void unstackNodeFLMPIndex(Node *node);
void stackSplitDepth(Node *tNode, unsigned int depth);
void unstackSplitDepth(Node *tNode);
#endif
