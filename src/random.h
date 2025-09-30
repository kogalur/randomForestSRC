#ifndef RF_RANDOM_H
#define RF_RANDOM_H
void randomStack(uint bSize, uint bnpSize);
void randomUnstack(uint bSize, uint bnpSize);
void randomSetChainParallel(uint b, int value);
void randomSetUChainParallel(uint b, int value);
void randomSetUChainParallelVimp(uint b, int value);
void randomSetChainParallelVimp(uint p, int value);
void randomSetChainSerial(uint b, int value);
void randomSetUChainSerial(uint b, int value);
void randomSetUChainSerialVimp(uint b, int value);
void randomSetChainSerialVimp(uint p, int value);
int randomGetChainParallel(uint b);
int randomGetUChainParallel(uint b);
int randomGetUChainParallelVimp(uint b);
int randomGetChainParallelVimp(uint p);
int randomGetChainSerial(uint b);
int randomGetUChainSerial(uint b);
int randomGetUChainSerialVimp(uint b);
int randomGetChainSerialVimp(uint p);
float randomChainParallel(uint b);
float randomUChainParallel(uint b);
float randomUChainParallelVimp(uint b);
float randomChainParallelVimp(uint p);
float randomChainSerial(uint b);
float randomUChainSerial(uint b);
float randomUChainSerialVimp(uint b);
float randomChainSerialVimp(uint p);
float ran1_generic(int *iy, int *iv, int *idum);
void lcgenerator(unsigned int *seed, unsigned char reset);
float ran1_original(int *idum);
#endif
