#ifndef RF_TERM_OPS_H
#define RF_TERM_OPS_H
#include "terminal.h"
Terminal *makeTerminal(void);
void freeTerminal(Terminal *parent);
void stackTermLMIIndex(Terminal *tTerm, unsigned int size);
void unstackTermLMIIndex(Terminal *tTerm);
void freeTerminalNodeLocalSurvivalStructures(Terminal *tTerm);
void freeTerminalNodeSurvivalStructuresIntermediate(Terminal *tTerm);
void freeTerminalNodeSurvivalStructuresFinal(Terminal *tTerm);
void freeTerminalNodeNonSurvivalStructures(Terminal *tTerm);
void stackAtRiskAndEventCount(Terminal *tTerm, unsigned int eTypeSize, unsigned int mTimeSize);
void unstackAtRiskAndEventCount(Terminal *tTerm);
void stackEventTimeIndex(Terminal *tTerm, unsigned int mTimeSize);
void unstackEventTimeIndex(Terminal *tTerm);
void stackLocalRatio(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalRatio(Terminal *tTerm);
void stackLocalRatioHazard(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalRatioHazard(Terminal *tTerm);
void stackLocalSurvival(Terminal *tTerm, unsigned int eTimeSize);
void unstackLocalSurvival(Terminal *tTerm);
void stackLocalNelsonAalen(Terminal *tTerm, unsigned int eTimeSize);
void unstackLocalNelsonAalen(Terminal *tTerm);
void stackLocalCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalCSH(Terminal *tTerm);
void stackLocalCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalCIF(Terminal *tTerm);
void stackNelsonAalen(Terminal *tTerm, unsigned int sTimeSize);
void unstackNelsonAalen(Terminal *tTerm);
void stackSurvival(Terminal *tTerm, unsigned int sTimeSize);
void unstackSurvival(Terminal *tTerm);
void stackCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize);
void unstackCSH(Terminal *tTerm);
void stackCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize);
void unstackCIF(Terminal *tTerm);
void stackMortality(Terminal *tTerm, unsigned int eTypeSize);
void unstackMortality(Terminal *tTerm);
void stackMultiClassProb(Terminal *tTerm, unsigned int rfCount, unsigned int *rfSize);
void stackMultiClassProbPartial(Terminal *tTerm, unsigned int rfCount);
void unstackMultiClassProb(Terminal *tTerm);
void stackMeanResponse(Terminal *tTerm, unsigned int rnfCount);
void unstackMeanResponse(Terminal *tTerm);
void stackMemberStream(Terminal *tTerm, unsigned int membrSize);
void unstackMemberStream(Terminal *tTerm);
#endif
