#ifndef RFSRCSPLITCUST_H
#define RFSRCSPLITCUST_H

/* 
   vvvvvvvv External Constants Below -- Do Not Change vvvvvvvv
*/

#define LEFT      0x01
#define RIGHT     0x00

#define CLAS_FAM     0
#define REGR_FAM     1
#define SURV_FAM     2
#define CRSK_FAM     3

/* 
   ^^^^^^^^ External Constants Above -- Do Not Change ^^^^^^^^
*/


/* 
   vvvvvvvv Do Not Touch These Delarations Below vvvvvvvv
*/

void registerCustomFunctions();

extern void registerThis (double (*func) (unsigned int    n,
                                          char           *membership,
                                          double         *time,
                                          double         *event,

                                          unsigned int    eventTypeSize,
                                          unsigned int    eventTimeSize,
                                          double         *eventTime,

                                          double         *response,
                                          double          mean,
                                          double          variance,
                                          unsigned int    maxLevel),

                          unsigned int family,
                          unsigned int slot);


/* 
   ^^^^^^^^ Do Not Touch These Delarations Above ^^^^^^^^
*/




/*
   Declare your custom funtions below:
*/

double getCustomSplitStatisticMultivariateRegression (unsigned int  n,
                                                      char         *membership,
                                                      double       *time,
                                                      double       *event,

                                                      unsigned int  eventTypeSize,
                                                      unsigned int  eventTimeSize,
                                                      double       *eventTime,

                                                      double       *response,
                                                      double        mean,
                                                      double        variance,
                                                      unsigned int  maxLevel);

double getCustomSplitStatisticMultivariateClassification (unsigned int  n,
                                                          char         *membership,
                                                          double       *time,
                                                          double       *event,

                                                          unsigned int  eventTypeSize,
                                                          unsigned int  eventTimeSize,
                                                          double       *eventTime,

                                                          double       *response,
                                                          double        mean,
                                                          double        variance,
                                                          unsigned int  maxLevel);

double getCustomSplitStatisticSurvival (unsigned int  n,
                                        char         *membership,
                                        double       *time,
                                        double       *event,

                                        unsigned int  eventTypeSize,
                                        unsigned int  eventTimeSize,
                                        double       *eventTime,

                                        double       *response,
                                        double        mean,
                                        double        variance,
                                        unsigned int  maxLevel);

double getCustomSplitStatisticCompetingRisk (unsigned int  n,
                                             char         *membership,
                                             double       *time,
                                             double       *event,

                                             unsigned int  eventTypeSize,
                                             unsigned int  eventTimeSize,
                                             double       *eventTime,

                                             double       *response,
                                             double        mean,
                                             double        variance,
                                             unsigned int  maxLevel);


unsigned int *alloc_uivector(unsigned int nh);
void          dealloc_uivector(unsigned int *v, unsigned int nh);

double       *alloc_dvector(double *v, unsigned int nh);
void          dealloc_dvector(double *v, unsigned int nh);

unsigned int **alloc_uimatrix(unsigned int n2h, unsigned int nh);
void          dealloc_uimatrix(unsigned int **v, unsigned int n2h, unsigned int nh);

/* RF_CRAN_BEG */
#endif
/* RF_CRAN_END */
