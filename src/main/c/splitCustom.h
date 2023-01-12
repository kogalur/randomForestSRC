/* 
   vvvvvvvv External Constants Below -- Do Not Change vvvvvvvv
*/

#define LEFT      0x01
#define RIGHT     0x02

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

void registerCustomFunctions(void);

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
                                          unsigned int    maxLevel,

                                          double        **feature,
                                          unsigned int    featureCount),

                          
                          unsigned int family,
                          unsigned int slot);


/* 
   ^^^^^^^^ Do Not Touch These Delarations Above ^^^^^^^^
*/




/*
   Declare your custom functions below:
*/

double getCustomSplitStatisticMultivariateRegression  (unsigned int  n,
                                                       char         *membership,
                                                       double       *time,
                                                       double       *event,
                                                       
                                                       unsigned int  eventTypeSize,
                                                       unsigned int  eventTimeSize,
                                                       double       *eventTime,
                                                       
                                                       double       *response,
                                                       double        mean,
                                                       double        variance,
                                                       unsigned int  maxLevel,
                                                       
                                                       double      **feature,
                                                       unsigned int  featureCount);


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
                                                          unsigned int  maxLevel,
                                                      
                                                          double      **feature,
                                                          unsigned int  featureCount);


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
                                        unsigned int  maxLevel,

                                        double      **feature,
                                        unsigned int  featureCount);


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
                                             unsigned int  maxLevel,

                                             double      **feature,
                                             unsigned int  featureCount);


unsigned int *alloc_uivector(unsigned int nh);
void          dealloc_uivector(unsigned int *v, unsigned int nh);

double       *alloc_dvector(double *v, unsigned int nh);
void          dealloc_dvector(double *v, unsigned int nh);

unsigned int **alloc_uimatrix(unsigned int n2h, unsigned int nh);
void          dealloc_uimatrix(unsigned int **v, unsigned int n2h, unsigned int nh);
