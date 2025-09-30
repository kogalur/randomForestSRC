#ifndef RF_TRACE_H
#define RF_TRACE_H
#define SUMM_DEF_TRACE  0x00000001
#define SUMM_LOW_TRACE  0x00000002
#define SUMM_MED_TRACE  0x00000004
#define SUMM_HGH_TRACE  0x00000008
#define SPLT_DEF_TRACE  0x00000010
#define SPLT_LOW_TRACE  0x00000020
#define SPLT_MED_TRACE  0x00000040
#define SPLT_HGH_TRACE  0x00000080
#define FORK_DEF_TRACE  0x00000100
#define MISS_LOW_TRACE  0x00000200
#define MISS_MED_TRACE  0x00000400
#define MISS_HGH_TRACE  0x00000800
#define OUTP_DEF_TRACE  0x00001000
#define NUMR_DEF_TRACE  0x00002000
#define FACT_LOW_TRACE  0x00004000
#define FACT_HGH_TRACE  0x00008000
#define ENSB_LOW_TRACE  0x00010000
#define ENSB_HGH_TRACE  0x00020000
#define BOOT_MED_TRACE  0x00040000
#define VIMP_LOW_TRACE  0x00080000
#define NODE_DEF_TRACE  0x00100000
#define TIME_DEF_TRACE  0x00200000
#define RAND_DEF_TRACE  0x00400000
#define QUAN_DEF_TRACE  0x00800000
#define SAMP_DEF_TRACE  0x01000000
#define SAMP_LOW_TRACE  0x02000000
#define TURN_OFF_TRACE  0x00000000
void setTraceFlag(unsigned int traceFlag, unsigned int tree);
unsigned int getTraceFlag(unsigned int tree);
unsigned int updateTimeStamp(unsigned int before);
void   setMaxMemoryAllocation(size_t value);
void   setMinMemoryAllocation(size_t value);
void   setProbeMemoryAllocation(size_t value);
size_t getMaxMemoryAllocation();
size_t getMinMemoryAllocation();
size_t getProbeMemoryAllocation();
void   increaseMemoryAllocation(size_t amount);
void   increaseProbeMemoryAllocation(size_t amount);
void   decreaseMemoryAllocation(size_t amount);
void   decreaseProbeMemoryAllocation(size_t amount);
void memoryCheck(void);
#endif
