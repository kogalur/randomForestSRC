#ifndef RF_STACK_OUTPUT_QQ_H
#define RF_STACK_OUTPUT_QQ_H
void stackTNQualitativeObjectsKnown(char     mode,
                                    uint   **pRF_RMBR_ID_,
                                    uint   **pRF_AMBR_ID_,
                                    uint   **pRF_TN_RCNT_,
                                    uint   **pRF_TN_ACNT_,
                                    uint   **pRF_OOB_SZ_,
                                    uint   **pRF_IBG_SZ_);
void stackTNQualitativeObjectsUnknown(char     mode,
                                      uint   **pRF_TN_RCNT_,
                                      uint   **pRF_TN_ACNT_,
                                      uint   **pRF_TN_OCNT_,
                                      uint   **pRF_TN_ICNT_);
void stackTNQuantitativeForestObjectsPtrOnly(char mode);
void unstackTNQuantitativeForestObjectsPtrOnly(char mode);
void stackTNQuantitativeTreeObjectsPtrOnly(uint treeID);
void unstackTNQuantitativeTreeObjectsPtrOnly(uint treeID);
void saveTNQuantitativeTreeObjects(uint treeID);
void stackTNQuantitativeForestObjectsOutput(char mode);
void writeTNQuantitativeForestObjectsOutput(char mode);
void stackTNQualitativeObjectsUnknownMembership(char   mode, uint **pRF_OMBR_ID_, uint **pRF_IMBR_ID_);
#endif
