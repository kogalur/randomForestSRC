
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "entryGeneric.h"
void processDefaultGrow(void) {
  RF_optHigh = RF_optHigh & (~OPT_MEMB_PRUN);
  RF_ptnCount             = 0;
  RF_optHigh = RF_optHigh & (~OPT_PART_PLOT);
  RF_partialLength = 0;
  RF_opt = RF_opt & (~OPT_VIMP_JOIN);
  RF_opt = RF_opt & (~OPT_OUTC_TYPE);
  RF_opt = RF_opt & (~OPT_COMP_RISK);
  RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
  RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
  RF_frSize = RF_fobservationSize = 0;
  RF_xMarginalSize = 0;
  if (RF_opt & OPT_IMPU_ONLY) {
    RF_opt                  = RF_opt & (OPT_IMPU_ONLY | OPT_BOOT_TYP1 | OPT_BOOT_TYP2);
    RF_optHigh              = RF_optHigh & (OPT_MISS_SKIP | OPT_BOOT_SWOR);
  }
  else {
  }
  RF_opt = RF_opt | OPT_MISS_OUT;
  RF_opt = RF_opt | OPT_LEAF;
  if ( !(RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2)) {
    RF_optHigh = RF_optHigh & (~OPT_BOOT_SWOR);
    RF_opt                  = RF_opt & (~OPT_OENS);
    RF_opt                  = RF_opt & (~OPT_PERF);
    RF_opt                  = RF_opt & (~OPT_VIMP);
    RF_optHigh              = RF_optHigh & (~OPT_CSE);
    RF_optHigh              = RF_optHigh & (~OPT_CSV);
    if (RF_opt & OPT_PROX) {
      RF_opt = RF_opt | OPT_PROX_IBG;
      RF_opt = RF_opt | OPT_PROX_OOB;
    }
    if (RF_optHigh & OPT_DIST) {
      RF_optHigh = RF_optHigh | OPT_DIST_IBG;
      RF_optHigh = RF_optHigh | OPT_DIST_OOB;
    }
    if (RF_optHigh & OPT_WGHT) {
      RF_optHigh = RF_optHigh | OPT_WGHT_IBG;
      RF_optHigh = RF_optHigh | OPT_WGHT_OOB;
    }
  }
  else {
  }
  if ((RF_splitRule == USPV_NRM) || (RF_splitRule == USPV_WT_OFF) || (RF_splitRule == USPV_WT_HVY)) {
    RF_opt                  = RF_opt & (~OPT_PERF);
    RF_opt                  = RF_opt & (~OPT_VIMP);
    RF_opt                  = RF_opt & (~OPT_OENS);
    RF_opt                  = RF_opt & (~OPT_FENS);
    RF_optHigh              = RF_optHigh & (~OPT_CSE);
    RF_optHigh              = RF_optHigh & (~OPT_CSV);
    RF_ySize = 0;
  }
  if (RF_opt & OPT_PERF) {
  }
  else {
    RF_opt = RF_opt & (~OPT_VIMP);
    RF_optHigh              = RF_optHigh & (~OPT_CSE);
    RF_optHigh              = RF_optHigh & (~OPT_CSV);
  }
  if (RF_opt & OPT_TREE) {
    RF_opt = RF_opt | OPT_SEED;
  }
  else {
    RF_opt = RF_opt & (~OPT_SEED);
  }
  RF_opt = RF_opt & (~OPT_EMPR_RISK);
  if ((RF_opt & OPT_OENS) || (RF_opt & OPT_FENS)) {
  }
  else {
    RF_optHigh = RF_optHigh & (~OPT_TERM_OUTG);
  }
  if (!(RF_opt & OPT_OENS)) {
    RF_opt = RF_opt & (~OPT_PERF);
    RF_optHigh              = RF_optHigh & (~OPT_CSE);
    RF_optHigh              = RF_optHigh & (~OPT_CSV);
  }
  if (RF_vtry > 0) {
    RF_opt = RF_opt & (~OPT_VIMP);
    RF_nImpute = 1;
  }
}
void processDefaultPredict(void) {
  char mode;
  RF_opt = RF_opt & (~OPT_IMPU_ONLY);
  RF_opt = RF_opt & (~OPT_TREE);
  RF_opt = RF_opt & (~OPT_SEED);
  RF_nImpute = 1;
  RF_opt = RF_opt | OPT_LEAF;
  RF_opt  = RF_opt | OPT_MISS_OUT;
  RF_optHigh = RF_optHigh & (~OPT_MEMB_OUTG);
  RF_optHigh = RF_optHigh & (~OPT_TERM_OUTG);
  RF_vtry                 = 0;
  RF_vtryArray            = NULL;
  if(RF_fobservationSize > 0) {
    mode = RF_PRED;
  }
  else {
    mode = RF_REST;
  }
  switch (mode) {
  case RF_PRED:
    RF_opt = RF_opt & (~OPT_OUTC_TYPE);
    RF_optHigh = RF_optHigh & (~OPT_PART_PLOT);
    RF_partialLength = RF_partialLength2 = 0;
    RF_opt = RF_opt & (~OPT_OENS);
    if (RF_ySize == 0) {
      RF_opt = RF_opt & (~OPT_PERF);
      RF_opt = RF_opt & (~OPT_VIMP);
      RF_optHigh              = RF_optHigh & (~OPT_CSE);
      RF_optHigh              = RF_optHigh & (~OPT_CSV);
      RF_opt = RF_opt & (~OPT_FENS);
    }
    else {
      if (RF_frSize == 0) {
        RF_opt                  = RF_opt & (~OPT_PERF);
        RF_opt                  = RF_opt & (~OPT_VIMP);
        RF_optHigh              = RF_optHigh & (~OPT_CSE);
        RF_optHigh              = RF_optHigh & (~OPT_CSV);
      }
    }
    if (RF_opt & OPT_PROX) {
      RF_opt = RF_opt | OPT_PROX_IBG;
      RF_opt = RF_opt | OPT_PROX_OOB;
    }
    if (RF_optHigh & OPT_DIST) {
      RF_optHigh = RF_optHigh | OPT_DIST_IBG;
      RF_optHigh = RF_optHigh | OPT_DIST_OOB;
    }
    if (RF_optHigh & OPT_WGHT) {
      RF_optHigh = RF_optHigh | OPT_WGHT_IBG;
      RF_optHigh = RF_optHigh | OPT_WGHT_OOB;
    }
    if (!(RF_opt & OPT_FENS)) {
      RF_opt                  = RF_opt & (~OPT_PERF);
      RF_optHigh              = RF_optHigh & (~OPT_CSE);
      RF_optHigh              = RF_optHigh & (~OPT_CSV);
    }
    RF_opt = RF_opt & (~OPT_SPLDPTH_1) & (~OPT_SPLDPTH_2);
    RF_opt = RF_opt & (~OPT_EMPR_RISK);
    if (RF_opt & OPT_FENS) {
    }
    else {
    }
    break;
  case RF_REST:
    RF_frSize = RF_fobservationSize = 0;
    if (RF_opt & OPT_OUTC_TYPE) {
      RF_optHigh = RF_optHigh & (~OPT_PART_PLOT);
      RF_partialLength = RF_partialLength2 = 0;
      RF_optHigh = RF_optHigh & (~OPT_MEMB_INCG);
      RF_optHigh = RF_optHigh & (~OPT_TERM_INCG);
    }
    else if (RF_optHigh & OPT_PART_PLOT) {
      RF_opt = RF_opt & (~OPT_OUTC_TYPE);
      RF_opt = RF_opt & (~OPT_PERF);
      RF_optHigh              = RF_optHigh & (~OPT_CSE);
      RF_optHigh              = RF_optHigh & (~OPT_CSV);
    }
    else {
    }
    if(RF_ySize == 0) {
      RF_opt = RF_opt & (~OPT_PERF);
      RF_opt = RF_opt & (~OPT_VIMP);
      RF_optHigh              = RF_optHigh & (~OPT_CSE);
      RF_optHigh              = RF_optHigh & (~OPT_CSV);
      RF_opt = RF_opt & (~OPT_OENS);
      RF_opt = RF_opt & (~OPT_FENS);
    }
    if ( !(RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2)) { 
      RF_opt                  = RF_opt & (~OPT_PERF);
      RF_opt                  = RF_opt & (~OPT_VIMP);
      RF_optHigh              = RF_optHigh & (~OPT_CSE);
      RF_optHigh              = RF_optHigh & (~OPT_CSV);
      RF_opt                  = RF_opt & (~OPT_OENS);
      if (RF_opt & OPT_PROX) {
        RF_opt = RF_opt | OPT_PROX_IBG;
        RF_opt = RF_opt | OPT_PROX_OOB;
      }
      if (RF_optHigh & OPT_DIST) {
        RF_optHigh = RF_optHigh | OPT_DIST_IBG;
        RF_optHigh = RF_optHigh | OPT_DIST_OOB;
      }
      if (RF_optHigh & OPT_WGHT) {
        RF_optHigh = RF_optHigh | OPT_WGHT_IBG;
        RF_optHigh = RF_optHigh | OPT_WGHT_OOB;
      }
    }
    RF_opt = RF_opt & (~OPT_EMPR_RISK);
    if ((RF_opt & OPT_OENS) || (RF_opt & OPT_FENS)) {
    }
    else {
    }
    if (!(RF_opt & OPT_OENS)) {
      RF_opt                  = RF_opt & (~OPT_PERF);
      RF_optHigh              = RF_optHigh & (~OPT_CSE);
      RF_optHigh              = RF_optHigh & (~OPT_CSV);
    }
    break;
  }
  if ( !(RF_opt & OPT_BOOT_TYP1) && (RF_opt & OPT_BOOT_TYP2) ) {
    RF_optHigh = RF_optHigh & (~OPT_BOOT_SWOR);
  }
  if (RF_ptnCount > 0) {
      RF_optHigh = RF_optHigh | OPT_MEMB_PRUN;
      RF_opt     = RF_opt & (~OPT_PERF);
      RF_opt     = RF_opt & (~OPT_VIMP);
      RF_optHigh              = RF_optHigh & (~OPT_CSE);
      RF_optHigh              = RF_optHigh & (~OPT_CSV);
      RF_opt     = RF_opt & (~OPT_PROX);
      RF_opt     = RF_opt & (~OPT_OENS);
      RF_opt     = RF_opt & (~OPT_FENS);
      RF_optHigh = RF_optHigh & (~OPT_DIST);
      RF_optHigh = RF_optHigh & (~OPT_WGHT);
  }
  else {
    RF_optHigh = RF_optHigh & (~OPT_MEMB_PRUN);
  }
  if (RF_xMarginalSize > 0) {
    RF_opt = RF_opt & (~OPT_PERF);
    RF_opt = RF_opt & (~OPT_VIMP);
    RF_optHigh              = RF_optHigh & (~OPT_CSE);
    RF_optHigh              = RF_optHigh & (~OPT_CSV);
    RF_opt = RF_opt & (~OPT_OENS);
    RF_opt = RF_opt & (~OPT_FENS);
  }
  if (RF_opt & OPT_PERF) {
  }
  else {
    RF_opt                  = RF_opt & (~OPT_VIMP);
    RF_optHigh              = RF_optHigh & (~OPT_CSE);
    RF_optHigh              = RF_optHigh & (~OPT_CSV);
  }
}  
