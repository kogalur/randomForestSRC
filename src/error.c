
// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***
#include           "global.h"
#include           "external.h"

// *** THIS HEADER IS AUTO GENERATED. DO NOT EDIT IT ***

      
    

#include "error.h"
const char vomit[] = "\nRF-SRC:  The application will now exit.\n";
void exit2R(void) {
  Rprintf("%s", vomit);
  error(NULL);
}
void printR(char *format, ...) {
  char *buffer;
  va_list aptr;
  buffer = (char *) malloc(sizeof(char) * 1023);
  va_start(aptr, format);
  vsnprintf(buffer, sizeof(char) * 1023, format, aptr);
  va_end(aptr);
  Rprintf("%s", buffer);
  free((char *) buffer);
}
