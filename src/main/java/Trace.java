package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;
import java.util.logging.Level;

final class Trace {

    static final int USR = 2^0;
    static final int LOW = 2^1;
    static final int MED = 2^2;
    static final int HGH = 2^3;

    private static int trace;

    private Trace() {
        trace = 0;
    }
    
    static void set(int val) {
        trace = val;
    }
    
    static boolean get(int val) {
        if ((trace & val) > 0) {
            return true;
        }
        else {
            return false;
        }
    }
}
