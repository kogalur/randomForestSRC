package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;
import java.util.logging.Level;

public final class Trace {

    public static final int USR = 2^0;
    public static final int LOW = 2^1;
    public static final int MED = 2^2;
    public  static final int HGH = 2^3;

    private static int trace;

    private Trace() {
        trace = 0;
    }
    
    public static void set(int val) {
        trace = val;
    }
    
    public static boolean get(int val) {
        if ((trace & val) > 0) {
            return true;
        }
        else {
            return false;
        }
    }
}
