package com.kogalur.randomforest;

import com.kogalur.randomforest.RFLogger;
import java.util.logging.Level;

class NativeOpt {

    String name;

    java.util.HashMap <String, Integer> option;

    NativeOpt(String optionName, String[] optionParm, int[] optionBit) {

        if (optionParm.length != optionBit.length) {
            RFLogger.log(Level.SEVERE, "Constructor parameters malformed and of unequal length for:  " + optionName);
            RFLogger.log(Level.SEVERE, "Please contact technical support.");
            throw new IllegalArgumentException();
        }

        name = optionName;
        option = new java.util.HashMap <String, Integer> (optionParm.length);

        for (int i = 0; i < optionParm.length; i++) {
            option.put(optionParm[i], optionBit[i]);
        }
    }


    int get(String optionParm) {
        if (!option.containsKey(optionParm)) {
            RFLogger.log(Level.SEVERE, "Unknown key:  " + optionParm);
            throw new IllegalArgumentException();
        }
        return option.get(optionParm).intValue();
    }


}
