package com.kogalur.randomforest;

import java.util.Random;

public class Simulation {

    public Simulation() {
    }

    
    public double[][] simmOne(int n, int p) {
        
        double[][] data;
        Random gen;
        
        data = new double[n][p];
        
        gen = new Random(19610804);
        
        for (int i=0; i<n; i++) {
            for (int j=0; j<p; j++) {
                data[i][j] = gen.nextGaussian(); 
            }
        }
        
        return data;
    }

    public double[][] regression(int n, int p, double[][] data) {

        double[]   coef;
        double[][] resp;
        double     result;
        
        coef = new double[p];
        resp = new double[n][1];

        for (int j=0; j<p; j++) {
            coef[j] = j+1;
        }

        for (int i=0; i<n; i++) {
            result = 0.0;
            for (int j=0; j<p; j++) {
                result = result + (coef[j] * data[i][j]);
            }
            resp[i][0] = result;
        }

        return resp;
    }
}

