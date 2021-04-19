/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testdata.h
 * Author: davidjwilliamson
 *
 * Created on 16 April 2021, 10:09
 */

#ifndef TESTDATA_H
#define TESTDATA_H
//
//struct TestData {
//    double Pos[3];
//    double OpticalDepth;
//};
//
//extern struct TestData *testpositions;
//
//extern int Nlocal,Ntot;

void generate_test_data(int N);
void dump_positions();

#endif /* TESTDATA_H */
