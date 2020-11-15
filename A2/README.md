# Assignment 2 

This program takes a reference genome file (ref_genome)with length of 1000000 characters and a sequnecing reads file (e.g. reads_0)with
length of 20000000 characters then recovers the underlying genome with a output file (e.g. assembly_0).


In order to run the program, you need to install Succinct Data Structure Library (SDSL) first. URL: https://github.com/simongog/sdsl-lite

Complie the program: g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib Aligner.cpp -o Aligner -lsdsl -ldivsufsort -ldivsufsort64
Run the program: ./Aligner


There is a accuracy test program for testing the result.

Complie the test program: g++ Accuracy_test.cpp -o Accuracy_test
Run the test program: ./Accuracy_test refer_genome assembly_0