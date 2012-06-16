/// \file nist.cpp usage example of Statistical_tests class

#include "StatisticalTests.h"
#include "ProgressTimer.h"

#define INPUT_FILE_NAME "data/logistic_key.bin"
#define OUTPUT_FILE_NAME "res/logistic_key.txt"

using namespace std;

int main(void)
{
    int nb_bits = 1000000;
    int nb_seq = 10;
    int byte_len = 8*sizeof(char);
    int nb_bytes = nb_bits/byte_len;
    if (nb_bits%byte_len) nb_bytes++;

    StatisticalTests test(nb_bits);
    test.setInputFileName(INPUT_FILE_NAME);
    test.setOutputFileName(OUTPUT_FILE_NAME);

    ProgressTimer timer;
    timer.setMax(nb_seq);
    timer.progress(0.0);
    for (int n=0;n<nb_seq;n++)
    {
        test.setBinaryVector(n*nb_bytes);
        test.frequency();
        test.blockFrequency();
        test.cumulativeSums();
        test.runs();
        test.longestRunOfOnes();
        test.rank();
        test.discreteFourierTransform();
        test.nonOverlappingTemplateMatchings();
        test.overlappingTemplateMatchings();
        test.universal();
        test.approximateEntropy();
        test.randomExcursions();
        test.randomExcursionsVariant();
        test.serial();
        test.lempelZivCompression();
        test.linearComplexity();
        timer.progress(n+1);
    }
    timer.tocPrint();
    //tests.showBinaryVector();
}
