/// \file nist.cpp usage example of Statistical_tests class

#include "Statistical_tests.cpp"
#include "../Timer/Progress_class.cpp"

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

    Statistical_Tests test(nb_bits);
    test.set_input_file_name(INPUT_FILE_NAME);
    test.set_output_file_name(OUTPUT_FILE_NAME);

    Progress_Timer timer;
    timer.set_max(nb_seq);
    timer.progress(0.0);
    for (int n=0;n<nb_seq;n++)
    {
        test.set_binary_vector(n*nb_bytes);
        test.Frequency();
        test.BlockFrequency();
        test.CumulativeSums();
        test.Runs();
        test.LongestRunOfOnes();
        test.Rank();
        test.DiscreteFourierTransform();
        test.NonOverlappingTemplateMatchings();
        test.OverlappingTemplateMatchings();
        test.Universal();
        test.ApproximateEntropy();
        test.RandomExcursions();
        test.RandomExcursionsVariant();
        test.Serial();
        test.LempelZivCompression();
        test.LinearComplexity();
        timer.progress(n+1);
    }
    timer.toc_print();
    //tests.show_binary_vector();
}
