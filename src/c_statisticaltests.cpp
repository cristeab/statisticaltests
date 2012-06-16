/// \file C_test.cpp Matlab interface for NIST tests

/*
 * Matlab interface for NIST tests (uses GSL library)
 * Matlab function declaration: [freq block_freq cum_sum runs longest_run rank dft ...
 *         non_overlap overlap univ app_en rand_ex rand_ex_v ...
 *         serial lz_compression complexity] = C_test(nb_bits, offset_bytes, input_file_name);
 * compilation command: mex -lgsl C_test.cpp
 */

#include <mex.h>
#include "StatisticalTests.h"

/*
 * The following functions are taken from IT++ library in order to ensure the interface
 * between C++ and Matlab.
 */
namespace itpp
{
    int mxArray2int(const mxArray *in)
    {
        int size;
        double* temp = (double*) mxGetPr(in);
        if (temp==0) mexErrMsgTxt("mxArray2int: Pointer to data is NULL");
        size = mxGetNumberOfElements(in);
        if (size!=1) mexErrMsgTxt("mxArray2int: Size of data is not equal to one");
        return (int) (*temp);
    }
    std::string mxArray2string(const mxArray *in)
    {
        if (in == 0) mexErrMsgTxt("mxArray2string: Pointer to data is NULL");
        std::string str = mxArrayToString(in);
        if (str.data() == 0) mexErrMsgTxt("mxArray2string: Could not convert mxArray to string");
        return str;
    }
    void Cvec2mxArray(double *in, mxArray *out)
    {
        double* temp = (double *) mxGetPr(out);
        if (temp==0) mexErrMsgTxt("Cvec2mxArray: Pointer to data is NULL");
        int size = mxGetNumberOfElements(out);
        if (size==0) mexErrMsgTxt("Cvec2mxArray: Size of data is zero");
        for(int i=0; i<size; i++) { *temp++= in[i]; }
    }
    void double2mxArray(const double &in, mxArray *out)
    {
        double* temp = (double *) mxGetPr(out);
        if (temp==0) mexErrMsgTxt("double2mxArray: Pointer to data is NULL"); 
        *temp= (double) in;
    }
}

// interface function between Matlab and C++ for NIST statistical tests
void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
{
    // Check the number of inputs and output arguments
    if(n_output!=16) mexErrMsgTxt("Sixteen outputs required!");
    if(n_input!=3) mexErrMsgTxt("Three inputs required!");
    
    // Convert input variables to IT++ format
    int nb_bits = itpp::mxArray2int(input[0]);
    int offset_bytes = itpp::mxArray2int(input[1]);
    std::string input_file_name = itpp::mxArray2string(input[2]);
    
    // ------------------ Start of routine ---------------------------
    StatisticalTests test(nb_bits);
    test.setInputFileName(input_file_name.c_str());
    double *P_value;
    
    test.setBinaryVector(offset_bytes);
    P_value = test.frequency();
    if (P_value)
    {
        output[0] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #0
        itpp::Cvec2mxArray(P_value, output[0]);// Convert the IT++ format to Matlab format for output #0
        delete[] P_value;
    }
    else
    {
        mexPrintf("Frequency test not applicable\n");
        output[0] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #0
        itpp::double2mxArray(0.0, output[0]);
    }
    
    P_value = test.blockFrequency();
    if (P_value)
    {
        output[1] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #1
        itpp::Cvec2mxArray(P_value, output[1]);// Convert the IT++ format to Matlab format for output #1
        delete[] P_value;
    }
    else
    {
        mexPrintf("BlockFrequency test not applicable\n");
        output[1] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #1
        itpp::double2mxArray(0.0, output[1]);
    }
    
    P_value = test.cumulativeSums();
    if (P_value)
    {
        output[2] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #2
        itpp::Cvec2mxArray(P_value, output[2]);// Convert the IT++ format to Matlab format for output #2
        delete[] P_value;
    }
    else
    {
        mexPrintf("CumulativeSums test not applicable\n");
        output[2] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #2
        itpp::double2mxArray(0.0, output[2]);
    }
    
    P_value = test.runs();
    if (P_value)
    {
        output[3] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #3
        itpp::Cvec2mxArray(P_value, output[3]);// Convert the IT++ format to Matlab format for output #3
        delete[] P_value;
    }
    else
    {
        mexPrintf("Runs test not applicable\n");
        output[3] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #3
        itpp::double2mxArray(0.0, output[3]);
    }
    
    P_value = test.longestRunOfOnes();
    if (P_value)
    {
        output[4] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #4
        itpp::Cvec2mxArray(P_value, output[4]);// Convert the IT++ format to Matlab format for output #4
        delete[] P_value;
    }
    else
    {
        mexPrintf("LongestRunOfOnes test not applicable\n");
        output[4] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #4
        itpp::double2mxArray(0.0, output[4]);
    }
    
    P_value = test.rank();
    if (P_value)
    {
        output[5] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #5
        itpp::Cvec2mxArray(P_value, output[5]);// Convert the IT++ format to Matlab format for output #5
        delete[] P_value;
    }
    else
    {
        mexPrintf("Rank test not applicable\n");
        output[5] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #5
        itpp::double2mxArray(0.0, output[5]);
    }
    
    P_value = test.discreteFourierTransform();
    if (P_value)
    {
        output[6] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #6
        itpp::Cvec2mxArray(P_value, output[6]);// Convert the IT++ format to Matlab format for output #6
        delete[] P_value;
    }
    else
    {
        mexPrintf("DiscreteFourierTransform test not applicable\n");
        output[6] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #6
        itpp::double2mxArray(0.0, output[6]);
    }
    
    P_value = test.nonOverlappingTemplateMatchings();
    if (P_value)
    {
        output[7] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #7
        itpp::Cvec2mxArray(P_value, output[7]);// Convert the IT++ format to Matlab format for output #7
        delete[] P_value;
    }
    else
    {
        mexPrintf("NonOverlappingTemplateMatchings test not applicable\n");
        output[7] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #7
        itpp::double2mxArray(0.0, output[7]);
    }
    
    P_value = test.overlappingTemplateMatchings();
    if (P_value)
    {
        output[8] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #8
        itpp::Cvec2mxArray(P_value, output[8]);// Convert the IT++ format to Matlab format for output #8
        delete[] P_value;
    }
    else
    {
        mexPrintf("OverlappingTemplateMatchings test not applicable\n");
        output[8] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #8
        itpp::double2mxArray(0.0, output[8]);
    }
    
    P_value = test.universal();
    if (P_value)
    {
        output[9] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #9
        itpp::Cvec2mxArray(P_value, output[9]);// Convert the IT++ format to Matlab format for output #9
        delete[] P_value;
    }
    else
    {
        mexPrintf("Universal test not applicable\n");
        output[9] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #9
        itpp::double2mxArray(0.0, output[9]);
    }
    
    P_value = test.approximateEntropy();
    if (P_value)
    {
        output[10] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #10
        itpp::Cvec2mxArray(P_value, output[10]);// Convert the IT++ format to Matlab format for output #10
        delete[] P_value;
    }
    else
    {
        mexPrintf("ApproximateEntropy test not applicable\n");
        output[10] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #10
        itpp::double2mxArray(0.0, output[10]);
    }
    
    P_value = test.randomExcursions();
    if (P_value)
    {
        output[11] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #11
        itpp::Cvec2mxArray(P_value, output[11]);// Convert the IT++ format to Matlab format for output #11
        delete[] P_value;
    }
    else
    {
        mexPrintf("RandomExcursions test not applicable\n");
        output[11] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #11
        itpp::double2mxArray(0.0, output[11]);
    }
    
    P_value = test.randomExcursionsVariant();
    if (P_value)
    {
        output[12] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #12
        itpp::Cvec2mxArray(P_value, output[12]);// Convert the IT++ format to Matlab format for output #12
        delete[] P_value;
    }
    else
    {
        mexPrintf("RandomExcursionsVariant test not applicable\n");
        output[12] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #12
        itpp::double2mxArray(0.0, output[12]);
    }
    
    P_value = test.serial();
    if (P_value)
    {
        output[13] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #13
        itpp::Cvec2mxArray(P_value, output[13]);// Convert the IT++ format to Matlab format for output #13
        delete[] P_value;
    }
    else
    {
        mexPrintf("RandomExcursionsVariant test not applicable\n");
        output[13] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #13
        itpp::double2mxArray(0.0, output[13]);
    }
    
    P_value = test.lempelZivCompression();
    if (P_value)
    {
        output[14] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #14
        itpp::Cvec2mxArray(P_value, output[14]);// Convert the IT++ format to Matlab format for output #14
        delete[] P_value;
    }
    else
    {
        mexPrintf("LempelZivCompression test not applicable\n");
        output[14] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #14
        itpp::double2mxArray(0.0, output[14]);
    }
    
    P_value = test.linearComplexity();
    if (P_value)
    {
        output[15] = mxCreateDoubleMatrix(1, test.getNbPValues(), mxREAL);// Create vector output #15
        itpp::Cvec2mxArray(P_value, output[15]);// Convert the IT++ format to Matlab format for output #15
        delete[] P_value;
    }
    else
    {
        mexPrintf("LinearComplexity test not applicable\n");
        output[15] = mxCreateDoubleMatrix(1, 1, mxREAL);// Create vector output #15
        itpp::double2mxArray(0.0, output[15]);
    }
    // ------------------ End of routine -----------------------------
}
