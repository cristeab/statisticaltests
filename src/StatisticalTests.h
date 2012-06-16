/// \file Statistical_tests.cpp Statistical_tests class

/** \mainpage A Statistical Test Suite for Pseudoradom Number Generators used in Cryptographic Applications
 *
 * \author Bogdan Cristea
 * \date June 16, 2012
 * \version 1.1
 * 
 * \section Description
 * The Statistical Tests Suite was proposed by the National Institute of Standards
 * and Technology (NIST) in 2001 and implemented in C.
 * The StatisticalTests class reimplements all tests in C++ using the code from
 * NIST (v1.5) and the GNU Scientific Library (GSL).
 * 
 * Each statistical test is implemented as a member function whose output is
 * the vector of computed P-values.
 * Computation details of each test can be written into an ASCII file.
 * 
 * A number of \f$16\f$ statistical tests are implemented:
 * -# Frequency (Monobit) Test (StatisticalTests::Frequency)
 * -# Frequency Test within a Block (StatisticalTests::BlockFrequency)
 * -# Cumulative Sums (Cusum) Test (StatisticalTests::CumulativeSums)
 * -# Runs Test (StatisticalTests::Runs)
 * -# Test for the Longest Run of Ones in a Block (StatisticalTests::LongestRunOfOnes)
 * -# Binary Matrix Rank Test (StatisticalTests::Rank)
 * -# Discrete Fourier Transform (Spectral) Test (StatisticalTests::DiscreteFourierTransform)
 * -# Non-overlapping Template Matching Test (StatisticalTests::NonOverlappingTemplateMatchings)
 * -# Overlapping Template Matching Test (StatisticalTests::OverlappingTemplateMatchings)
 * -# Maurer's ''Universal Statistical`` Test (StatisticalTests::Universal)
 * -# Approximate Entropy Test (StatisticalTests::ApproximateEntropy)
 * -# Random Excursions Test (StatisticalTests::RandomExcursions)
 * -# Random Excursions Variant Test (StatisticalTests::RandomExcursionsVariant)
 * -# Serial Test (StatisticalTests::Serial)
 * -# Lempel-Ziv Compression Test (StatisticalTests::LempelZivCompression)
 * -# Linear Complexity Test (StatisticalTests::LinearComplexity)
 * 
 * For further details see \ref label_details "implementation issues".
 * 
 * \section label_usage Usage
 * In order to use StatisticalTests class, the GSL must be installed.
 * For example, assuming that the GSL is installed in your path, the compilation
 * command used on an Athlon64 machine is:
 * 
 *  g++ -Wall -O3 -pipe -march=athlon64 -lgsl nist.cpp -o nist
 * 
 * \section label_interface Matlab interface
 * It is possible to call StatisticalTests class member functions from Matlab
 * using an interface function between Matlab and C++
 * written in C++. The Matlab function must be called as follows:
 * 
 * [freq block_freq cum_sum runs longest_run rank dft non_overlap overlap univ
 * app_en rand_ex rand_ex_v serial lz_compression complexity] = C_test(nb_bits, offset_bytes, input_file_name);
 *  
 * The compilation command for the interface function is: 
 * 
 * mex -lgsl C_test.cpp
 */

/// Statistical Tests class

/** 
 * Each statistical test is implemented as a member function whose output is the
 * vector of computed P-values. Computation details
 * of each test can be written into an ASCII file.
 * 
 * \anchor label_details
 * The obtained results were compared against the results obtained with the program
 * from NIST using a modified version of the function
 * who reads bits from binary file. In our implementation, the least significative
 * bit (LSB) is read first as opposed to the
 * original function of NIST who reads first the most significative bit (MSB).
 * Also the original function was changed to correctly
 * take into account the size of "unsigned long" type.
 * 
 * The following tests give different results with respect to the original program:
 * - StatisticalTests::DiscreteFourierTransform (the length of the FFT transform is incorrect)
 * - StatisticalTests::RandomExcursions (when processing the cycles, cycleStart and
 * cycleStop variables are incorrectly updated)
 */

#ifndef STATISTICAL_TESTS
#define STATISTICAL_TESTS

#include <iostream>
#include <fstream>

#include <cmath>
#include <gsl/gsl_sf_erf.h> //for erfc
#include <gsl/gsl_sf_gamma.h> //for gamma_inc_Q
#include <gsl/gsl_cdf.h> //for ugaussian_P (standard normal cumulative probability distribution function)
#include <gsl/gsl_fft_real.h>

class StatisticalTests
{
public:
	/// %StatisticalTests constructor
	/** Internal variables are initialized:
	 * - number of bits is set (default 1 million bits)
	 * - memory is reserved for the vector of bits
	 * - no computation details
	 * - significance level is set to 0.01
	 * - test result vector of strings is initialized with "FAILURE" and "SUCCESS"
	 */
	explicit StatisticalTests(int nb=1000000) :
        input_file_name(NULL),
		nb_bits(nb),
		alpha(0.01),
        nb_P_values(0),
		epsilon(new Bit[nb_bits]),
		details(false)
    {
		test_result[0] = "FAILURE";// P_value<alpha
		test_result[1] = "SUCCESS";// P_value>=alpha
	}
	/// Sets input file name
	/**
	 * The input file stores the binary sequences as binary characters, which
	 * are then read by StatisticalTests::set_binary_vector.
	 */
	void setInputFileName(const char* ifn)
	{
		input_file_name = ifn;
	}
	/// Sets output file name
	/**
	 * Test results are written to output file in ASCII format.
	 */
	void setOutputFileName(const char* output_file_name)
	{
		output_stream.open(output_file_name);
		if (output_stream==NULL)
		{
			std::cout << "set_output_file_name: cannot open file for writing.\n";
			return;
		}
		details = true;
	}
	/// Sets binary vector
	/**
	 * Its input represents the offset, in bytes, which is used to read bits from
	 * file. Default is \f$0\f$ (the beginning
	 * of the file). This allows to subsequently read from input file and then
	 * test several binary sequences of the same length.
	 */
	void setBinaryVector(int nb_offset_bytes=0);
	/// Shows binary vector
	/**
	 * Prints binary vector in ASCII format to output file. This function could
	 * be useful to visualize the binary sequence
	 * after it is initialized by StatisticalTests::set_binary_vector.
	 */
	void showBinaryVector(void)
	{
		for (int n=0;n<nb_bits;n++)
			output_stream << (epsilon[n].b?'1':'0');
		output_stream << std::endl;
	}
	/// Sets details
	/**
	 * When its input is true, computation details are written to output file
	 * (whose name is set by StatisticalTests::set_output_file_name)
	 * By default, no computation details are given.
	 */
	void setDetails(bool in)
	{
		details = in;
	}
	/// Sets significance level
	/**
	 * If P-value >= alpha (significance level) then the binary sequence is considered random.
	 */
	void setSignificanceLevel(double a)
	{
		alpha = a;
	}
	/// Gets the number of P-values
	/**
	 * Allows to obtain the number P-values returned by a test function applied
	 * to the binary sequence. Thus further processing
	 * can be applied to the obtained vector of P-values.
	 */
	int getNbPValues(void)
	{
		return nb_P_values;
	};
	/// Frequency (Monobit) Test
	/**
	 * Test the proportion of zeroes and ones for the entire sequence. The purpose
	 * of this test is to determine whether the number of ones and zeros in a
	 * sequence are approximately the same as would be expected for a truly random
	 * sequence. The test assesses the closeness of the fraction of ones to \f$1/2\f$,
	 * that is, the number of ones and zeroes in a sequence should be about the same.
	 * All subsequent tests depend on the passing of this test; there is no
	 * evidence to indicate that the tested sequence is non-random.
	 */
	double* frequency(void);
	/// Frequency Test within a Block
	/**
	 * The focus of the test is the proportion of ones within \f$block\_len\f$-bit
	 * blocks. The purpose of this test is to
	 * determine whether the frequency of ones in an \f$block\_len\f$-bit block
	 * is approximately \f$block\_len/2\f$, as would be
	 * expected under an assumption of randomness. For block size \f$block\_len=1\f$,
	 * this test degenerates to the Frequency (Monobit) test.
	 */
	double* blockFrequency(int block_len=20000);
	/// Cumulative Sums (Cusum) Test
	/**
	 * The focus of this test is the maximal excursion (from zero) of the random
	 * walk defined by the cumulative sum of adjusted (\f$-1\f$, \f$+1\f$) digits
	 * in the sequence. The purpose of the test is to determine whether the cumulative
	 * sum of the partial sequences occurring in the tested sequence is too large
	 * or too small relative to the expected behavior of that cumulative sum for random
	 * sequences. This cumulative sum may be considered as a random walk. For a
	 * random sequence, the excursions of the random walk should be near zero.
	 * For certain types of non-random sequences, the excursions of this random
	 * walk from zero will be large. The test is applied forward and backward
	 * through the sequence (\f$2\f$ P-values are generated).
	 */
	double* cumulativeSums(void);
	/// Runs Test
	/**
	 * The focus of this test is the total number of runs in the sequence, where
	 * a run is an uninterrupted sequence of identical bits. A run of length \f$k\f$
	 * consists of exactly \f$k\f$ identical bits and is bounded before and after
	 * with a bit of the opposite value. The purpose of the runs test is to determine
	 * whether the number of runs of ones and zeros of various lengths is as expected
	 * for a random sequence. In particular, this test determines whether the
	 * oscillation between such zeros and ones is too fast or too slow.
	 */
	double* runs(void);
	/// Test for the Longest Run of Ones in a Block
	/**
	 * The focus of the test is the longest run of ones within \f$long\_runs\_case\f$-bit
	 * blocks. The purpose of this test is to determine whether the length of
	 * the longest run of ones within the tested sequence is consistent
	 * with the length of the longest run of ones that would be expected in a
	 * random sequence. Note that an irregularity in the expected length of the
	 * longest run of ones implies that there is also an irregularity in the expected
	 * length of the longest run of zeroes. Therefore, only a test for ones is
	 * necessary. The test code has been pre-set to accommodate three
	 * values for \f$long\_runs\_case\f$: \f$long\_runs\_case = 8\f$, \f$long\_runs\_case = 128\f$
	 * and \f$long\_runs\_case = 10^4\f$ in accordance with the number of bits
	 * (\f$128\f$, \f$6272\f$ and \f$750000\f$ respectively).
	 */
	double* longestRunOfOnes(int long_runs_case=10000);
	/// Binary Matrix Rank Test
	/**
	 * The focus of the test is the rank of disjoint sub-matrices of the entire
	 * sequence. The purpose of this test is to check for linear dependence among
	 * fixed length substrings of the original sequence. The matrix size is set
	 * to \f$32\times 32\f$.
	 */
	double* rank(void);
	/// Discrete Fourier Transform (Spectral) Test
	/**
	 * The focus of this test is the peak heights in the Discrete Fourier Transform
	 * of the sequence. The purpose of this test is to detect periodic features
	 * (i.e., repetitive patterns that are near each other) in the tested sequence
	 * that would indicate a deviation from the assumption of randomness. The
	 * intention is to detect whether the number of peaks exceeding the \f$95\%\f$
	 * threshold is significantly different than \f$5\%\f$. The Discrete Fourier
	 * Transform size equals the binary sequence length.
	 */
	double* discreteFourierTransform(void);
	/// Non-overlapping Template Matching Test
	/**
	 * The focus of this test is the number of occurrences of pre-specified target
	 * strings. The purpose of this test is to detect generators that produce too
	 * many occurrences of a given non-periodic (aperiodic) pattern. For this test
	 * and for the Overlapping Template Matching test, an \f$template\_len\f$-bit
	 * window is used to search for a specific \f$template\_len\f$-bit pattern.
	 * If the pattern is not found, the window slides one bit position. If the pattern
	 * is found, the window is reset to the bit after the found pattern, and the
	 * search resumes. It is recommended that \f$template\_len = 9\f$ or
	 * \f$template\_len = 10\f$ be specified to obtain meaningful results.
	 */
	double* nonOverlappingTemplateMatchings(int template_len=10);
	/// Overlapping Template Matching Test
	/**
	 * The focus of the Overlapping Template Matching test is the number of occurrences
	 * of prespecified target strings. Both this test and the Non-overlapping
	 * Template Matching test use an \f$template\_len\f$-bit window to search for
	 * a specific \f$template\_len\f$-bit pattern. As with the Non-overlapping Template
	 * Matching test, if the pattern is not found, the window slides one bit position.
	 * The difference between this test and the Non-overlapping Template Matching
	 * test is that when the pattern is found, the window slides
	 * only one bit before resuming the search. Various values of \f$template\_len\f$
	 * may be selected, but for the time being, NIST recommends \f$template\_len = 9\f$
	 * or \f$template\_len = 10\f$.
	 */
	double* overlappingTemplateMatchings(int template_len=10);
	/// Maurer's ''Universal Statistical`` Test
	/**
	 * The focus of this test is the number of bits between matching patterns
	 * (a measure that is related to the length of a compressed sequence).
	 * The purpose of the test is to detect whether or not the sequence can be
	 * significantly compressed without loss of information. A significantly
	 * compressible sequence is considered to be non-random. This test requires
	 * a long sequence of bits (\f$nb\_bits \geq (Q + K)L\f$) which are divided
	 * into two segments of \f$L\f$-bit blocks, where \f$L\f$ should be chosen so
	 * that \f$6 \leq L \leq 16\f$. The first segment consists of \f$Q\f$
	 * initialization blocks, where \f$Q\f$ should be chosen so that \f$Q = 10\cdot 2^L\f$.
	 * The second segment consists of \f$K\f$ test blocks, where \f$K = \lceil L\rceil - Q \approx 1000\cdot 2^L\f$.
	 * \f$L\f$ is chosen using \f$nb\_bits\f$ and then \f$Q\f$ is computed.
	 */
	double* universal(void);
	/// Approximate Entropy Test
	/**
	 * As with the Serial test, the focus of this test is the frequency of all
	 * possible overlapping \f$block\_len\f$-bit patterns across the entire sequence.
	 * The purpose of the test is to compare the frequency of overlapping blocks
	 * of two consecutive/adjacent lengths (\f$block\_len\f$ and \f$block\_len+1\f$) against the
	 * expected result for a random sequence. Choose \f$block\_len\f$ and \f$nb\_bits\f$
	 * such that \f$block\_len < \lfloor log_2 (nb\_bits)\rfloor -2\f$.
	 * As the test is implemented, the default value is the recommended value
	 * for \f$1\f$ milion bits.
	 */
	double* approximateEntropy(int first_block_len=14);
	/// Random Excursions Test
	/**
	 * The focus of this test is the number of cycles having exactly \f$K\f$ visits
	 * in a cumulative sum random walk. The cumulative sum random walk is derived
	 * from partial sums after the \f$(0,1)\f$ sequence is transferred to the
	 * appropriate \f$(-1, +1)\f$ sequence. A cycle of a random walk consists
	 * of a sequence of steps of unit length taken at random that begin at and
	 * return to the origin. The purpose of this test is to determine if the number
	 * of visits to a particular state within a cycle deviates from what one would
	 * expect for a random sequence. This test is actually a series of
	 * eight tests (and conclusions), one test and conclusion for each of the states:
	 * \f$-4\f$, \f$-3\f$, \f$-2\f$, \f$-1\f$ and \f$+1\f$,
	 * \f$+2\f$, \f$+3\f$, \f$+4\f$.
	 */
	double* randomExcursions(void);
	/// Random Excursions Variant Test
	/**
	 * The focus of this test is the total number of times that a particular
	 * state is visited (i.e., occurs) in a cumulative sum random walk. The purpose
	 * of this test is to detect deviations from the expected number of visits to
	 * various states in the random walk. This test is actually a series of
	 * eighteen tests (and conclusions), one test and conclusion for each of
	 * the states: \f$-9\f$, \f$-8\f$, ..., \f$-1\f$ and
	 * \f$+1\f$, \f$+2\f$, ..., \f$+9\f$.
	 */
	double* randomExcursionsVariant(void);
	/// Serial Test
	/**
	 * The focus of this test is the frequency of all possible overlapping
	 * \f$pattern\_len\f$-bit patterns across the entire sequence. The purpose of
	 * this test is to determine whether the number of occurrences of the \f$2^{pattern\_len}\f$
	 * \f$pattern\_len\f$-bit overlapping patterns is approximately the same as
	 * would be expected for a random sequence. Random sequences have uniformity;
	 * that is, every \f$pattern\_len\f$-bit pattern has the same chance
	 * of appearing as every other \f$pattern\_len\f$-bit pattern. Note that for
	 * \f$pattern\_len = 1\f$, the Serial test is equivalent to the
	 * Frequency test. Choose \f$pattern\_len\f$ and \f$nb\_bits\f$ such that
	 * \f$pattern\_len < \lfloor log_2 (nb\_bits)\rfloor -2\f$.
	 * The default value is chosen for \f$1\f$ milion bits.
	 */
	double* serial(int pattern_len=16);
	/// Lempel-Ziv Compression Test
	/**
	 * The focus of this test is the number of cumulatively distinct patterns
	 * (words) in the sequence. The purpose of the test is to determine how far
	 * the tested sequence can be compressed. The sequence is considered to be
	 * non-random if it can be significantly compressed. A random
	 * sequence will have a characteristic number of distinct patterns.
	 * Note that since no known theory is available to determine the exact
	 * values of \f$\mu\f$ and \f$\sigma^2\f$ (needed to compute the P-value), these
	 * values were computed (under an assumption of randomness) using SHA-1. The Blum-
	 * Blum-Shub generator will give similar results for \f$\mu\f$ and \f$\sigma^2\f$.
	 * The code adapts automatically the number of bits in order to use a value
	 * for which \f$\mu\f$ and \f$\sigma^2\f$ are precomputed.
	 */
	double* lempelZivCompression(void);
	/// Linear Complexity Test
	/**
	 * The focus of this test is the length of a linear feedback shiftregister (LFSR).
	 * The purpose of this test is to determine whether or not the sequence is
	 * complex enough to be considered random. Random sequences are characterized
	 * by longer LFSRs. An LFSR that is too short implies nonrandomness.
	 * The value of \f$block\_len\f$ must be in the range \f$500\leq block\_len\leq 5000\f$,
	 * and \f$nb\_blocks\geq 200\f$ for the \f$\chi^2\f$ result to be valid. Note
	 * that for 1 milion bits, the default value \f$block\_len=5000\f$ implies
	 * that \f$nb\_blocks=200\f$.
	 */
	double* linearComplexity(int block_len=5000);
	/// %StatisticalTests destructor
	/**
	 * Memory is freed for the vector of bits and the output file stream is closed
	 */
	~StatisticalTests()
	{
		delete[] epsilon;
		output_stream.close();
	}
private:
	/// input file name from who the binary sequence is read
	const char* input_file_name;
	/// ouput file stream where the test results are written
	std::ofstream output_stream;
	/// number of bits
	int nb_bits;
	/// significance level
	double alpha;
	/// test result strings (FAILURE or SUCCESS)
	const char* test_result[2];
	/// number of P-values generated by each test
	int nb_P_values;
	/// %bit structure used to store the binary sequence
	struct Bit
	{
		unsigned char b:1;
	}*epsilon;// BIT FIELD STRUCTURE
	/// if true computation details are written to output file
	bool details;
	/// internal function of Rank test
	void defMatrix(Bit *matrix, int k);
	/// internal function of Rank test
	int computeRank(Bit* matrix);
	/// internal function of Rank test
	void performElementaryRowOperations(Bit* A, int flag, int i);
	/// internal function of Rank test
	int findUnitElementAndSwap(Bit* A, int flag, int i);
	/// internal function of Rank test
	int swapRows(Bit* A, int i, int index);
	/// internal function of Rank test
	int determineRank(const Bit* A);
	/// internal function of Rank test
	bool generateNonperiodicTemplate(int index, Bit* pattern, int template_len, 
            unsigned int* workspace, int len); //internal function of NonOverlappingTemplateMatchings
	/// internal function of OverlappingTemplateMatchings test
	double Pr(int u, double eta);
	/// internal function of Serial test
	double psi2(int m);
};

#endif
