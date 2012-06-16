/// \file Statistical_tests.cpp Statistical_tests class

/** \mainpage A Statistical Test Suite for Pseudoradom Number Generators used in Cryptographic Applications
 *
 * \author Bogdan Cristea
 * \date November 18, 2007
 * \version 1.1
 * 
 * \section Description
 * The Statistical Tests Suite was proposed by the National Institute of Standards
 * and Technology (NIST) in 2001 and implemented in C.
 * The Statistical_Tests class reimplements all tests in C++ using the code from
 * NIST (v1.5) and the GNU Scientific Library (GSL).
 * 
 * Each statistical test is implemented as a member function whose output is
 * the vector of computed P-values.
 * Computation details of each test can be written into an ASCII file.
 * 
 * A number of \f$16\f$ statistical tests are implemented:
 * -# Frequency (Monobit) Test (Statistical_Tests::Frequency)
 * -# Frequency Test within a Block (Statistical_Tests::BlockFrequency)
 * -# Cumulative Sums (Cusum) Test (Statistical_Tests::CumulativeSums)
 * -# Runs Test (Statistical_Tests::Runs)
 * -# Test for the Longest Run of Ones in a Block (Statistical_Tests::LongestRunOfOnes)
 * -# Binary Matrix Rank Test (Statistical_Tests::Rank)
 * -# Discrete Fourier Transform (Spectral) Test (Statistical_Tests::DiscreteFourierTransform)
 * -# Non-overlapping Template Matching Test (Statistical_Tests::NonOverlappingTemplateMatchings)
 * -# Overlapping Template Matching Test (Statistical_Tests::OverlappingTemplateMatchings)
 * -# Maurer's ''Universal Statistical`` Test (Statistical_Tests::Universal)
 * -# Approximate Entropy Test (Statistical_Tests::ApproximateEntropy)
 * -# Random Excursions Test (Statistical_Tests::RandomExcursions)
 * -# Random Excursions Variant Test (Statistical_Tests::RandomExcursionsVariant)
 * -# Serial Test (Statistical_Tests::Serial)
 * -# Lempel-Ziv Compression Test (Statistical_Tests::LempelZivCompression)
 * -# Linear Complexity Test (Statistical_Tests::LinearComplexity)
 * 
 * For further details see \ref label_details "implementation issues".
 * 
 * \section label_usage Usage
 * In order to use Statistical_Tests class, the GSL must be installed.
 * For example, assuming that the GSL is installed in your path, the compilation
 * command used on an Athlon64 machine is:
 * 
 *  g++ -Wall -O3 -pipe -march=athlon64 -lgsl nist.cpp -o nist
 * 
 * \section label_interface Matlab interface
 * It is possible to call Statistical_Tests class member functions from Matlab
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
 * - Statistical_Tests::DiscreteFourierTransform (the length of the FFT transform is incorrect)
 * - Statistical_Tests::RandomExcursions (when processing the cycles, cycleStart and
 * cycleStop variables are incorrectly updated)
 */

#ifndef STATISTICAL_TESTS
#define STATISTICAL_TESTS

#include <iostream>
#include <fstream>

#include <cmath>
#include "gsl/gsl_sf_erf.h" //for erfc
#include "gsl/gsl_sf_gamma.h" //for gamma_inc_Q
#include "gsl/gsl_cdf.h" //for ugaussian_P (standard normal cumulative probability distribution function)
#include "gsl/gsl_fft_real.h"

class Statistical_Tests
{
public:
	/// %Statistical_Tests constructor
	/** Internal variables are initialized:
	 * - number of bits is set (default 1 million bits)
	 * - memory is reserved for the vector of bits
	 * - no computation details
	 * - significance level is set to 0.01
	 * - test result vector of strings is initialized with "FAILURE" and "SUCCESS"
	 */
	Statistical_Tests(const int &nb=1000000)
	{
		nb_bits = nb;// number of bits
		epsilon = new Bit[nb_bits];// vector of bits
		details = false;// if true write test details to file
		alpha = 0.01;// significance level
		test_result[0] = "FAILURE";// P_value<alpha
		test_result[1] = "SUCCESS";// P_value>=alpha
	};
	/// Sets input file name
	/**
	 * The input file stores the binary sequences as binary characters, which
	 * are then read by Statistical_Tests::set_binary_vector.
	 */
	void set_input_file_name(const char* ifn)
	{
		input_file_name = ifn;
	};
	/// Sets output file name
	/**
	 * Test results are written to output file in ASCII format.
	 */
	void set_output_file_name(const char* output_file_name)
	{
		output_stream.open(output_file_name);
		if (output_stream==NULL)
		{
			std::cout << "set_output_file_name: cannot open file for writing.\n";
			return;
		}
		details = true;
	};
	/// Sets binary vector
	/**
	 * Its input represents the offset, in bytes, which is used to read bits from
	 * file. Default is \f$0\f$ (the beginning
	 * of the file). This allows to subsequently read from input file and then
	 * test several binary sequences of the same length.
	 */
	void set_binary_vector(const int &nb_offset_bytes=0);
	/// Shows binary vector
	/**
	 * Prints binary vector in ASCII format to output file. This function could
	 * be useful to visualize the binary sequence
	 * after it is initialized by Statistical_Tests::set_binary_vector.
	 */
	void show_binary_vector(void)
	{
		for (int n=0;n<nb_bits;n++)
			output_stream << (epsilon[n].b?'1':'0');
		output_stream << std::endl;
	};
	/// Sets details
	/**
	 * When its input is true, computation details are written to output file
	 * (whose name is set by Statistical_Tests::set_output_file_name)
	 * By default, no computation details are given.
	 */
	void set_details(const bool &in)
	{
		details = in;
	}
	/// Sets significance level
	/**
	 * If P-value >= alpha (significance level) then the binary sequence is considered random.
	 */
	void set_significance_level(const double &a)
	{
		alpha = a;
	};
	/// Gets the number of P-values
	/**
	 * Allows to obtain the number P-values returned by a test function applied
	 * to the binary sequence. Thus further processing
	 * can be applied to the obtained vector of P-values.
	 */
	int get_nb_P_values(void)
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
	double* Frequency(void);
	/// Frequency Test within a Block
	/**
	 * The focus of the test is the proportion of ones within \f$block\_len\f$-bit
	 * blocks. The purpose of this test is to
	 * determine whether the frequency of ones in an \f$block\_len\f$-bit block
	 * is approximately \f$block\_len/2\f$, as would be
	 * expected under an assumption of randomness. For block size \f$block\_len=1\f$,
	 * this test degenerates to the Frequency (Monobit) test.
	 */
	double* BlockFrequency(const int &block_len=20000);
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
	double* CumulativeSums(void);
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
	double* Runs(void);
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
	double* LongestRunOfOnes(const int &long_runs_case=10000);
	/// Binary Matrix Rank Test
	/**
	 * The focus of the test is the rank of disjoint sub-matrices of the entire
	 * sequence. The purpose of this test is to check for linear dependence among
	 * fixed length substrings of the original sequence. The matrix size is set
	 * to \f$32\times 32\f$.
	 */
	double* Rank(void);
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
	double* DiscreteFourierTransform(void);
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
	double* NonOverlappingTemplateMatchings(const int &template_len=10);
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
	double* OverlappingTemplateMatchings(const int &template_len=10);
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
	double* Universal(void);
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
	double* ApproximateEntropy(const int &first_block_len=14);
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
	double* RandomExcursions(void);
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
	double* RandomExcursionsVariant(void);
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
	double* Serial(const int &pattern_len=16);
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
	double* LempelZivCompression(void);
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
	double* LinearComplexity(const int &block_len=5000);
	/// %Statistical_Tests destructor
	/**
	 * Memory is freed for the vector of bits and the output file stream is closed
	 */
	~Statistical_Tests()
	{
		delete[] epsilon;
		output_stream.close();
	};
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
	void def_matrix(Bit *matrix, const int &k);
	/// internal function of Rank test
	int computeRank(Bit* matrix);
	/// internal function of Rank test
	void perform_elementary_row_operations(Bit* A, const int &flag, const int &i);
	/// internal function of Rank test
	int find_unit_element_and_swap(Bit* A, const int &flag, const int &i);
	/// internal function of Rank test
	int swap_rows(Bit* A, const int &i, const int &index);
	/// internal function of Rank test
	int determine_rank(const Bit* A);
	/// internal function of Rank test
	bool generate_nonperiodic_template(const int &index, Bit* pattern,
			const int &template_len, unsigned int* workspace, const int &len); //internal function of NonOverlappingTemplateMatchings
	/// internal function of OverlappingTemplateMatchings test
	double Pr(const int &u, const double &eta);
	/// internal function of Serial test
	double psi2(const int &m);
};

//start reading beginning with the LSB
void Statistical_Tests::set_binary_vector(const int &nb_offset_bytes)
{
	std::ifstream input_stream(input_file_name, std::ios::binary);
	if (input_stream==NULL)
	{
		std::cout << "set_binary_vector: cannot open file for reading.\n";
		return;
	}
	input_stream.seekg(nb_offset_bytes, std::ios::beg);
	char temp;
	int byte_len = 8*sizeof(char);
	char mask = 0x1;//LSB
	register int n,k;
	for (n=0;n<nb_bits/byte_len;n++)
	{
		input_stream.get(temp);
		for (k=0;k<byte_len;k++)
		{
			epsilon[k+byte_len*n].b = ((temp&mask)?1:0);
			temp >>= 1;
		}
		if (input_stream.eof())
		{
			nb_bits = byte_len*n;
			output_stream << "set_binary_vector: eof reached. " << nb_bits << " bits were read.\n";
			input_stream.close();
			return;
		}
	}
	if (nb_bits%byte_len) //read last byte if any
	{
		input_stream.get(temp);
		int remaining_bits = nb_bits-byte_len*(nb_bits/byte_len);
		for (n=(nb_bits-remaining_bits);n<nb_bits;n++)
		{
			epsilon[n].b = ((temp&mask)?1:0);
			temp >>= 1;
		}
	}
	input_stream.close();
}

double* Statistical_Tests::Frequency(void)
{
	register double sum = 0.0;
	register int i;
	for (i=0;i<nb_bits;i++)
		sum += 2*(double)epsilon[i].b-1;
	double s_obs = fabs(sum)/sqrt(nb_bits);
	nb_P_values = 1;
	double* P_value = new double[1];
	P_value[0] = gsl_sf_erfc(s_obs/sqrt(2.0));

	//show computation details
	if (details)
	{
		int nb_ones;
		int nb_zeros;
		if (sum>=0)
		{
			nb_ones = nb_bits/2+int(sum);
			nb_zeros = nb_bits-nb_ones;
		}
		else
		{
			nb_zeros = nb_bits/2+int(sum);
			nb_ones = nb_bits-nb_zeros;
		}
		output_stream << "Frequency test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "Nb. of ones: " << nb_ones << std::endl;
		output_stream << "Nb. of zeros: " << nb_zeros << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::BlockFrequency(const int &block_len)
{
	int N = nb_bits/block_len; 		// # OF SUBSTRING BLOCKS
	double sum = 0.0;
	double pi,blockSum;
	register int i,j;
	for (i=0;i<N;i++)                 // N=10000 FOR EACH SUBSTRING BLOCK
	{
		pi = 0.0;
		blockSum = 0.0;
		for (j=0;j<block_len;j++)      // m=100 COMPUTE The "i"th Pi Value
			blockSum += (double)epsilon[j+i*block_len].b;
		pi = blockSum/(double)block_len;
		sum += pow(pi-0.5, 2.0);
	}
	double chi_squared = double(4*block_len)*sum;

	double arg1 = (double)N/2.0;
	double arg2 = chi_squared/2.0;
	nb_P_values = 1;
	double* P_value = new double[1];
	P_value[0] = gsl_sf_gamma_inc_Q(arg1,arg2);

	if (details)
	{
		output_stream << "Block frequency test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "Nb. of blocks: " << N << std::endl;
		output_stream << "Block length: " << block_len << std::endl;
		output_stream << "Nb. of discarded bits: " << nb_bits-N*block_len << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::CumulativeSums(void)
{
	int i,k,start,finish,mode;
	double cusum,z,sum,sum1,sum2;

	nb_P_values = 2;
	double* P_value = new double[2];
	for (mode = 0; mode < 2; mode++) // mode = {0,1}  => {forward,reverse}
	{
		sum = 0.0;
		cusum = 1;
		if (mode == 0)
			for (i = 0; i < nb_bits; i++)
			{
				sum += 2*(int)epsilon[i].b - 1;
				cusum = std::max(cusum, fabs(sum));
			}
		else if (mode == 1)
			for (i = nb_bits-1; i >= 0; i--)
			{
				sum += 2*(int)epsilon[i].b - 1;
				cusum = std::max(cusum, fabs(sum));
			}
		z = cusum;

		sum1 = 0.0;
		start = (-nb_bits/(int)z+1)/4;
		finish = (nb_bits/(int)z-1)/4;
		for (k=start;k<=finish;k++)
			sum1 += (gsl_cdf_ugaussian_P((4*k+1)*z/sqrt(nb_bits))-gsl_cdf_ugaussian_P((4*k-1)*z/sqrt(nb_bits)));

		sum2 = 0.0;
		start = (-nb_bits/(int)z-3)/4;
		finish = (nb_bits/(int)z-1)/4;
		for (k=start;k<=finish;k++)
			sum2 += (gsl_cdf_ugaussian_P((4*k+3)*z/sqrt(nb_bits))-gsl_cdf_ugaussian_P((4*k+1)*z/sqrt(nb_bits)));
		P_value[mode] = (1.0 - sum1 + sum2);
	}

	if (details)
	{
		output_stream << "Cumulative sums test\t" << test_result[P_value[0]>=alpha]
		                  << "\t" << test_result[P_value[1]>=alpha] << std::endl;
		output_stream << "forward P-value: " << P_value[0] << std::endl;
		output_stream << "reverse P-value: " << P_value[1] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::Runs(void)
{
	register int i;
	double argument, pi, V_n_obs, tau;
	double product;
	nb_P_values = 1;
	double* P_value = new double[1];

	int* r = new int[nb_bits];

	double sum = 0.0;
	for (i = 0; i < nb_bits; i++)
		sum += (int)epsilon[i].b;
	pi = sum/nb_bits;
	tau = 2.0/sqrt(nb_bits);
	if (fabs(pi - 0.5) < tau)
	{
		for (i = 0; i < nb_bits-1; i++)
		{
			if ((int)epsilon[i].b == (int)epsilon[i+1].b)
				r[i] = 0;
			else
				r[i] = 1;
		}
		V_n_obs = 0;
		for (i = 0; i < nb_bits-1; i++)
			V_n_obs += r[i];
		V_n_obs++;
		product = pi * (1.e0 - pi);
		argument = fabs(V_n_obs - 2.0*nb_bits*product)/(2.0*sqrt(2.0*nb_bits)*product);
		P_value[0] = gsl_sf_erfc(argument);
		delete[] r;
	}
	else
	{
		output_stream << "Runs: ESTIMATOR CRITERIA NOT MET! PI = " << pi << std::endl;
		delete[] r;
		nb_P_values = 0;
		return NULL;
	}

	if (details)
	{
		output_stream << "Runs test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::LongestRunOfOnes(const int &long_runs_case)
{
	double *pi=NULL;
	int K=0,M=0;
	if (long_runs_case==8)
	{
		pi = new double[4];
		pi[0] = 0.21484375;
		pi[1] = 0.3671875;
		pi[2] = 0.23046875;
		pi[3] = 0.1875;
		K = 3;
		M = 8;
	}
	else if (long_runs_case==128)
	{
		pi = new double[6];
		pi[0] = 0.1174035788;
		pi[1] = 0.242955959;
		pi[2] = 0.249363483;
		pi[3] = 0.17517706;
		pi[4] = 0.102701071;
		pi[5] = 0.112398847;
		K = 5;
		M = 128;
	}
	else if (long_runs_case==10000)
	{
		pi = new double[7];
		pi[0] = 0.0882;
		pi[1] = 0.2092;
		pi[2] = 0.2483;
		pi[3] = 0.1933;
		pi[4] = 0.1208;
		pi[5] = 0.0675;
		pi[6] = 0.0727;
		K = 6;
		M = 10000;
	}

	double run, v_n_obs, sum, chi2;
	nb_P_values = 1;
	double* P_value = new double[1];
	register int i,j;

	unsigned int *nu = new unsigned int[K+1];
	for (i=0;i<=K;i++)
		nu[i] = 0;

	int N = nb_bits/M;
	for (i = 0; i < N; i++)
	{
		v_n_obs = 0.e0;
		run = 0.e0;
		for (j = i*M; j < (i+1)*M; j++)
		{
			if ((int)epsilon[j].b == 1)
			{
				run++;
				v_n_obs = std::max(v_n_obs, run);
			}
			else
				run = 0.e0;
		}
		if (long_runs_case==8)
		{
			if ((int)v_n_obs == 0 || (int)v_n_obs == 1)
				nu[0]++;
			else if ((int)v_n_obs == 2)
				nu[1]++;
			else if ((int)v_n_obs == 3)
				nu[2]++;
			else if ((int)v_n_obs >= 4)
				nu[3]++;
		}
		else if (long_runs_case==128)
		{
			if ((int)v_n_obs <= 4)
				nu[0]++;
			else if ((int)v_n_obs == 5)
				nu[1]++;
			else if ((int)v_n_obs == 6)
				nu[2]++;
			else if ((int)v_n_obs == 7)
				nu[3]++;
			else if ((int)v_n_obs == 8)
				nu[4]++;
			else if ((int)v_n_obs >= 9)
				nu[5]++;
		}
		else if (long_runs_case==10000)
		{
			if ((int)v_n_obs <= 10)
				nu[0]++;
			else if ((int)v_n_obs == 11)
				nu[1]++;
			else if ((int)v_n_obs == 12)
				nu[2]++;
			else if ((int)v_n_obs == 13)
				nu[3]++;
			else if ((int)v_n_obs == 14)
				nu[4]++;
			else if ((int)v_n_obs == 15)
				nu[5]++;
			else if ((int)v_n_obs >= 16)
				nu[6]++;
		}
	}
	chi2 = 0.0;					// Compute Chi Square
	sum = 0;
	for (i = 0; i < K+1; i++)
	{
		chi2 += pow((double)nu[i] - (double)N*pi[i],2)/((double)N*pi[i]);
		sum += nu[i];
	}
	P_value[0] = gsl_sf_gamma_inc_Q(K/2.0, chi2/2.0);

	delete[] pi;
	delete[] nu;

	if (details)
	{
		output_stream << "LongestRunOfOnes test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::Rank(void)
{
	int N = nb_bits/1024;// NUMBER OF MATRICES
	if (N==0)
	{
		output_stream << "Rank: Insufficient # Of Bits To Define An 32x32 matrix\n";
		nb_P_values = 0;
		return NULL;
	}

	// COMPUTE PROBABILITIES
	register int        i, k;
	double product = 1;
	for (i = 0; i <= 31; i++)
		product *= ((1.0-pow(2.0, double(i-32)))*(1.0-pow(2.0, double(i-32))))/(1.0-pow(2.0, double(i-32)));
	double p_32 = product;

	product = 1;
	for (i = 0; i <= 30; i++)
		product *= ((1.0-pow(2.0, double(i-32)))*(1.0-pow(2, double(i-32))))/(1.0-pow(2, double(i-31)));
	double p_31 = product/2.0;

	double p_30 = 1.0 - (p_32+p_31);

	int R;
	int F_31 = 0, F_32 = 0;
	Bit* matrix = new Bit[1024];
	for (k = 0; k < N; k++)// FOR EACH 32x32 MATRIX
			{
		def_matrix(matrix, k);
		R = computeRank(matrix);
		if (R == 32) F_32++;			// DETERMINE FREQUENCIES
		if (R == 31) F_31++;
			}
	delete[] matrix;
	int F_30 = N - (F_32+F_31);

	double chi_squared =(pow(double(F_32) - double(N)*p_32, 2.0)/(double(N)*p_32)) +
			pow(double(F_31) - double(N)*p_31, 2.0)/(double(N)*p_31) +
			pow(double(F_30) - double(N)*p_30, 2.0)/(double(N)*p_30);

	nb_P_values = 1;
	double* P_value = new double[1];
	P_value[0] = exp(-chi_squared/2.0);

	if (details)
	{
		output_stream << "Rank test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "Number of matrices 32x32: " << N << std::endl;
		output_stream << "Number of matrices of full rank: " << F_32 << std::endl;
		output_stream << "Number of matrices of full rank - 1: " << F_31 << std::endl;
		output_stream << "Number of remaining matrices: " << F_30 << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

void Statistical_Tests::def_matrix(Bit *matrix, const int &k)
{
	register int i,j;
	for (i = 0; i < 32; i++)
		for (j = 0; j < 32; j++)
			matrix[i+32*j].b = epsilon[k*1024+j+i*32].b;
}

int Statistical_Tests::computeRank(Bit* matrix)
{
	register int i;

	// FORWARD APPLICATION OF ELEMENTARY ROW OPERATIONS

	for (i = 0; i < 31; i++)
		if (matrix[i+32*i].b == 1)
			perform_elementary_row_operations(matrix, 0, i);
		else
		{ 	// matrix[i][i] = 0
			if (find_unit_element_and_swap(matrix, 0, i) == 1)
				perform_elementary_row_operations(matrix, 0, i);
		}

	// BACKWARD APPLICATION OF ELEMENTARY ROW OPERATIONS
	for (i = 31; i > 0; i--)
		if (matrix[i+32*i].b == 1)
			perform_elementary_row_operations(matrix, 1, i);
		else
		{ 	// matrix[i][i] = 0
			if (find_unit_element_and_swap(matrix, 1, i) == 1)
				perform_elementary_row_operations(matrix, 1, i);
		}
	return determine_rank(matrix);
}

void Statistical_Tests::perform_elementary_row_operations(Bit* A, const int &flag, const int &i)
{
	register int j,k;

	switch (flag)
	{
	case 0:
		for (j = i+1; j < 32;  j++)
			if (A[j+32*i].b == 1)
				for (k = i; k < 32; k++)
					A[j+32*k].b = (A[j+32*k].b + A[i+32*k].b) % 2;
		break;

	case 1:
		for (j = i-1; j >= 0;  j--)
			if (A[j+32*i].b == 1)
				for (k = 0; k < 32; k++)
					A[j+32*k].b = (A[j+32*k].b + A[i+32*k].b) % 2;
		break;
	}
}

int Statistical_Tests::find_unit_element_and_swap(Bit* A, const int &flag, const int &i)
{
	register int index;
	int row_op = 0;

	switch (flag)
	{
	case 0:
		index = i+1;
		while ((index < 32) && (A[index+32*i].b == 0))
			index++;
		if (index < 32)
			row_op = swap_rows(A, i, index);
		break;
	case 1:
		index = i-1;
		while ((index >= 0) && (A[index+32*i].b == 0))
			index--;
		if (index >= 0)
			row_op = swap_rows(A, i, index);
		break;
	}
	return row_op;
}

int Statistical_Tests::swap_rows(Bit* A, const int &i, const int &index)
{
	register int p;
	Bit temp;
	for (p = 0; p < 32; p++)
	{
		temp.b = A[i+32*p].b;
		A[i+32*p].b = A[index+32*p].b;
		A[index+32*p].b = temp.b;
	}
	return 1;
}

int Statistical_Tests::determine_rank(const Bit* A)
{
	// DETERMINE RANK, THAT IS, COUNT THE NUMBER OF NONZERO ROWS
	int allZeroes;
	int rank = 32;
	register int i, j;
	for (i = 0; i < 32; i++)
	{
		allZeroes = 1;
		for (j=0; j < 32; j++)
		{
			if (A[i+32*j].b == 1)
			{
				allZeroes = 0;
				break;
			}
		}
		if (allZeroes == 1) rank--;
	}
	return rank;
}

double* Statistical_Tests::DiscreteFourierTransform(void)
{
	//generate data
	double *data = new double[nb_bits];
	register int i;
	for (i = 0; i < nb_bits; i++)
		data[i] = 2.0*(double)epsilon[i].b - 1.0;

	//compute FFT
	gsl_fft_real_wavetable *real;
	gsl_fft_real_workspace *work;
	real = gsl_fft_real_wavetable_alloc(nb_bits);
	work = gsl_fft_real_workspace_alloc(nb_bits);
	gsl_fft_real_transform(data, 1, nb_bits, real, work);
	gsl_fft_real_wavetable_free(real);//free memory
	gsl_fft_real_workspace_free(work);

	// COMPUTE MAGNITUDE
	int m_len;
	double *m = NULL;
	if ((nb_bits%2)==0)//nb_bits is even
	{
		m_len = nb_bits/2+1;
		m = new double[m_len];
		m[0] = fabs(data[0]);//sqrt(data.^2)
		for (i = 0; i < (m_len-2); i++)
			m[i+1] = sqrt(pow(data[2*i+1], 2.0)+pow(data[2*i+2], 2.0));
		m[nb_bits/2] = data[nb_bits-1];
	}
	else //nb_bits is odd
	{
		m_len = (nb_bits-1)/2+1;
		m =  new double[m_len];
		m[0] = fabs(data[0]);//sqrt(data.^2)
		for (i=0;i<(m_len-1);i++)
			m[i+1] = sqrt(pow(data[2*i+1], 2.0)+pow(data[2*i+2], 2.0));
	}
	delete[] data;//free memory

	// CONFIDENCE INTERVAL
	int count = 0;
	double upperBound = sqrt(3*nb_bits);
	for (i = 0; i < m_len; i++)
		if (m[i] < upperBound) count++;
	delete[] m;//free memory

	double N_l = double(count);       // number of peaks less than h = sqrt(3*n)
	double N_o = 0.95*double(m_len); // theoretical value
	double d = (N_l - N_o)/sqrt(m_len*0.95*0.05);

	nb_P_values = 1;
	double* P_value = new double[1];
	P_value[0] = gsl_sf_erfc(fabs(d)/sqrt(2.0));

	if (details)
	{
		output_stream << "DiscreteFourierTransform test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "Upper bound: " << upperBound << std::endl;
		output_stream << "Number of peaks less than upper bound: " << count << std::endl;
		output_stream << "Theoretical number of peaks less than upper bound: " << N_o << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::NonOverlappingTemplateMatchings(const int &template_len)
{
	int workspace_len = 8*sizeof(unsigned int);
	if (workspace_len<template_len)
	{
		output_stream << "NonOverlappingTemplateMatchings: template length must be at most "
				<< workspace_len << std::endl;
		nb_P_values = 0;
		return NULL;
	}

	int N = 8;//number of independent blocks
	int M = nb_bits/N;//block length

	double lambda = double(M-template_len+1)/pow(2.0, double(template_len));
	if (lambda<=0)
	{
		output_stream << "NonOverlappingTemplateMatchings: Lambda (" << lambda << ") is not positive!\n";
		nb_P_values = 0;
		return NULL;
	}

	unsigned int   W_obs;
	double         chi2;
	int            i, j, k, match;

	unsigned int *Wj = new unsigned int[N];
	double varWj = double(M)*(1.0/pow(2.0, double(template_len)) -
			double(2*template_len-1)/pow(2.0, double(2*template_len)));

	long long int num = (long long int)1 << template_len;
	Bit* pattern = new Bit[template_len];
	unsigned int *workspace = new unsigned int[workspace_len];
	double* P_value = new double[num];
	nb_P_values = 0;
	int nb_passes = 0;
	long long int jj;
	for (jj = 0; jj < num; jj++)
	{
		//generate template
		if (generate_nonperiodic_template(jj, pattern, template_len, workspace, workspace_len))//if template is found
		{
			//compare template with substring
			for (i = 0; i < N; i++)
			{
				W_obs = 0;
				for (j = 0; j < M-template_len+1; j++)
				{
					match = 1;
					for (k = 0; k < template_len; k++)
					{
						if ((int)pattern[k].b != (int)epsilon[i*M+j+k].b)
						{
							match = 0;
							break;
						}
					}
					if (match == 1) W_obs++;
				}
				Wj[i] = W_obs;
			}
			chi2 = 0.0;                                   // Compute Chi Square
			for (i = 0; i < N; i++)
				chi2 += pow(((double)Wj[i] - lambda)/pow(varWj,0.5),2);
			P_value[nb_P_values] = gsl_sf_gamma_inc_Q(N/2.0, chi2/2.0);
			if (P_value[nb_P_values]>=alpha) nb_passes++;
			nb_P_values++;
		}
	}
	//free memory
	delete[] Wj;
	delete[] pattern;
	delete[] workspace;

	if (details)
	{
		output_stream << "NonOverlappingTemplateMatchings test\t" << 100*nb_passes/nb_P_values;
		output_stream << "% " << test_result[1] << "\t" << 100-100*nb_passes/nb_P_values;
		output_stream << "% " << test_result[0] << std::endl;
		output_stream << "Number of independent blocks: " << N << std::endl;
		output_stream << "Block length: " << M << std::endl;
		output_stream << "Number of templates: " << nb_P_values << std::endl;
		for (i = 0; i < nb_P_values; i++)
			output_stream << "#" << i << " P-value: " << P_value[i] << std::endl;
	}

	//return vector of p_values for each template
	return P_value;
}

bool Statistical_Tests::generate_nonperiodic_template(const int &index,
		Bit* pattern, const int &template_len, unsigned int* workspace, const int &len)
{
	int displayMask = 1 << (len-1);
	int value = index;
	register int i, match, c;
	for (c = 0; c < len; c++)
	{
		workspace[c] = (value & displayMask)?1:0;
		value <<= 1;
	}
	for (i = 1; i < template_len; i++)
	{
		match = 1;
		if ((workspace[len-template_len] !=  workspace[len-1]) &&
				((workspace[len-template_len] !=  workspace[len-2])||
						(workspace[len-template_len+1] != workspace[len-1])))
		{
			for (c = len-template_len; c <= (len-1)-i; c++)
			{
				if (workspace[c] != workspace[c+i])
				{
					match = 0;
					break;
				}
			}
		}
		if (match==1)
			return false;//value does not generate a nonperiodic pattern
	}
	for (c = len-template_len; c <= (len-1); c++)
		pattern[(c+template_len-len)].b = workspace[c];
	return true;//found a nonperiodic pattern
}

double* Statistical_Tests::OverlappingTemplateMatchings(const int &template_len)
{
	int    i, k, match;
	double W_obs, eta, sum, chi2;
	double lambda;
	nb_P_values = 1;
	double* P_value = new double[1];
	int    M = 1032, N, j, K = 5;
	unsigned int nu[6] = {0,0,0,0,0,0};
	double pi[6] = {0.143783,0.139430,0.137319,0.124314,0.106209,0.348945};

	N = nb_bits/M;

	Bit* pattern = new Bit[template_len];
	for (i = 0; i < template_len; i++)
		pattern[i].b = 1;

	lambda = (double)(M-template_len+1)/pow(2.0,double(template_len));
	eta = lambda/2.0;
	sum = 0.0;
	for (i = 0; i < K; i++)
	{			// Compute Probabilities
		pi[i] = Pr(i,eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;

	for (i = 0; i < N; i++)
	{
		W_obs = 0;
		for (j = 0; j < M-template_len+1; j++)
		{
			match = 1;
			for (k = 0; k < template_len; k++)
			{
				if ((int)pattern[k].b != (int)epsilon[i*M+j+k].b)
					match = 0;
			}
			if (match == 1) W_obs++;
		}
		if (W_obs <= 4)
			nu[(int)W_obs]++;
		else
			nu[K]++;
	}
	delete[] pattern;

	chi2 = 0.0;                                   // Compute Chi Square
	for (i = 0; i < K+1; i++)
		chi2 += pow((double)nu[i] - (double)N*pi[i],2.0)/((double)N*pi[i]);
	P_value[0] = gsl_sf_gamma_inc_Q(double(K)/2.0, chi2/2.0);

	if (details)
	{
		output_stream << "OverlappingTemplateMatchings test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "Number of independent blocks: " << N << std::endl;
		output_stream << "Block length: " << M << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

double Statistical_Tests::Pr(const int &u, const double &eta)
{
	int    l;
	double sum, p;

	if (u == 0)
		p = exp(-eta);
	else
	{
		sum = 0.0;
		for (l = 1; l <= u; l++)
			sum += exp(-eta-u*log(2)+l*log(eta)-gsl_sf_lngamma(double(l+1))+
					gsl_sf_lngamma(double(u))-gsl_sf_lngamma(double(l))-gsl_sf_lngamma(double(u-l+1)));
		p = sum;
	}
	return p;
}

double* Statistical_Tests::Universal(void)
{
	int block_len, nb_ini_steps;
	int      i, j, p, K;
	double   arg, sigma, phi, sum, c;
	int*    T, decRep;
	double   expected_value[17] = {0,0,0,0,0,0,5.2177052,6.1962507,7.1836656,
			8.1764248, 9.1723243, 10.170032, 11.168765,
			12.168070, 13.167693, 14.167488, 15.167379};
	double   variance[17] = {0,0,0,0,0,0,2.954,3.125,3.238,3.311,3.356,3.384,
			3.401,3.410,3.416,3.419,3.421};

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * THE FOLLOWING REDEFINES L, SHOULD THE CONDITION:     n >= 1010*2^L*L       *
	 * NOT BE MET, FOR THE BLOCK LENGTH L.                                        *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	if (nb_bits >= 387840)     block_len = 6;
	if (nb_bits >= 904960)     block_len = 7;
	if (nb_bits >= 2068480)    block_len = 8;
	if (nb_bits >= 4654080)    block_len = 9;
	if (nb_bits >= 10342400)   block_len = 10;
	if (nb_bits >= 22753280)   block_len = 11;
	if (nb_bits >= 49643520)   block_len = 12;
	if (nb_bits >= 107560960)  block_len = 13;
	if (nb_bits >= 231669760)  block_len = 14;
	if (nb_bits >= 496435200)  block_len = 15;
	if (nb_bits >= 1059061760) block_len = 16;

	nb_ini_steps = 10*(1 << block_len);
	K = nb_bits/block_len-nb_ini_steps;	 		    // BLOCKS TO TEST

	if ((block_len < 6) || (block_len > 16) || ((double)nb_ini_steps < 10*pow(2,block_len)))
	{
		output_stream << "Universal: block_len IS OUT OF RANGE.\n";
		output_stream << "Universal: nb_ini_steps IS LESS THAN " << 10*(1 << block_len) << std::endl;
		return NULL;
	}
	// COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper

	c = 0.7 - 0.8/double(block_len) + (4.0 + 32.0/double(block_len))*pow(double(K), -3.0/double(block_len))/15.0;
	sigma = c * sqrt(variance[block_len]/(double)K);
	sum = 0.0;
	p = (1 << block_len);
	T = new int[p];
	for (i = 0; i < p; i++) T[i] = 0;
	for (i = 1; i <= nb_ini_steps; i++)
	{		// INITIALIZE TABLE
		decRep = 0;
		for (j = 0; j < block_len; j++) decRep += epsilon[(i-1)*block_len+j].b * (1 << (block_len-1-j));
		T[decRep] = i;
	}
	for (i = nb_ini_steps+1; i <= nb_ini_steps+K; i++)
	{ 	// PROCESS BLOCKS
		decRep = 0;
		for (j = 0; j < block_len; j++) decRep += epsilon[(i-1)*block_len+j].b * (1 << (block_len-1-j));
		sum += log(i - T[decRep])/log(2);
		T[decRep] = i;
	}
	delete[] T;//free memory
	phi = sum/double(K);

	arg = fabs(phi-expected_value[block_len])/(sqrt(2.0) * sigma);

	nb_P_values = 1;
	double* P_value = new double[1];
	P_value[0] = gsl_sf_erfc(arg);

	if (details)
	{
		output_stream << "Universal test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "Block length: " << block_len << std::endl;
		output_stream << "Number of initialisation steps: " << nb_ini_steps << std::endl;
		output_stream << "Number of blocks to test: " << K << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	//return p_value
	return P_value;
}

double* Statistical_Tests::ApproximateEntropy(const int &first_block_len)
{
	if (first_block_len > (int)(log(nb_bits)/log(2.0)-5))
	{
		output_stream << "ApproximateEntropy: first_block_len = " << first_block_len
				<< " exceeds recommended values = " << std::max(1,(int)(log(nb_bits)/log(2)-5)) << std::endl;
		nb_P_values = 0;
		return NULL;
	}

	int           r, blockSize;
	double        ApEn[2];

	r = 0;
	for (blockSize = first_block_len; blockSize <= first_block_len+1; blockSize++)
	{
		if (blockSize == 0)
		{
			ApEn[0] = 0.00;
			r++;
		}
		else
		{
			int powLen = (1 << (blockSize+1))-1;
			unsigned int* P = new unsigned int[powLen];
			register int i,j,k;
			for (i = 0; i < powLen; i++) P[i] = 0;
			for (i = 0; i < nb_bits; i++)
			{ /* COMPUTE FREQUENCY */
				k = 1;
			for (j = 0; j < blockSize; j++)
			{
				if ((int)epsilon[(i+j)%nb_bits].b == 0)
					k *= 2;
				else if ((int)epsilon[(i+j)%nb_bits].b == 1)
					k = 2*k+1;
			}
			P[k-1]++;
			}
			double sum = 0.0;
			int len = (1 << blockSize);
			int index = len-1;
			for (i = 0; i < len; i++)
			{
				if (P[index] > 0) sum += double(P[index])*log(double(P[index])/double(nb_bits));
				index++;
			}
			ApEn[r] = sum/double(nb_bits);
			r++;
			delete[] P;
		}
	}
	double apen = ApEn[0] - ApEn[1];

	double chi_squared = 2.0*double(nb_bits)*(log(2.0) - apen);

	nb_P_values = 1;
	double* P_value = new double[1];
	P_value[0] = gsl_sf_gamma_inc_Q(double(1 << (first_block_len-1)), chi_squared/2.0);

	if (details)
	{
		output_stream << "ApproximateEntropy test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "First block length: " << first_block_len << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::RandomExcursions(void)
{
	int    b, i, j, k, J, x, constraint;
	int    stateX[8] = {-4,-3,-2,-1,1,2,3,4};
	int    counter[8] = {0,0,0,0,0,0,0,0};
	double sum, nu[6][8];
	double pi[10][6] = {{0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000},
			{.5000000000,.2500000000,.1250000000,.06250000000,.03125000000,.0312500000},
			{.7500000000,.06250000000,.04687500000,.03515625000,.02636718750,.0791015625},
			{.8333333333,.02777777778,.02314814815,.01929012346,.01607510288,.0803755143},
			{.8750000000,.01562500000,.01367187500,.01196289063,.01046752930,.0732727051}};

	int *S_k = new int[nb_bits];
	int cycle_len = std::max(1000, nb_bits/128);
	int *cycle = new int[cycle_len];

	J = 0; 					/* DETERMINE CYCLES */
	S_k[0] = 2*(int)epsilon[0].b - 1;
	for (i = 1; i < nb_bits; i++)
	{
		S_k[i] = S_k[i-1] + 2*(int)epsilon[i].b - 1;
		if (S_k[i] == 0)
		{
			J++;
			if (J >= cycle_len)
			{
				output_stream << "RandomExcursions: exceeding the max" << cycle_len << std::endl;
				delete[] S_k;
				delete[] cycle;
				nb_P_values = 0;
				return NULL;
			}
			cycle[J] = i;//cycle stop
		}
	}
	if (S_k[nb_bits-1] != 0) // unfinished cycle
	{
		J++;
		cycle[J] = nb_bits-1;
	}

	nb_P_values = 8;
	double *P_value = new double[8];
	int nb_passes = 0;
	constraint = std::max(int(0.005*sqrt(nb_bits)), 500);
	if (J < constraint)
	{
		output_stream << "RandomExcursions: test not applicable. Insufficient number of cycles\n";
		output_stream << "Number of cycles " << J << " less than rejection constraint "
				<< constraint << "." << std::endl;
		delete[] S_k;
		delete[] cycle;
		delete[] P_value;
		nb_P_values = 0;
		return NULL;
	}
	else
	{
		for (k = 0; k < 6; k++)
			for (i = 0; i < 8; i++)
				nu[k][i] = 0.0;
		int cycleStart = 0;
		int cycleStop  = cycle[1];
		for (j = 1; j <= J; j++)
		{                           // FOR EACH CYCLE
			for (i = 0; i < 8; i++) counter[i] = 0;
			for (i = cycleStart; i < cycleStop; i++)
			{
				if ((S_k[i] >= 1 && S_k[i] <= 4)||(S_k[i] >= -4 && S_k[i] <= -1))
				{
					if (S_k[i] < 0)
						b = 4;
					else
						b = 3;
					counter[S_k[i]+b]++;
				}
			}
			for (i = 0; i < 8; i++)
			{
				if (counter[i] >= 0 && counter[i] <= 4)
					nu[counter[i]][i]++;
				else if (counter[i] >= 5)
					nu[5][i]++;
			}
			if (j==J) break;//avoid to read from uninitialized location
			cycleStart = cycle[j]+1;
			cycleStop = cycle[j+1];
		}
		for (i = 0; i < 8; i++)
		{
			x = stateX[i];
			sum = 0.0;
			for (k = 0; k < 6; k++)
				sum += pow(nu[k][i] - double(J)*pi[(int)fabs(x)][k], 2.0)/double(J*pi[(int)fabs(x)][k]);
			P_value[i] = gsl_sf_gamma_inc_Q(2.5, sum/2.0);
			if (P_value[i]>=alpha) nb_passes++;
		}
	}
	delete[] S_k;
	delete[] cycle;

	if (details)
	{
		output_stream << "RandomExcursions test\t" << 100*nb_passes/nb_P_values;
		output_stream << "% " << test_result[1] << "\t" << 100-100*nb_passes/nb_P_values;
		output_stream << "% " << test_result[0] << std::endl;
		output_stream << "Number of cycles: " << J << std::endl;
		output_stream << "Rejection constraint " << constraint << std::endl;
		for (i = 0; i < 8; i++)
			output_stream << "#" << i << " P-value: " << P_value[i] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::RandomExcursionsVariant(void)
{
	int    i, p, J, x;
	int    stateX[18] = {-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9};
	int    count;
	int*   S_k = new int[nb_bits];

	J = 0;
	S_k[0] = 2*(int)epsilon[0].b - 1;
	for (i = 1; i < nb_bits; i++)
	{
		S_k[i] = S_k[i-1] + 2*(int)epsilon[i].b - 1;
		if (S_k[i] == 0) J++;
	}
	if (S_k[nb_bits-1] != 0) J++;

	int constraint = std::max(int(0.005*sqrt(nb_bits)), 500);
	if (J < constraint)
	{
		output_stream << "RandomExcursionsVariant: test not applicable. Insufficient number of cycles\n";
		output_stream << "Number of cycles " << J << " less than rejection constraint "
				<< constraint << "." << std::endl;
		delete[] S_k;
		nb_P_values = 0;
		return NULL;
	}

	nb_P_values = 18;
	int nb_passes = 0;
	double *P_value = new double[18];
	for (p = 0; p <= 17; p++)
	{
		x = stateX[p];
		count = 0;
		for (i = 0; i < nb_bits; i++)
			if (S_k[i] == x) count++;
		P_value[p] = gsl_sf_erfc(fabs(count-J)/(sqrt(2.0*J*(4.0*fabs(x)-2))));
		nb_passes++;
	}
	delete[] S_k;

	if (details)
	{
		output_stream << "RandomExcursionsVariant test\t" << 100*nb_passes/nb_P_values;
		output_stream << "% " << test_result[1] << "\t" << 100-100*nb_passes/nb_P_values;
		output_stream << "% " << test_result[0] << std::endl;
		output_stream << "Number of cycles: " << J << std::endl;
		output_stream << "Rejection constraint " << constraint << std::endl;
		for (p=0;p<nb_P_values;p++)
			output_stream << "#" << p << " P-value: " << P_value[p] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::Serial(const int &pattern_len)
{
	double psim0 = psi2(pattern_len);
	double psim1 = psi2(pattern_len-1);
	double psim2 = psi2(pattern_len-2);
	double del1 = psim0 - psim1;
	double del2 = psim0 - 2.0*psim1 + psim2;
	nb_P_values = 2;
	double *P_value = new double[2];
	int pow2 = (1 << (pattern_len-2));
	P_value[0] = gsl_sf_gamma_inc_Q(double(pow2), del1/2.0);
	P_value[1] = gsl_sf_gamma_inc_Q(double(pow2/2), del2/2.0);

	if (details)
	{
		output_stream << "Serial test\t" << test_result[P_value[0]>=alpha] << "\t"
				<< test_result[P_value[1]>=alpha] << std::endl;
		for (int i=0;i<2;i++)
			output_stream << "#" << i << " P-value: " << P_value[i] << std::endl;
	}

	return P_value;
}

double Statistical_Tests::psi2(const int &m)
{
	if ((m == 0) || (m == -1)) return 0.0;

	int powLen = (1 << (m+1))-1;
	int* P = new int[powLen];

	register int i,j,k;
	for (i = 0; i < powLen; i++) P[i] = 0;	  /* INITIALIZE NODES */
	for (i = 0; i < nb_bits; i++)
	{		 /* COMPUTE FREQUENCY */
		k = 1;
		for (j = 0; j < m; j++)
		{
			if (epsilon[(i+j)%nb_bits].b == 0)
				k *= 2;
			else if (epsilon[(i+j)%nb_bits].b == 1)
				k = 2*k+1;
		}
		P[k-1]++;
	}

	double sum = 0.0;
	powLen = (1 << m);
	for (i = powLen-1; i < 2*powLen-1; i++)
		sum += pow(double(P[i]), 2.0);
	delete[] P;

	return (sum * double(powLen)/double(nb_bits)) - double(nb_bits);
}

double* Statistical_Tests::LempelZivCompression(void)
{
	int internal_nb_bits=0;
	double mean=0.0,variance=0.0;
	if (nb_bits>=1000000)
	{
		internal_nb_bits = 1000000;
		mean = 69588.20190000;
		variance = 73.23726011;
	}
	else if (nb_bits>=800000)
	{
		internal_nb_bits = 800000;
		mean = 56821.9500;
		variance = 67.4184;
	}
	else if (nb_bits>=600000)
	{
		internal_nb_bits = 600000;
		mean = 43787.5000;
		variance = 70.1579;
	}
	else if (nb_bits>=400000)
	{
		internal_nb_bits = 400000;
		mean = 30361.9500;
		variance = 58.7868;
	}
	else if (nb_bits>=200000)
	{
		internal_nb_bits = 200000;
		mean = 16292.1000;
		variance = 21.4632;
	}
	else if (nb_bits>=100000)
	{
		internal_nb_bits = 100000;
		mean = 8782.500000;
		variance = 11.589744;
	}
	else
	{
		output_stream << "LempelZivCompression: number of bits must be at least 100000. Abording test\n";
		nb_P_values = 0;
		return NULL;
	}
	if (nb_bits != internal_nb_bits)
	{
		output_stream << "LempelZivCompression: no precomputed values available for mean and variance for "
				<< nb_bits << " bits.\n";
		output_stream << "Choosing " << internal_nb_bits << " bits.\n";
	}	

	int k = (int)(log(internal_nb_bits)/log(2)+6);
	int powLen = (1 << k);
	Bit* P = new Bit[powLen];
	register int i,j;
	for (i = 0; i < powLen; i++)
		P[i].b = 0;

	int W = 0;//Number of words
	int prev_I;
	i = 0;
	int max_len = 0;
	int done;
	while (i <= internal_nb_bits-1)
	{
		done = 0;
		j = 1;
		prev_I = i;
		while (!done)
		{
			if (2*j+1 <= powLen)
			{
				if ((int)epsilon[i].b == 0)
				{
					if (P[2*j].b == 1)
					{
						j *= 2;
					}
					else
					{
						P[2*j].b = 1;
						done = 1;
					}
				}
				else
				{
					if (P[2*j+1].b == 1)
					{
						j = 2*j+1;
					}
					else
					{
						P[2*j+1].b = 1;
						done = 1;
					}
				}
				i++;
				if (i > internal_nb_bits-1)
				{
					done = 1;
				}
				if (done)
				{
					W++;
					max_len = std::max(max_len,i-prev_I);
				}
			}
			else
			{
				output_stream << "LempelZivCompression: segmentation violation imminent. Test is terminated\n";
				done = 1;
				i = internal_nb_bits;
			}
		}
	}
	delete[] P;

	nb_P_values = 1;
	double* P_value = new double[1];
	P_value[0] = 0.5*gsl_sf_erfc((mean-W)/sqrt(2.0*variance));

	if (details)
	{
		output_stream << "LempelZivCompression test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "Number of words: " << W << std::endl;
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

double* Statistical_Tests::LinearComplexity(const int &block_len)
{
	int       i, ii, j, d, N;
	int       L, m, N_, parity, sign;
	double    T_, mean;
	int       K = 6;
	double    pi[7]={0.01047,0.03125,0.12500,0.50000,0.25000,0.06250,0.020833};
	double    nu[7];

	N = nb_bits/block_len;
	Bit *B_ = new Bit[block_len];
	Bit *C = new Bit[block_len];
	Bit *P = new Bit[block_len];
	Bit *T = new Bit[block_len];

	for (i = 0; i < K+1; i++) nu[i] = 0.00;
	for (ii = 0; ii < N; ii++)
	{
		for (i = 0; i < block_len; i++)
		{
			B_[i].b = 0;
			C[i].b = 0;
			T[i].b = 0;
			P[i].b = 0;
		}
		L = 0;
		m = -1;
		d = 0;
		C[0].b = 1;
		B_[0].b = 1;
		// DETERMINE LINEAR COMPLEXITY
		N_ = 0;
		while (N_ < block_len)
		{
			d = (int)epsilon[ii*block_len+N_].b;
			for (i = 1; i <= L; i++)
				d += (int)C[i].b*(int)epsilon[ii*block_len+N_-i].b;
			d = d%2;
			if (d == 1)
			{
				for (i = 0; i < block_len; i++)
				{
					T[i].b = C[i].b;
					P[i].b = 0;
				}
				for (j = 0; j < block_len; j++)
					if (B_[j].b == 1) P[j+N_-m].b = 1;
				for (i = 0; i < block_len; i++)
					C[i].b = (C[i].b + P[i].b)%2;
				if (L <= N_/2)
				{
					L = N_ + 1 - L;
					m = N_;
					for (i = 0; i < block_len; i++)
						B_[i].b = T[i].b;
				}
			}
			N_++;
		}
		if ((parity = (block_len+1)%2) == 0)
			sign = -1;
		else
			sign = 1;
		mean = block_len/2. + (9.+sign)/36. - 1./pow(2,block_len) * (block_len/3. + 2./9.);
		if ((parity = block_len%2) == 0)
			sign = 1;
		else
			sign = -1;
		T_ = sign * (L - mean) + 2./9.;

		if (T_ <= -2.5)
			nu[0]++;
		else if (T_ > -2.5 && T_ <= -1.5)
			nu[1]++;
		else if (T_ > -1.5 && T_ <= -0.5)
			nu[2]++;
		else if (T_ > -0.5 && T_ <= 0.5)
			nu[3]++;
		else if (T_ > 0.5 && T_ <= 1.5)
			nu[4]++;
		else if (T_ > 1.5 && T_ <= 2.5)
			nu[5]++;
		else
			nu[6]++;
	}
	delete[] B_;
	delete[] P;
	delete[] C;
	delete[] T;

	double chi2 = 0.00;
	for (i = 0; i < K+1; i++)
		chi2 += pow(nu[i]-N*pi[i],2)/(N*pi[i]);

	nb_P_values = 1;
	double* P_value = new double[1];
	P_value[0] = gsl_sf_gamma_inc_Q(K/2.0,chi2/2.0);

	if (details)
	{
		output_stream << "LinearComplexity test\t" << test_result[P_value[0]>=alpha] << std::endl;
		output_stream << "Block length: " << block_len << std::endl;
		output_stream << "Number of blocks: " << N << std::endl;
		output_stream << "Discarded bits: " << nb_bits-N*block_len << std::endl; 
		output_stream << "P-value: " << P_value[0] << std::endl;
	}

	return P_value;
}

#endif
