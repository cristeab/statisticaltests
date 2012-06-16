#include "StatisticalTests.h"

//start reading beginning with the LSB
void StatisticalTests::setBinaryVector(int nb_offset_bytes)
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

double* StatisticalTests::frequency(void)
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

double* StatisticalTests::blockFrequency(int block_len)
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

double* StatisticalTests::cumulativeSums(void)
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

double* StatisticalTests::runs(void)
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

double* StatisticalTests::longestRunOfOnes(int long_runs_case)
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

double* StatisticalTests::rank(void)
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
		defMatrix(matrix, k);
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

void StatisticalTests::defMatrix(Bit *matrix, int k)
{
	register int i,j;
	for (i = 0; i < 32; i++)
		for (j = 0; j < 32; j++)
			matrix[i+32*j].b = epsilon[k*1024+j+i*32].b;
}

int StatisticalTests::computeRank(Bit* matrix)
{
	register int i;

	// FORWARD APPLICATION OF ELEMENTARY ROW OPERATIONS

	for (i = 0; i < 31; i++)
		if (matrix[i+32*i].b == 1)
			performElementaryRowOperations(matrix, 0, i);
		else
		{ 	// matrix[i][i] = 0
			if (findUnitElementAndSwap(matrix, 0, i) == 1)
				performElementaryRowOperations(matrix, 0, i);
		}

	// BACKWARD APPLICATION OF ELEMENTARY ROW OPERATIONS
	for (i = 31; i > 0; i--)
		if (matrix[i+32*i].b == 1)
			performElementaryRowOperations(matrix, 1, i);
		else
		{ 	// matrix[i][i] = 0
			if (findUnitElementAndSwap(matrix, 1, i) == 1)
				performElementaryRowOperations(matrix, 1, i);
		}
	return determineRank(matrix);
}

void StatisticalTests::performElementaryRowOperations(Bit* A, int flag, int i)
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

int StatisticalTests::findUnitElementAndSwap(Bit* A, int flag, int i)
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
			row_op = swapRows(A, i, index);
		break;
	case 1:
		index = i-1;
		while ((index >= 0) && (A[index+32*i].b == 0))
			index--;
		if (index >= 0)
			row_op = swapRows(A, i, index);
		break;
	}
	return row_op;
}

int StatisticalTests::swapRows(Bit* A, int i, int index)
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

int StatisticalTests::determineRank(const Bit* A)
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

double* StatisticalTests::discreteFourierTransform(void)
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

double* StatisticalTests::nonOverlappingTemplateMatchings(int template_len)
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
		if (generateNonperiodicTemplate(jj, pattern, template_len, workspace, workspace_len))//if template is found
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

bool StatisticalTests::generateNonperiodicTemplate(int index, Bit* pattern, int template_len, 
        unsigned int* workspace, int len)
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

double* StatisticalTests::overlappingTemplateMatchings(int template_len)
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

double StatisticalTests::Pr(int u, double eta)
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

double* StatisticalTests::universal(void)
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

double* StatisticalTests::approximateEntropy(int first_block_len)
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

double* StatisticalTests::randomExcursions(void)
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

double* StatisticalTests::randomExcursionsVariant(void)
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

double* StatisticalTests::serial(int pattern_len)
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

double StatisticalTests::psi2(int m)
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

double* StatisticalTests::lempelZivCompression(void)
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

double* StatisticalTests::linearComplexity(int block_len)
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

