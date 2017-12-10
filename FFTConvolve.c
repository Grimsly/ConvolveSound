/**
**FFT algortihm credited to Leonard Manzara's test.c
**/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

char ChunkID[4];
int ChunkSize;
char RIFFType[4];

char SubChunk1ID[4];
int SubChunk1Size;

short AudioFormat;
short NumChannels;

int SampleRate;
int ByteRate;

short BlockAlign;
short BitsPerSample;

char SubChunk2ID[4];
int SubChunk2Size;

int bytes;
int sampleNum;

short* data;

short* dataIR;
int sampleNumIR;

double* newData;
float maxSample = -1;
float MAX_VAL = 32767.f;

#define PI 3.141592653589793
#define TWO_PI (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

void four1(double data[], int nn, int isign){
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	n = nn << 1;
	j = 1;

	for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j + 1], data[i + 1]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax = 2;
	while (n > mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j + 1];
				tempi = wr * data[j + 1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j + 1] = data[i + 1] - tempi;
				data[i] += tempr;
				data[i + 1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}

int inWave(char* file)
{
	FILE* input = fopen(file, "rb");
	if (input != NULL)
	{
		printf("Reading %s.\n", file);
		fread(ChunkID, 1, 4, input);
		fread(&ChunkSize, 1, 4, input);
		fread(RIFFType, 1, 4, input);

		fread(SubChunk1ID, 1, 4, input);
		fread(&SubChunk1Size, 1, 4, input);
		fread(&AudioFormat, 1, 2, input);
		fread(&NumChannels, 1, 2, input);
		fread(&SampleRate, 1, 4, input);
		fread(&ByteRate, 1, 4, input);
		fread(&BlockAlign, 1, 2, input);
		fread(&BitsPerSample, 1, 2, input);

		if (SubChunk1Size == 18)	//for the extra bytes
		{
			short empty;
			fread(&empty, 1, 2, input);
		}

		fread(SubChunk2ID, 1, 4, input);
		fread(&SubChunk2Size, 1, 4, input);

		bytes = BitsPerSample / 8;
		sampleNum = SubChunk2Size / bytes;
		data = (short*)malloc(sizeof(short) * sampleNum);

		int index = 0;
		short sample = 0;
		while (fread(&sample, 1, bytes, input) == bytes)
		{
			data[index++] = sample;
			sample = 0;
		}

		fclose(input);
		printf("Closing %s.\n", file);
	}
	else
	{
		printf("Couldn't open specified file\n");
		return 0;
	}
	return 1;
}

int inIR(char* file){
	FILE* input = fopen(file, "rb");
	int subChunk2SizeIR;

	if (input != NULL){
		printf("Reading %s.\n", file);

		int IRsubChunk1Size;
		fseek(input, 16, SEEK_SET);
		fread(&IRsubChunk1Size, 1, 4, input);

		if(IRsubChunk1Size == 18){
			short empty;
			fread(&empty, 1, 2, input);
		}

		short BitsPerSampleIR;
		fseek(input, 34, SEEK_SET);
		fread(&BitsPerSampleIR, 1, 2, input);

		if (IRsubChunk1Size == 18){
			fseek(input, 42, SEEK_SET);
			fread(&subChunk2SizeIR, 1, 4, input);
			fseek(input, 46, SEEK_SET);
		}
		else{
			fseek(input, 40, SEEK_SET);
			fread(&subChunk2SizeIR, 1, 4, input);
			fseek(input, 44, SEEK_SET);
		}

		int BytesPerSampleIR = BitsPerSampleIR / 8;
		sampleNumIR = subChunk2SizeIR / BytesPerSampleIR;

		dataIR = (short*) malloc(sizeof(short) * sampleNumIR);
		int i = 0;
		short sample = 0;
		while (fread(&sample, 1, BytesPerSampleIR, input) == BytesPerSampleIR){
			dataIR[i++] = sample;
			sample = 0;
		}
        fclose(input);
		printf("Closing %s.\n", file);
	}
	else{
		printf("File not available.\n");
		return 0;
	}
	return 1;
}

int outWave(char* file, double* newInput, int newSize)
{
	FILE* outFile = fopen(file, "wb");
	if (outFile != NULL)
	{
		printf("\n Now Writing %s.\n", file);

		fwrite(ChunkID, 1, 4, outFile);
		fwrite(&ChunkSize, 1, 4, outFile);
		fwrite(RIFFType, 1, 4, outFile);

		fwrite(SubChunk1ID, 1, 4, outFile);
		fwrite(&SubChunk1Size, 1, 4, outFile);
		fwrite(&AudioFormat, 1, 2, outFile);
		fwrite(&NumChannels, 1, 2, outFile);
		fwrite(&SampleRate, 1, 4, outFile);
		fwrite(&ByteRate, 1, 4, outFile);
		fwrite(&BlockAlign, 1, 2, outFile);
		fwrite(&BitsPerSample, 1, 2, outFile);

		if (SubChunk1Size == 18)
		{
			short empty = 0;
			fwrite(&empty, 1, 2, outFile);
		}

		fwrite(SubChunk2ID, 1, 4, outFile);
		fwrite(&SubChunk2Size, 1, 4, outFile);
        
        for (int i = 0; i < newSize; i++)
		{
			if (i == 0)
			{
				maxSample = newInput[0];
			}
			else if (newInput[i] > maxSample)
			{
				maxSample = newInput[i];
			}
		}

		for (int i = 0; i < newSize; ++i)
		{
			newInput[i] = (newInput[i] / maxSample);
			short sample = (short)(newInput[i] * MAX_VAL);
			fwrite(&sample, 1, bytes, outFile);
		}

		free(newData);
		fclose(outFile);
		printf("Now Closing %s.\n", file);
	}
	else
	{
		printf("Couldn't open file\n");
		return 0;
	}
	return 1;
}

int main(int argc, char* argv[])
{
	clock_t begin, end;
	begin = clock();

	if (argc != 4)
	{
		printf("Need 4 arguments.");
		exit(-1);
	}

	char* inputFile = argv[1];
	inWave(inputFile);

	char* IRfile = argv[2];
	inIR(IRfile);
    
    int maxSize;
	if (sampleNum > sampleNumIR)
		maxSize = sampleNum;
	else
		maxSize = sampleNumIR;

	int maxSizePow2;
	maxSizePow2 = (int)log2(maxSize) + 1;
	maxSizePow2 = pow(2, maxSizePow2);
	int doubleMaxSize = 2 * maxSizePow2;

	//allocate memory for complex input and IR using double max size
	double *complexWave;
	double *complexIR;
	complexWave = (double*)malloc(sizeof(double) * doubleMaxSize);
	complexIR = (double*)malloc(sizeof(double) * doubleMaxSize);
    
    double maxWave;
    double maxIR;
    
    for (int i = 0; i < sampleNum; i++)
	{
			if (i == 0)
			{
				maxWave = data[0];
			}
			else if (data[i] > maxWave)
			{
				maxWave = data[i];
			}
	}
    
    for (int i = 0; i < sampleNumIR; i++)
	{
			if (i == 0)
			{
				maxIR = dataIR[0];
			}
			else if (dataIR[i] > maxIR)
			{
				maxIR = dataIR[i];
			}
	}
    
    for (int i = 0; i < doubleMaxSize; i++)
	{
		complexIR[i] = 0.0;
		complexWave[i] = 0.0;
	}

	double MAX_VAL = 32767.f;

	for (int i = 0; i < sampleNum; i++)
	{
		complexWave[2 * i] = (((double)data[i]) / maxWave) * MAX_VAL;
	}

	for (int i = 0; i < sampleNumIR; i++)
	{
		complexIR[2 * i] = (((double)dataIR[i]) / maxIR) * MAX_VAL;
	}


	four1(complexWave-1, maxSizePow2, 1);
	four1(complexIR-1, maxSizePow2, 1);

	double *complexOutput;
	complexOutput = (double*)malloc(sizeof(double) *doubleMaxSize);

	for (int i = 0; i < maxSizePow2; i++)
	{
        double wave1 = complexWave[i];
        double wave2 = complexWave[i + 1];
        double IR1 = complexIR[i];
        double IR2 = complexIR[i + 1];
		complexOutput[i * 2] = wave1 * IR1 - wave2 * IR2;
		complexOutput[i * 2 + 1] = wave2 * IR1 + wave1 * IR2;
	}
    
	four1(complexOutput-1, maxSizePow2, 1);
    
    outWave(argv[3], complexOutput, doubleMaxSize);
	end = clock();
	printf("Time: %f seconds", (end - begin) / (double)CLOCKS_PER_SEC);
	free(data);
}