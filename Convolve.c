/**
**Convolve algortihm credited to Leonard Manzara's convolve.c
**/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


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

float* newData;
float max_Sample = -1;
float MAX_VAL = 32767.f;

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

int outWave(char* file, float* newInput, int newSize)
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

		for (int i = 0; i < newSize; ++i)
		{
			newInput[i] = (newInput[i] / max_Sample);
			short sample = (short)(newInput[i] * MAX_VAL);
			fwrite(&sample, 1, bytes, outFile);
		}

		//clean up
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

void convolve(short *x, int N, short *h, int M, float *y, int P){
    for (int i = 0; i < N; i++){
			for (int j = 0; j < M; j++)
			{
				y[i + j] += ((float)x[i] / MAX_VAL) * ((float)h[j] / MAX_VAL);
			}
			if (i == 0)
			{
				max_Sample = y[0];
			}
			else if (y[i] > max_Sample)
			{
				max_Sample = y[i];
			}
	}
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
    
    newData = (float*)malloc(sizeof(float) * (sampleNum + sampleNumIR - 1));
    
    convolve(data, sampleNum, dataIR, sampleNumIR, newData, sampleNum + sampleNumIR - 1);
	
    outWave(argv[3], newData, sampleNum + sampleNumIR - 1);
	end = clock();
	printf("Time spent: %f seconds", (end - begin) / (double)CLOCKS_PER_SEC);
	free(data);
}