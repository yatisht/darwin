#pragma once

#define WORD_SIZE 128

class DRAM
{
public:
	char* buffer;
	size_t size;

	size_t referenceSize;
	size_t bufferPosition;
public:
	DRAM();
	~DRAM();
};

