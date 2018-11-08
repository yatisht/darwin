#include "DRAM.h"
#include <memory>

DRAM::DRAM()
	: size(4ull * 1024ull * 1024ull * 1024ull),
	referenceSize(0)
{

	buffer = (char*)_aligned_malloc(size, 64);
}


DRAM::~DRAM()
{
}
