#include "Index.h"

#include "Processor.h"

std::vector<std::string>  Darwin::Index::chr_id;
std::vector<uint32_t>  Darwin::Index::chr_len;
std::vector<uint32_t>  Darwin::Index::chr_len_unpadded;
std::vector<uint32_t>  Darwin::Index::chr_coord;

void Darwin::Index::init()
{
	memset(g_DRAM->buffer, 'N', WORD_SIZE);

	g_DRAM->referenceSize += WORD_SIZE;

	chr_coord.push_back(g_DRAM->referenceSize);
}

bond::blob Darwin::Index::add_chr(std::string name, bond::blob seq)
{
	chr_id.push_back(name);

	size_t seq_len = seq.size();

	chr_len.push_back(seq_len);

	std::memcpy(g_DRAM->buffer + g_DRAM->referenceSize, seq.data(), seq_len);

	bond::blob ref(g_DRAM->buffer + g_DRAM->referenceSize, seq_len);

	size_t extra = seq_len % WORD_SIZE;
	if (extra != 0)
	{
		extra = WORD_SIZE - extra;
		memset(g_DRAM->buffer + g_DRAM->referenceSize + seq_len, 'N', extra);

		seq_len += extra;
	}

	g_DRAM->referenceSize += seq_len;

	chr_coord.push_back(g_DRAM->referenceSize);

	return ref;
}
