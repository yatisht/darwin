#pragma once
#include <stdio.h>
#include <string>
#include <vector>

#include "bond/core/blob.h"

namespace Darwin
{
	namespace Index
	{
		extern std::vector<std::string> chr_id;
		extern std::vector<uint32_t> chr_coord;
		extern std::vector<uint32_t> chr_len;
		extern std::vector<uint32_t> chr_len_unpadded;

		void init();
		bond::blob add_chr(std::string name, bond::blob seq);
	}
}

