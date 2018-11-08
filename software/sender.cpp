#include "graph.h"
tbb::reader_writer_lock fpga_writer_lock;

seeder_input reference_sender_body::operator()(seeder_input input)
{
	auto &reads = get<0>(input);

	size_t token = get<1>(input);

	size_t word_size = WORD_SIZE;
	size_t max_char_to_send = MAX_CHAR_TO_SEND;
	
	size_t num_bytes_to_send = 0;

	//fpga_writer_lock.lock();
	for (size_t i = 0; i < reads.size(); i++)
	{
		Read &read = reads[i];

		const size_t read_len = read.seq.size();
		char *read_char = (char *)read.seq.data();

		size_t dram_start_addr = read_char - g_DRAM->buffer;

        // SEND READ
		for (size_t j = 0; j < read_len;) {
			num_bytes_to_send = std::min(max_char_to_send, read_len - j);
			if (num_bytes_to_send % word_size != 0) {
				num_bytes_to_send += word_size - (num_bytes_to_send	% word_size);
			}

			assert(num_bytes_to_send % word_size == 0);

			InitializeDRAMMessage init_dram_message;
			GenerateInitializeMemoryMessage(init_dram_message, read_char, j, dram_start_addr, num_bytes_to_send);

			{
				InitializeDRAMMessageResponse init_dram_response;
				g_InitializeReferenceMemory(0, g_DRAM->buffer, init_dram_message, init_dram_response);
			}

			dram_start_addr += num_bytes_to_send;
			j += num_bytes_to_send;
		}
	}
	//fpga_writer_lock.unlock();

	return input;
}


seeder_input read_sender_body::operator()(seeder_input input)
{
	auto &reads = get<0>(input);

	size_t token = get<1>(input);

	size_t word_size = WORD_SIZE;
	size_t max_char_to_send = MAX_CHAR_TO_SEND;

	size_t num_bytes_to_send = 0;

	//fpga_writer_lock.lock();
	for (size_t i = 0; i < reads.size(); i++)
	{
		Read &read = reads[i];

		const size_t read_len = read.seq.size();
		char *read_char = (char *)read.seq.data();

		size_t dram_start_addr = read_char - g_DRAM->buffer;

		// SEND READ
		for (size_t j = 0; j < read_len;) {
			num_bytes_to_send = std::min(max_char_to_send, read_len - j);
			if (num_bytes_to_send % word_size != 0) {
				num_bytes_to_send += word_size - (num_bytes_to_send	% word_size);
			}

			assert(num_bytes_to_send % word_size == 0);

			InitializeDRAMMessage init_dram_message;
			GenerateInitializeMemoryMessage(init_dram_message, read_char, j, dram_start_addr, num_bytes_to_send);

			{
				InitializeDRAMMessageResponse init_dram_response;
				g_InitializeReadMemory(0, g_DRAM->buffer, init_dram_message, init_dram_response);
			}

			dram_start_addr += num_bytes_to_send;
			j += num_bytes_to_send;
		}
	}
	//fpga_writer_lock.unlock();

	return input;
}
