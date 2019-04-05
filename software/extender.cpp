#include <atomic>
#include "graph.h"

std::atomic<int> extender_body::num_extend_tiles(0);
std::atomic<int> extender_body::num_active_tiles(0);
std::atomic<int> extender_body::num_large_tiles(0);

//printer_input extender_body::operator()(extender_input input, extender_node::output_ports_type &op)
void extender_body::operator()(extender_input input, extender_node::output_ports_type &op)
	{
		auto &payload = get<0>(input);

		auto &reads = get<0>(payload);

		auto &data = get<1>(payload);

		size_t token = get<1>(input);

        extend_data output;

		// Forward reads
		{
			auto& f_extend_locations = data.fwLocations;

			int num_extensions_remaining = f_extend_locations.size();
			int first_extension_batch_size = std::min(num_extensions_remaining, cfg.batch_size);

			std::vector<ExtendAlignments> extend_alignments;
			extend_alignments.clear();

			auto loc = f_extend_locations.cbegin();

			// first batch
			for (int b = 0; b < first_extension_batch_size; b++)
			{
				//Add to batch here
				ExtendAlignments extend_alignment = makeForwardAlignment(reads, loc);
                loc++;

				extend_alignments.push_back(extend_alignment);
			}

			int tiles_active = first_extension_batch_size;

			while (num_extensions_remaining > 0) {

				BatchAlignmentInputFieldsDRAM request_batch_dram;
				BatchAlignmentResultDRAM result_batch_dram;
				request_batch_dram.requests.resize(first_extension_batch_size);
				request_batch_dram.do_traceback = 1;

				for (int b = 0; b < first_extension_batch_size; b++) {
					// fill request batch dram here
					ExtendAlignments* e = &extend_alignments[b];

					request_batch_dram.requests[b].index = b;

					if (e->left_extension_done == 0) {
                        int r_tile_size = cfg.tile_size;
                        int q_tile_size = cfg.tile_size;
                        if (e->used_large_tile) {
                            uint64_t hit_offset = e->left_hit_offsets.back();
                            uint64_t h1 = e->reference_start_addr + e->curr_reference_offset;
                            uint64_t o1 = e->curr_query_offset;
                            
                            uint64_t h2 = (hit_offset >> 32);
                            uint64_t o2 = ((hit_offset << 32) >> 32);

                            if ((h1 - h2) > (o1 - o2)) {
                                r_tile_size = 1984;
                                q_tile_size = 960;
                            }
                            else {
                                r_tile_size = 960;
                                q_tile_size = 1984;
                            }
                            num_large_tiles++;
                        }
                        
//                        int num_pop = 0;
//                        for (int k = e->left_hit_offsets.size()-1; k >= 0; k--) {
//                            uint64_t hit_offset = e->left_hit_offsets[k];
//                            uint64_t hit = (hit_offset >> 32);
//                            uint64_t offset = ((hit_offset << 32) >> 32);
//
//                            if ((hit < e->reference_start_addr + e->curr_reference_offset) && (offset < e->curr_query_offset)) {
//                                break;
//                            }
//                            else {
//                                num_pop++;
//                            }
//                        }
//
//                        for (int k = 0; k < num_pop; k++) {
//                            e->left_hit_offsets.pop_back();
//                        }

//                        if (e->left_hit_offsets.size() > 0) {
//                            uint64_t hit_offset = e->left_hit_offsets.back();
//                            
//                            uint64_t h1 = e->reference_start_addr + e->curr_reference_offset;
//                            uint64_t o1 = e->curr_query_offset;
//                            
//                            uint64_t h2 = (hit_offset >> 32);
//                            uint64_t o2 = ((hit_offset << 32) >> 32);
//                            
//                            assert (h1 >= h2);
//                            assert (o1 >= o2);
//                            uint32_t match = ((h1 - h2) < (o1 - o2)) ? (h1 - h2) : (o1 - o2);
//                            uint32_t gap = ((h1 - h2) > (o1 - o2)) ? ((h1 - h2) - (o1 - o2)) : ((o1 - o2) - (h1 - h2));
//							assert(match < 9600);
//							assert(gap < 9600);
//
//                            if (gap > 100) {
//                                tile_size = FIND_MIN(2*(gap+match), 1984);
//								//tile_size = 1400;
//                                std::cerr << "Setting tile size to "<< tile_size << "\n";
//                            }
//                        }

						request_batch_dram.requests[b].ref_size = std::min(e->curr_reference_offset + 1, uint32_t(r_tile_size));
						request_batch_dram.requests[b].query_size = std::min(e->curr_query_offset + 1, uint32_t(q_tile_size));

						request_batch_dram.requests[b].ref_bases_start_addr = e->reference_start_addr + ((e->curr_reference_offset >= r_tile_size) ? (e->curr_reference_offset - r_tile_size + 1) : 0);
						request_batch_dram.requests[b].query_bases_start_addr = e->query_start_addr + ((e->curr_query_offset >= q_tile_size) ? (e->curr_query_offset - q_tile_size + 1) : 0);

						request_batch_dram.requests[b].max_tb_steps = 2*cfg.tile_size;
						request_batch_dram.requests[b].score_threshold = 0;

						request_batch_dram.requests[b].align_fields = start_end;

					}
					else if (e->right_extension_done == 0) {
                        int r_tile_size = cfg.tile_size;
                        int q_tile_size = cfg.tile_size;
                        if (e->used_large_tile) {
                            uint64_t hit_offset = e->right_hit_offsets.back();
                            uint64_t h1 = e->reference_start_addr + e->curr_reference_offset;
                            uint64_t o1 = e->curr_query_offset;
                            
                            uint64_t h2 = (hit_offset >> 32);
                            uint64_t o2 = ((hit_offset << 32) >> 32);

                            if ((h2 - h1) > (o2 - o1)) {
                                r_tile_size = 1984;
                                q_tile_size = 960;
                            }
                            else {
                                r_tile_size = 960;
                                q_tile_size = 1984;
                            }
                            num_large_tiles++;
                        }
                        
//                        int num_pop =  0;
//                        for (int k = e->right_hit_offsets.size()-1; k >= 0; k--) {
//                            uint64_t hit_offset = e->right_hit_offsets[k];
//                            uint64_t hit = (hit_offset >> 32);
//                            uint64_t offset = ((hit_offset << 32) >> 32);
//
//                            if ((hit > e->reference_start_addr + e->curr_reference_offset) && (offset > e->curr_query_offset)) {
//                                break;
//                            }
//                            else {
//                                num_pop++;
//                            }
//                        }
//                        
//                        for (int k = 0; k < num_pop; k++) {
//                            e->right_hit_offsets.pop_back();
//                        }

//                        if (e->right_hit_offsets.size() > 0) {
//                            uint64_t hit_offset = e->right_hit_offsets.back();
//                            
//                            uint64_t h1 = (hit_offset >> 32);
//                            uint64_t o1 = ((hit_offset << 32) >> 32);
//                            
//                            uint64_t h2 = e->reference_start_addr + e->curr_reference_offset;
//                            uint64_t o2 = e->curr_query_offset;
//                            
//                            assert (h1 >= h2);
//                            assert (o1 >= o2);
//                            uint32_t match = ((h1 - h2) < (o1 - o2)) ? (h1 - h2) : (o1 - o2);
//                            uint32_t gap = ((h1 - h2) > (o1 - o2)) ? ((h1 - h2) - (o1 - o2)) : ((o1 - o2) - (h1 - h2));
//							assert(match < 9600);
//							assert(gap < 9600);
//
//                            if (gap > 100) {
//                                tile_size = FIND_MIN(2*(gap+match), 1984);
//								//tile_size = 1400;
//                                std::cerr << "Setting tile size to "<< tile_size << "\n";
//                            }
//                        }
//

						request_batch_dram.requests[b].ref_size = std::min(e->reference_length - e->curr_reference_offset, uint32_t(r_tile_size));
						request_batch_dram.requests[b].query_size = std::min(e->query_length - e->curr_query_offset, uint32_t(q_tile_size));

						request_batch_dram.requests[b].ref_bases_start_addr = e->reference_start_addr + e->curr_reference_offset;
						request_batch_dram.requests[b].query_bases_start_addr = e->query_start_addr + e->curr_query_offset;

						request_batch_dram.requests[b].max_tb_steps = 2*cfg.tile_size;
						request_batch_dram.requests[b].score_threshold = 0;

						request_batch_dram.requests[b].align_fields = reverse_ref + reverse_query + start_end;
					}
					// Completed extension
					else {
						request_batch_dram.requests[b].ref_size = 8;
						request_batch_dram.requests[b].query_size = 8;

						request_batch_dram.requests[b].ref_bases_start_addr = 0;
						request_batch_dram.requests[b].query_bases_start_addr = 0;

						request_batch_dram.requests[b].max_tb_steps = 2*cfg.tile_size;
						request_batch_dram.requests[b].score_threshold = 0;

						request_batch_dram.requests[b].align_fields = 0;
					}
				}

				// process batch 
				{
                    //fpga_writer_lock.lock_read();

					g_BatchAlignmentSIMD(token, g_DRAM->buffer, request_batch_dram, result_batch_dram);

                    //fpga_writer_lock.unlock();
				}

				num_extend_tiles += first_extension_batch_size;
				num_active_tiles += tiles_active;

				for (int b = 0; b < first_extension_batch_size; b++) {
					AlignmentResult result = result_batch_dram.results[b];
					int idx = result.index;
					ExtendAlignments e = extend_alignments[idx];

					int read_num = e.read_num;
					Read &read = reads[read_num];
					const size_t read_len = read.seq.size();
					char *read_char = (char *)read.seq.data();

					// check if extension is done
					if ((e.left_extension_done == 1) && (e.right_extension_done == 1)) {
						//do nothing
					}
					else {
						//read result. modify vector.
						if (e.left_extension_done == 0) {
							size_t left_extension_size = result.TB_pointers.size() * 32;
							char* reference_buffer = (char*)malloc(left_extension_size);
							char* query_buffer = (char*)malloc(left_extension_size);
							size_t reference_pos = left_extension_size;
							size_t query_pos = left_extension_size;

//                            int tile_size = cfg.tile_size;
                            int r_tile_size = cfg.tile_size;
                            int q_tile_size = cfg.tile_size;
                            if (e.used_large_tile && (cfg.do_overlap == 0)) {
                                uint64_t hit_offset = e.left_hit_offsets.back();
                                uint64_t h1 = e.reference_start_addr + e.curr_reference_offset;
                                uint64_t o1 = e.curr_query_offset;

                                uint64_t h2 = (hit_offset >> 32);
                                uint64_t o2 = ((hit_offset << 32) >> 32);

                                if ((h1 - h2) > (o1 - o2)) {
                                    r_tile_size = 1984;
                                    q_tile_size = 960;
                                }
                                else {
                                    r_tile_size = 960;
                                    q_tile_size = 1984;
                                }
                            }
                            int num_tb_steps = 0;

							for (int t = 0; t < result.total_TB_pointers; t += 32) {
                                uint64_t tb = result.TB_pointers[t/32];
                                int num_p = std::min(result.total_TB_pointers - t, 32);
                                for (int p = 0; p < num_p; p++) {
                                    int dir = ((tb >> 2 * p) & TB_MASK);
                                    switch (dir) {
                                        case M: {
                                                    reference_buffer[--reference_pos] = g_DRAM->buffer[e.reference_start_addr + e.curr_reference_offset];
                                                    query_buffer[--query_pos] = read_char[e.curr_query_offset];
                                                    if (e.curr_reference_offset > 0) {
                                                        e.curr_reference_offset--;
                                                    }
                                                    else {
                                                        e.reference_start_offset = 0;
                                                    }
                                                    if (e.curr_query_offset > 0) {
                                                        e.curr_query_offset--;
                                                    }
                                                    else {
                                                        e.query_start_offset = 0;
                                                    }
                                                    break;
                                                }
                                        case I: {
                                                    reference_buffer[--reference_pos] = '-';
                                                    query_buffer[--query_pos] = read_char[e.curr_query_offset];
                                                    if (e.curr_query_offset > 0) {
                                                        e.curr_query_offset--;
                                                    }
                                                    else {
                                                        e.query_start_offset = 0;
                                                    }
                                                    break;
                                                }
                                        case D: {
                                                    reference_buffer[--reference_pos] = g_DRAM->buffer[e.reference_start_addr + e.curr_reference_offset];
                                                    query_buffer[--query_pos] = '-';
                                                    if (e.curr_reference_offset > 0) {
                                                        e.curr_reference_offset--;
                                                    }
                                                    else {
                                                        e.reference_start_offset = 0;
                                                    }
                                                    break;
                                                }
                                    }
                                    num_tb_steps++;
                                    if ((num_tb_steps >= std::min(r_tile_size, q_tile_size) - cfg.tile_overlap) && (dir == M)) {
                                        break;
                                    }
                                }
                            }

                            e.aligned_reference_str = std::string(reference_buffer + reference_pos, left_extension_size - reference_pos) + e.aligned_reference_str;
                            e.aligned_query_str = std::string(query_buffer + query_pos, left_extension_size - query_pos) + e.aligned_query_str;

                            int num_pop = 0;
                            for (int k = e.left_hit_offsets.size()-1; k >= 0; k--) {
                                uint64_t hit_offset = e.left_hit_offsets[k];
                                uint64_t hit = (hit_offset >> 32);
                                uint64_t offset = ((hit_offset << 32) >> 32);

                                if ((hit < e.reference_start_addr + e.curr_reference_offset) && (offset < e.curr_query_offset)) {
                                    break;
                                }
                                else {
                                    num_pop++;
                                }
                            }
                            for (int k = 0; k < num_pop; k++) {
                                e.left_hit_offsets.pop_back();
                            }

							if ((result.total_TB_pointers == 0) || (e.left_hit_offsets.size() == 0) || (e.reference_start_offset == 0) || (e.query_start_offset == 0)) {
                                if ((e.used_large_tile) || (e.left_hit_offsets.size() == 0) || (e.reference_start_offset == 0) || (e.query_start_offset == 0)) {
                                    e.left_extension_done = 1;
                                    if (e.reference_start_offset > 0) {
                                        e.reference_start_offset = e.curr_reference_offset + 1;
                                    }
                                    if (e.query_start_offset > 0) {
                                        e.query_start_offset = e.curr_query_offset + 1;
                                    }
                                    extend_alignments[idx] = e;
                                    if ((e.curr_reference_offset + 1 < e.reference_length) && (e.curr_query_offset + 1 < e.query_length) && (e.right_extension_done == 0)) {
                                        e.curr_reference_offset = e.reference_end_offset + 1;
                                        e.curr_query_offset = e.query_end_offset + 1;
                                        extend_alignments[idx] = e;
                                    }
                                    else {
                                        //fetch next
                                        e.right_extension_done = 1;
                                        num_extensions_remaining--;
                                        extend_alignments[idx] = e;

                                        // check if more elements in queue
                                        if (loc != f_extend_locations.cend()) {
                                            extend_alignments[idx] = makeForwardAlignment(reads, loc);
                                            loc++;
                                        }
                                        else {
                                            tiles_active--;
                                        }
                                    }
                                }
                                else {
                                    e.used_large_tile = true;
                                    //std::cerr << "Using large tile \n";
                                    //std::cerr << e.read_num << " " << e.reference_start_offset + e.reference_start_addr << " " << result.total_TB_pointers << " " << e.curr_query_offset << " " << ((e.left_hit_offsets.front() << 32) >> 32) << std::endl;
                                    extend_alignments[idx] = e;
                                }
							}
							else {
								e.used_large_tile = false;
								extend_alignments[idx] = e;
							}

							free(reference_buffer);
							free(query_buffer);
						}
						else if (e.right_extension_done == 0) {
							size_t right_extension_size = result.TB_pointers.size() * 32;
							char* reference_buffer = (char*)malloc(right_extension_size);
							char* query_buffer = (char*)malloc(right_extension_size);
							size_t reference_pos = 0;
							size_t query_pos = 0;

                            int r_tile_size = cfg.tile_size;
                            int q_tile_size = cfg.tile_size;
                            if (e.used_large_tile && (cfg.do_overlap == 0)) {
                                uint64_t hit_offset = e.right_hit_offsets.back();
                                uint64_t h1 = e.reference_start_addr + e.curr_reference_offset;
                                uint64_t o1 = e.curr_query_offset;

                                uint64_t h2 = (hit_offset >> 32);
                                uint64_t o2 = ((hit_offset << 32) >> 32);

                                if ((h2 - h1) > (o2 - o1)) {
                                    r_tile_size = 1984;
                                    q_tile_size = 960;
                                }
                                else {
                                    r_tile_size = 960;
                                    q_tile_size = 1984;
                                }
                            }
                            int num_tb_steps = 0;

							for (int t = 0; t < result.total_TB_pointers; t += 32) {
								uint64_t tb = result.TB_pointers[t/32];
								int num_p = std::min(result.total_TB_pointers - t, 32);
								for (int p = 0; p < num_p; p++) {
									int dir = ((tb >> 2 * p) & TB_MASK);
									switch (dir) {
									case M: {
										reference_buffer[reference_pos++] = g_DRAM->buffer[e.reference_start_addr + e.curr_reference_offset];
										query_buffer[query_pos++] = read_char[e.curr_query_offset];
										if (e.curr_reference_offset < e.reference_length) {
											e.curr_reference_offset++;
										}
										if (e.curr_query_offset < e.query_length) {
											e.curr_query_offset++;
										}
										break;
									}
									case I: {
										reference_buffer[reference_pos++] = '-';
										query_buffer[query_pos++] = read_char[e.curr_query_offset];
										if (e.curr_query_offset < e.query_length) {
											e.curr_query_offset++;
										}
										break;
									}
									case D: {
										reference_buffer[reference_pos++] = g_DRAM->buffer[e.reference_start_addr + e.curr_reference_offset];
										query_buffer[query_pos++] = '-';
										if (e.curr_reference_offset < e.reference_length) {
											e.curr_reference_offset++;
										}
										break;
									}
									}
                                    num_tb_steps++;
                                    if ((num_tb_steps >= std::min(r_tile_size, q_tile_size) - cfg.tile_overlap) && (dir == M)) {
                                        break;
                                    }
								}
							}


							e.aligned_reference_str = e.aligned_reference_str + std::string(reference_buffer, reference_pos);
							e.aligned_query_str = e.aligned_query_str + std::string(query_buffer, query_pos);
                        
                            int num_pop =  0;
                            for (int k = e.right_hit_offsets.size()-1; k >= 0; k--) {
                                uint64_t hit_offset = e.right_hit_offsets[k];
                                uint64_t hit = (hit_offset >> 32);
                                uint64_t offset = ((hit_offset << 32) >> 32);

                                if ((hit > e.reference_start_addr + e.curr_reference_offset) && (offset > e.curr_query_offset)) {
                                    break;
                                }
                                else {
                                    num_pop++;
                                }
                            }

                            for (int k = 0; k < num_pop; k++) {
                                e.right_hit_offsets.pop_back();
                            }

							if ((result.total_TB_pointers == 0) || (e.curr_reference_offset == e.reference_length) || (e.curr_query_offset == e.query_length)) {
                                if ((e.used_large_tile) || (e.right_hit_offsets.size() == 0)|| (e.curr_reference_offset == e.reference_length) || (e.curr_query_offset == e.query_length)) {
                                    e.reference_end_offset = e.curr_reference_offset - 1;
                                    e.query_end_offset = e.curr_query_offset - 1;

                                    // Send alignment to printer
//                                    get<0>(op).try_put(printer_input(read, e));
//                                    get<0>(op).try_put(printer_input(token));
                                    int score = AlignmentScore(e.aligned_reference_str, e.aligned_query_str);
                                    e.score = score;
                                    output.extend_alignments.push_back(e);

                                    e.right_extension_done = 1;
                                    num_extensions_remaining--;
                                    extend_alignments[idx] = e;

                                    // check if more elements in queue
                                    if (loc != f_extend_locations.cend()) {
                                        extend_alignments[idx] = makeForwardAlignment(reads, loc);
                                        loc++;
                                    }
                                    else{
                                        tiles_active--;
                                    }
                                }
                                else {
                                    e.used_large_tile = true;
                                    //std::cerr << "Using large tile \n";
                                    extend_alignments[idx] = e;
                                }
                            }
							else {
								e.used_large_tile = false;
								extend_alignments[idx] = e;
							}
							free(reference_buffer);
							free(query_buffer);
						}
					}
				}
			}

			assert(tiles_active == 0);
		}

		// Reverse-complement read
		{
			auto& rc_extend_locations = data.rcLocations;

			int num_extensions_remaining = rc_extend_locations.size();
			int first_extension_batch_size = std::min(num_extensions_remaining, cfg.batch_size);

			std::vector<ExtendAlignments> extend_alignments;
			extend_alignments.clear();

			auto loc = rc_extend_locations.cbegin();

			// first batch
			for (int b = 0; b < first_extension_batch_size; b++) {
				//Add to batch here
				ExtendAlignments extend_alignment = makeBackwardAlignment(reads, loc);
                loc++;
				extend_alignments.push_back(extend_alignment);
			}

			int tiles_active = first_extension_batch_size;

			while (num_extensions_remaining > 0) {

				BatchAlignmentInputFieldsDRAM request_batch_dram;
				BatchAlignmentResultDRAM result_batch_dram;
				request_batch_dram.requests.resize(first_extension_batch_size);
				request_batch_dram.do_traceback = 1;

				for (int b = 0; b < first_extension_batch_size; b++) {
					// fill request batch dram here
					ExtendAlignments* e = &extend_alignments[b];

					int read_num = e->read_num;
					Read &read = reads[read_num];
					const size_t read_len = read.seq.size();

					request_batch_dram.requests[b].index = b;
					if (e->left_extension_done == 0) {

                        int r_tile_size = cfg.tile_size;
                        int q_tile_size = cfg.tile_size;
                        if (e->used_large_tile) {
                            uint64_t hit_offset = e->left_hit_offsets.back();
                            uint64_t h1 = e->reference_start_addr + e->curr_reference_offset;
                            uint64_t o1 = e->curr_query_offset;
                            
                            uint64_t h2 = (hit_offset >> 32);
                            uint64_t o2 = ((hit_offset << 32) >> 32);

                            if ((h1 - h2) > (o1 - o2)) {
                                r_tile_size = 1984;
                                q_tile_size = 960;
                            }
                            else {
                                r_tile_size = 960;
                                q_tile_size = 1984;
                            }
                            num_large_tiles++;
                        }
                        
//                        int num_pop = 0;
//                        for (int k = e->left_hit_offsets.size()-1; k >= 0; k--) {
//                            uint64_t hit_offset = e->left_hit_offsets[k];
//                            uint64_t hit = (hit_offset >> 32);
//                            uint64_t offset = ((hit_offset << 32) >> 32);
//
//                            if ((hit < e->reference_start_addr + e->curr_reference_offset) && (offset < e->curr_query_offset)) {
//                                break;
//                            }
//                            else {
//                                num_pop++;
//                            }
//                        }
//
//                        for (int k = 0; k < num_pop; k++) {
//                            e->left_hit_offsets.pop_back();
//                        }

//                        if (e->left_hit_offsets.size() > 0) {
//                            uint64_t hit_offset = e->left_hit_offsets.back();
//                            
//                            uint64_t h1 = e->reference_start_addr + e->curr_reference_offset;
//                            uint64_t o1 = e->curr_query_offset;
//                            
//                            uint64_t h2 = (hit_offset >> 32);
//                            uint64_t o2 = ((hit_offset << 32) >> 32);
//                            
//                            assert (h1 >= h2);
//                            assert (o1 >= o2);
//                            uint32_t match = ((h1 - h2) < (o1 - o2)) ? (h1 - h2) : (o1 - o2);
//                            uint32_t gap = ((h1 - h2) > (o1 - o2)) ? ((h1 - h2) - (o1 - o2)) : ((o1 - o2) - (h1 - h2));
//							assert(match < 9600);
//							assert(gap < 9600);
//
//                            if (gap > 100) {
//                                tile_size = FIND_MIN(2*(gap+match), 1984);
//								//tile_size = 1400;
//                                std::cerr << "Setting tile size to "<< tile_size << "\n";
//                            }
//                        }

						request_batch_dram.requests[b].ref_size = std::min(e->curr_reference_offset + 1, uint32_t(r_tile_size));
						request_batch_dram.requests[b].query_size = std::min(e->curr_query_offset + 1,  uint32_t(q_tile_size));

						request_batch_dram.requests[b].ref_bases_start_addr = e->reference_start_addr + ((e->curr_reference_offset >= r_tile_size) ? (e->curr_reference_offset - r_tile_size + 1) : 0);
						request_batch_dram.requests[b].query_bases_start_addr = e->query_start_addr + read_len - request_batch_dram.requests[b].query_size - ((e->curr_query_offset >= q_tile_size) ? (e->curr_query_offset - q_tile_size + 1) : 0);

						request_batch_dram.requests[b].max_tb_steps = 2*cfg.tile_size;
						request_batch_dram.requests[b].score_threshold = 0;

						request_batch_dram.requests[b].align_fields = reverse_query + complement_query + start_end;

					}
					else if (e->right_extension_done == 0) {
                        int r_tile_size = cfg.tile_size;
                        int q_tile_size = cfg.tile_size;
                        if (e->used_large_tile) {
                            uint64_t hit_offset = e->right_hit_offsets.back();
                            uint64_t h1 = e->reference_start_addr + e->curr_reference_offset;
                            uint64_t o1 = e->curr_query_offset;
                            
                            uint64_t h2 = (hit_offset >> 32);
                            uint64_t o2 = ((hit_offset << 32) >> 32);

                            if ((h2 - h1) > (o2 - o1)) {
                                r_tile_size = 1984;
                                q_tile_size = 960;
                            }
                            else {
                                r_tile_size = 960;
                                q_tile_size = 1984;
                            }
                            num_large_tiles++;
                        }
                        
//                        int num_pop = 0;
//                        for (int k = e->right_hit_offsets.size()-1; k >= 0; k--) {
//                            uint64_t hit_offset = e->right_hit_offsets[k];
//                            uint64_t hit = (hit_offset >> 32);
//                            uint64_t offset = ((hit_offset << 32) >> 32);
//
//                            if ((hit > e->reference_start_addr + e->curr_reference_offset) && (offset > e->curr_query_offset)) {
//                                break;
//                            }
//                            else {
//                                num_pop++;
//                            }
//                        }
//
//
//                        for (int k = 0; k < num_pop; k++) {
//                            e->right_hit_offsets.pop_back();
//                        }

//                        if (e->right_hit_offsets.size() > 0) {
//                            uint64_t hit_offset = e->right_hit_offsets.back();
//                            
//                            uint64_t h1 = (hit_offset >> 32);
//                            uint64_t o1 = ((hit_offset << 32) >> 32);
//                            
//                            uint64_t h2 = e->reference_start_addr + e->curr_reference_offset;
//                            uint64_t o2 = e->curr_query_offset;
//                            
//                            assert (h1 >= h2);
//                            assert (o1 >= o2);
//                            uint32_t match = ((h1 - h2) < (o1 - o2)) ? (h1 - h2) : (o1 - o2);
//                            uint32_t gap = ((h1 - h2) > (o1 - o2)) ? ((h1 - h2) - (o1 - o2)) : ((o1 - o2) - (h1 - h2));
//							assert(match < 9600);
//							assert(gap < 9600);
//                            if (gap > 100) {
//                                tile_size = FIND_MIN(2*(gap+match), 1984);
//								//tile_size = 1400; 
//								std::cerr << "Setting tile size to "<< tile_size << "\n";
//                            }
//                        }

						request_batch_dram.requests[b].ref_size = std::min(e->reference_length - e->curr_reference_offset, uint32_t(r_tile_size));
						request_batch_dram.requests[b].query_size = std::min(e->query_length - e->curr_query_offset, uint32_t(q_tile_size));

						request_batch_dram.requests[b].ref_bases_start_addr = e->reference_start_addr + e->curr_reference_offset;
						request_batch_dram.requests[b].query_bases_start_addr = e->query_start_addr + read_len - request_batch_dram.requests[b].query_size - e->curr_query_offset;

						request_batch_dram.requests[b].max_tb_steps = 2*cfg.tile_size;
						request_batch_dram.requests[b].score_threshold = 0;

						request_batch_dram.requests[b].align_fields = reverse_ref + complement_query + start_end;
					}
					// Completed extension
					else {
						request_batch_dram.requests[b].ref_size = 8;
						request_batch_dram.requests[b].query_size = 8;

						request_batch_dram.requests[b].ref_bases_start_addr = 0;
						request_batch_dram.requests[b].query_bases_start_addr = 0;

						request_batch_dram.requests[b].max_tb_steps = 2*cfg.tile_size;
						request_batch_dram.requests[b].score_threshold = 0;

						request_batch_dram.requests[b].align_fields = 0;
					}
				}

				// process batch 
				{
                    //fpga_writer_lock.lock_read();

					g_BatchAlignmentSIMD(token, g_DRAM->buffer, request_batch_dram, result_batch_dram);

                    //fpga_writer_lock.unlock();
				}

				num_extend_tiles += first_extension_batch_size;
				num_active_tiles += tiles_active;

				for (int b = 0; b < first_extension_batch_size; b++) {
					AlignmentResult result = result_batch_dram.results[b];
					int idx = result.index;
					ExtendAlignments e = extend_alignments[idx];

					int read_num = e.read_num;
					Read &read = reads[read_num];
					const size_t read_len = read.seq.size();
					char *rev_read_char = (char *)read.rc_seq.data();

					// check if extension is done
					if ((e.left_extension_done == 1) && (e.right_extension_done == 1)) {
						//do nothing
					}
					else {
						//read result. modify vector.
						if (e.left_extension_done == 0) {
							size_t left_extension_size = result.TB_pointers.size() * 32;
							char* reference_buffer = (char*)malloc(left_extension_size);
							char* query_buffer = (char*)malloc(left_extension_size);
							size_t reference_pos = left_extension_size;
							size_t query_pos = left_extension_size;

                            int r_tile_size = cfg.tile_size;
                            int q_tile_size = cfg.tile_size;
                            if (e.used_large_tile && (cfg.do_overlap == 0)) {
                                uint64_t hit_offset = e.left_hit_offsets.back();
                                uint64_t h1 = e.reference_start_addr + e.curr_reference_offset;
                                uint64_t o1 = e.curr_query_offset;

                                uint64_t h2 = (hit_offset >> 32);
                                uint64_t o2 = ((hit_offset << 32) >> 32);

                                if ((h1 - h2) > (o1 - o2)) {
                                    r_tile_size = 1984;
                                    q_tile_size = 960;
                                }
                                else {
                                    r_tile_size = 960;
                                    q_tile_size = 1984;
                                }
                            }
                            int num_tb_steps = 0;

							for (int t = 0; t < result.total_TB_pointers; t += 32) {
								uint64_t tb = result.TB_pointers[t/32];
								int num_p = std::min(result.total_TB_pointers - t, 32);
								for (int p = 0; p < num_p; p++) {
									int dir = ((tb >> 2 * p) & TB_MASK);
									switch (dir) {
									case M: {
										reference_buffer[--reference_pos] = g_DRAM->buffer[e.reference_start_addr + e.curr_reference_offset];
										query_buffer[--query_pos] = rev_read_char[e.curr_query_offset];
										if (e.curr_reference_offset > 0) {
											e.curr_reference_offset--;
										}
										else {
											e.reference_start_offset = 0;
										}
										if (e.curr_query_offset > 0) {
											e.curr_query_offset--;
										}
										else {
											e.query_start_offset = 0;
										}
										break;
									}
									case I: {
										reference_buffer[--reference_pos] = '-';
										query_buffer[--query_pos] = rev_read_char[e.curr_query_offset];
										if (e.curr_query_offset > 0) {
											e.curr_query_offset--;
										}
										else {
											e.query_start_offset = 0;
										}
										break;
									}
									case D: {
										reference_buffer[--reference_pos] = g_DRAM->buffer[e.reference_start_addr + e.curr_reference_offset];
										query_buffer[--query_pos] = '-';
										if (e.curr_reference_offset > 0) {
											e.curr_reference_offset--;
										}
										else {
											e.reference_start_offset = 0;
										}
										break;
									}
									}
                                    num_tb_steps++;
                                    if ((num_tb_steps >= std::min(r_tile_size, q_tile_size) - cfg.tile_overlap) && (dir == M)) {
                                        break;
                                    }
								}
							}

							e.aligned_reference_str = std::string(reference_buffer + reference_pos, left_extension_size - reference_pos) + e.aligned_reference_str;
							e.aligned_query_str = std::string(query_buffer + query_pos, left_extension_size - query_pos) + e.aligned_query_str;

                            int num_pop = 0;
                            for (int k = e.left_hit_offsets.size()-1; k >= 0; k--) {
                                uint64_t hit_offset = e.left_hit_offsets[k];
                                uint64_t hit = (hit_offset >> 32);
                                uint64_t offset = ((hit_offset << 32) >> 32);

                                if ((hit < e.reference_start_addr + e.curr_reference_offset) && (offset < e.curr_query_offset)) {
                                    break;
                                }
                                else {
                                    num_pop++;
                                }
                            }
                            for (int k = 0; k < num_pop; k++) {
                                e.left_hit_offsets.pop_back();
                            }

							if ((result.total_TB_pointers == 0) || (e.reference_start_offset == 0) || (e.query_start_offset == 0)) {
                                if ((e.used_large_tile) || (e.left_hit_offsets.size() == 0)|| (e.reference_start_offset == 0) || (e.query_start_offset == 0)) {
                                    e.left_extension_done = 1;
                                    if (e.reference_start_offset > 0) {
                                        e.reference_start_offset = e.curr_reference_offset + 1;
                                    }
                                    if (e.query_start_offset > 0) {
                                        e.query_start_offset = e.curr_query_offset + 1;
                                    }
                                    extend_alignments[idx] = e;
                                    if ((e.curr_reference_offset + 1 < e.reference_length) && (e.curr_query_offset + 1 < e.query_length) && (e.right_extension_done == 0)) {
                                        e.curr_reference_offset = e.reference_end_offset + 1;
                                        e.curr_query_offset = e.query_end_offset + 1;
                                        extend_alignments[idx] = e;
                                    }
                                    else {
                                        // Send alignment to printer
//                                        get<0>(op).try_put(printer_input(read, e));
//                                        get<0>(op).try_put(printer_input(token));
                                        int score = AlignmentScore(e.aligned_reference_str, e.aligned_query_str);
                                        e.score = score;
                                        output.extend_alignments.push_back(e);

                                        e.right_extension_done = 1;
                                        num_extensions_remaining--;
                                        extend_alignments[idx] = e;

                                        // check if more elements in queue
                                        if (loc != rc_extend_locations.cend()) {
                                            // modify request vector
                                            extend_alignments[idx] = makeBackwardAlignment(reads, loc);
                                            loc++;
                                        }
                                        else {
                                            tiles_active--;
                                        }
                                    }
                                }
                                else {
                                    e.used_large_tile = true;
                                    //std::cerr << "Using large tile \n";
                                    extend_alignments[idx] = e;
                                }
							}
							else {
								e.used_large_tile = false;
								extend_alignments[idx] = e;
							}

							free(reference_buffer);
							free(query_buffer);
						}
						else if (e.right_extension_done == 0) {
							size_t right_extension_size = result.TB_pointers.size() * 32;
							char* reference_buffer = (char*)malloc(right_extension_size);
							char* query_buffer = (char*)malloc(right_extension_size);
							size_t reference_pos = 0;
							size_t query_pos = 0;

                            int r_tile_size = cfg.tile_size;
                            int q_tile_size = cfg.tile_size;
                            if (e.used_large_tile && (cfg.do_overlap == 0)) {
                                uint64_t hit_offset = e.right_hit_offsets.back();
                                uint64_t h1 = e.reference_start_addr + e.curr_reference_offset;
                                uint64_t o1 = e.curr_query_offset;

                                uint64_t h2 = (hit_offset >> 32);
                                uint64_t o2 = ((hit_offset << 32) >> 32);

                                if ((h2 - h1) > (o2 - o1)) {
                                    r_tile_size = 1984;
                                    q_tile_size = 960;
                                }
                                else {
                                    r_tile_size = 960;
                                    q_tile_size = 1984;
                                }
                            }
                            int num_tb_steps = 0;

							for (int t = 0; t < result.total_TB_pointers; t += 32) {
								uint64_t tb = result.TB_pointers[t/32];
								int num_p = std::min(result.total_TB_pointers - t, 32);
								for (int p = 0; p < num_p; p++) {
									int dir = ((tb >> 2 * p) & TB_MASK);
									switch (dir) {
									case M: {
										reference_buffer[reference_pos++] = g_DRAM->buffer[e.reference_start_addr + e.curr_reference_offset];
										query_buffer[query_pos++] = rev_read_char[e.curr_query_offset];
										if (e.curr_reference_offset < e.reference_length) {
											e.curr_reference_offset++;
										}
										if (e.curr_query_offset < e.query_length) {
											e.curr_query_offset++;
										}
										break;
									}
									case I: {
										reference_buffer[reference_pos++] = '-';
										query_buffer[query_pos++] = rev_read_char[e.curr_query_offset];
										if (e.curr_query_offset < e.query_length) {
											e.curr_query_offset++;
										}
										break;
									}
									case D: {
										reference_buffer[reference_pos++] = g_DRAM->buffer[e.reference_start_addr + e.curr_reference_offset];
										query_buffer[query_pos++] = '-';
										if (e.curr_reference_offset < e.reference_length) {
											e.curr_reference_offset++;
										}
										break;
									}
									}
                                    num_tb_steps++;
                                    if ((num_tb_steps >= std::min(r_tile_size, q_tile_size) - cfg.tile_overlap) && (dir == M)) {
                                        break;
                                    }
								}
							}


							e.aligned_reference_str = e.aligned_reference_str + std::string(reference_buffer, reference_pos);
							e.aligned_query_str = e.aligned_query_str + std::string(query_buffer, query_pos);
                        
                            int num_pop =  0;
                            for (int k = e.right_hit_offsets.size()-1; k >= 0; k--) {
                                uint64_t hit_offset = e.right_hit_offsets[k];
                                uint64_t hit = (hit_offset >> 32);
                                uint64_t offset = ((hit_offset << 32) >> 32);

                                if ((hit > e.reference_start_addr + e.curr_reference_offset) && (offset > e.curr_query_offset)) {
                                    break;
                                }
                                else {
                                    num_pop++;
                                }
                            }

                            for (int k = 0; k < num_pop; k++) {
                                e.right_hit_offsets.pop_back();
                            }

                            if ((result.total_TB_pointers == 0) || (e.curr_reference_offset == e.reference_length) || (e.curr_query_offset == e.query_length)) {
                                if ((e.used_large_tile) || (e.right_hit_offsets.size() == 0)|| (e.curr_reference_offset == e.reference_length) || (e.curr_query_offset == e.query_length)) {
                                    e.reference_end_offset = e.curr_reference_offset - 1;
                                    e.query_end_offset = e.curr_query_offset - 1;

                                    // Send alignment to printer
//                                    get<0>(op).try_put(printer_input(read, e));
//                                    get<0>(op).try_put(printer_input(token));
                                    int score = AlignmentScore(e.aligned_reference_str, e.aligned_query_str);
                                    e.score = score;
                                    output.extend_alignments.push_back(e);

                                    e.right_extension_done = 1;
                                    num_extensions_remaining--;
                                    extend_alignments[idx] = e;

                                    // check if more elements in queue
                                    if (loc != rc_extend_locations.cend()) {
                                        // modify request vector
                                        extend_alignments[idx] = makeBackwardAlignment(reads, loc);
                                        loc++;
                                    }
                                    else {
                                        tiles_active--;
                                    }
                                }
                                else {
                                    e.used_large_tile = true;
                                    //std::cerr << "Using large tile \n";
                                    extend_alignments[idx] = e;
                                }
                            }
                            else {
                                e.used_large_tile = false;
                                extend_alignments[idx] = e;
                            }
                            free(reference_buffer);
                            free(query_buffer);
                        }
                    }
                }
            }

            assert(tiles_active == 0);
        }

//        for (size_t i = 0; i < reads.size(); i++)
//        {
//            Read &read = reads[i];
//
//        }

        get<1>(op).try_put(token);
        get<0>(op).try_put(printer_input(printer_payload(reads, output), token));
        //		return printer_input(printer_payload(reads, output), token);
    }

ExtendAlignments  extender_body::makeForwardAlignment(std::vector<Read> &batch, std::vector<ExtendLocations>::const_iterator &loc)
{
    int read_num = loc->read_num;

    Read &read = batch[read_num];

    const size_t read_len = read.seq.size();
    char *read_char = (char *)read.seq.data();

    int chr_id = loc->chr_id;
    uint32_t chr_start = Index::chr_coord[chr_id];

    ExtendAlignments extend_alignment;

    extend_alignment.read_num = read_num;
    extend_alignment.chr_id = chr_id;
    extend_alignment.curr_reference_offset = loc->reference_pos - chr_start;
    extend_alignment.curr_query_offset = loc->query_pos;
    extend_alignment.reference_start_offset = loc->reference_pos - chr_start;
    extend_alignment.query_start_offset = loc->query_pos;
    extend_alignment.reference_end_offset = loc->reference_pos - chr_start;
    extend_alignment.query_end_offset = loc->query_pos;

    extend_alignment.reference_start_addr = Index::chr_coord[chr_id];
    extend_alignment.query_start_addr = (read_char - g_DRAM->buffer);
    extend_alignment.reference_length = Index::chr_len[chr_id];
    extend_alignment.query_length = read_len;

    extend_alignment.left_hit_offsets.assign(loc->left_hit_offsets.begin(), loc->left_hit_offsets.end());
    extend_alignment.right_hit_offsets.assign(loc->right_hit_offsets.begin(), loc->right_hit_offsets.end());

    extend_alignment.left_extension_done = 0;
    extend_alignment.right_extension_done = 0;

    extend_alignment.do_print = true;
    extend_alignment.used_large_tile = false;

    extend_alignment.aligned_reference_str = "";
    extend_alignment.aligned_query_str = "";

    extend_alignment.score = 0;

    extend_alignment.strand = '+';

    return extend_alignment;
}

ExtendAlignments  extender_body::makeBackwardAlignment(std::vector<Read> &batch, std::vector<ExtendLocations>::const_iterator &loc)
{
    int read_num = loc->read_num;

    Read &read = batch[read_num];

    const size_t read_len = read.seq.size();
    char *read_char = (char *)read.seq.data();

    int chr_id = loc->chr_id;
    uint32_t chr_start = Index::chr_coord[chr_id];

    ExtendAlignments extend_alignment;

    extend_alignment.read_num = read_num;
    extend_alignment.chr_id = chr_id;
    extend_alignment.curr_reference_offset = loc->reference_pos - chr_start;
    extend_alignment.curr_query_offset = loc->query_pos;
    extend_alignment.reference_start_offset = loc->reference_pos - chr_start;
    extend_alignment.query_start_offset = loc->query_pos;
    extend_alignment.reference_end_offset = loc->reference_pos - chr_start;
    extend_alignment.query_end_offset = loc->query_pos;

    extend_alignment.reference_start_addr = chr_start;
    extend_alignment.query_start_addr = (read_char - g_DRAM->buffer);
    extend_alignment.reference_length = Index::chr_len[chr_id];
    extend_alignment.query_length = read_len;

    extend_alignment.left_hit_offsets.assign(loc->left_hit_offsets.begin(), loc->left_hit_offsets.end());
    extend_alignment.right_hit_offsets.assign(loc->right_hit_offsets.begin(), loc->right_hit_offsets.end());

    extend_alignment.left_extension_done = 0;
    extend_alignment.right_extension_done = 0;

    extend_alignment.do_print = true;
    extend_alignment.used_large_tile = false;

    extend_alignment.aligned_reference_str = "";
    extend_alignment.aligned_query_str = "";

    extend_alignment.score = 0;

    extend_alignment.strand = '-';

    return extend_alignment;
}

int extender_body::AlignmentScore(std::string ref, std::string query) {
    int score = 0;
    int index = 0;
    int open = 0;
    int short_gap_penalty = 0;
    int long_gap_penalty = 0;

    int mat_offset[] = { 0, 1, 3, 6 };

    for (int l = 0; l < ref.length(); l++) {
        char r = ref[l];
        char q = query[l];
        if ((r == '-') || (q == '-')) {
            short_gap_penalty += (open) ? cfg.gap_extend : cfg.gap_open;
            long_gap_penalty += (open) ? cfg.long_gap_extend : cfg.long_gap_open;
            open = 1;
        }
        else {
            int r_nt = NtChar2Int(r);
            int q_nt = NtChar2Int(q);
            if ((r_nt <= 3) && (q_nt <= 3)) {
                if (r_nt > q_nt) {
                    index = q_nt * 4 + r_nt - mat_offset[q_nt];
                }
                else {
                    index = r_nt * 4 + q_nt - mat_offset[r_nt];
                }
                score += cfg.gact_sub_mat[index];
            }
            else {
                score += cfg.gact_sub_mat[10];
            }
            score += (long_gap_penalty < short_gap_penalty) ? short_gap_penalty : long_gap_penalty;
            open = 0;
            short_gap_penalty = 0;
            long_gap_penalty = 0;
        }
    }
    return score;
}

