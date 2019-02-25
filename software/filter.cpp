#include <atomic>
#include "graph.h"

std::atomic<int> filter_body::num_filter_tiles(0);
std::atomic<int> filter_body::num_extend_requests(0);
std::atomic<int> filter_body::num_slope_filtered(0);

extender_input filter_body::operator()(filter_input input)
	{
		auto &payload = get<0>(input);

		auto &reads = get<0>(payload);

		auto &data = get<1>(payload);

		size_t token = get<1>(input);

		filter_data output;

		std::deque<ExtendLocations> read_extend_locations;
		size_t word_size = WORD_SIZE;

		const size_t max_requests = cfg.first_tile_batch_size;

		// Forward reads
		auto& f_anchors = data.fwAnchors;

		for (size_t b = 0; b < f_anchors.size(); b += max_requests)
		{
			size_t first_hit = b;
			size_t last_hit = std::min(b + max_requests, f_anchors.size());

			size_t num_requests = last_hit - first_hit;
			if (num_requests > 0) {

				BatchAlignmentInputFieldsDRAM  request_batch_dram;
				BatchAlignmentResultDRAM result_batch_dram;
				request_batch_dram.requests.resize(num_requests);

				request_batch_dram.do_traceback = 0; //traceback;

				for (size_t c = first_hit, r = 0; c < last_hit; c++, r++)
				{
					uint32_t hit = (f_anchors[c].hit_offset >> 32);
					uint32_t offset = ((f_anchors[c].hit_offset << 32) >> 32);

					size_t chr_id = std::upper_bound(Index::chr_coord.cbegin(), Index::chr_coord.cend(), hit) - Index::chr_coord.cbegin() - 1;
					uint32_t chr_start = Index::chr_coord[chr_id];
					uint32_t chr_end = chr_start + Index::chr_len[chr_id];

					size_t read_num = std::upper_bound(data.fwAnchorBuckets.cbegin(), data.fwAnchorBuckets.cend(), c) - data.fwAnchorBuckets.cbegin() - 1;
					Read &read = reads[read_num];
					const size_t read_len = read.seq.size();
					char *read_char = (char *)read.seq.data();

					uint32_t ref_tile_start = (hit + cfg.first_tile_size < chr_end) ? hit : ((chr_end > cfg.first_tile_size) ? chr_end - cfg.first_tile_size : 0);
					uint32_t query_tile_start = (offset + cfg.first_tile_size < read_len) ? offset : ((read_len > cfg.first_tile_size) ? read_len - cfg.first_tile_size : 0);
					uint32_t ref_tile_size = std::min(uint32_t(cfg.first_tile_size), (chr_end - chr_start));
					uint32_t query_tile_size = std::min(size_t(cfg.first_tile_size), read_len);

					request_batch_dram.requests[r].index = c - b;
					request_batch_dram.requests[r].ref_size = ref_tile_size;
					request_batch_dram.requests[r].query_size = query_tile_size;

					request_batch_dram.requests[r].ref_bases_start_addr = ref_tile_start;
					request_batch_dram.requests[r].query_bases_start_addr = (read_char - g_DRAM->buffer) + query_tile_start;

					request_batch_dram.requests[r].max_tb_steps = 2 * cfg.first_tile_size;
					request_batch_dram.requests[r].score_threshold = 0;

					request_batch_dram.requests[r].align_fields = 0;
				}


                //fpga_writer_lock.lock_read();

				g_BatchAlignmentSIMD(token, g_DRAM->buffer, request_batch_dram, result_batch_dram);

                //fpga_writer_lock.unlock();


				num_filter_tiles += num_requests;


				for (size_t r = 0; r < num_requests; r++) {
					AlignmentResult result = result_batch_dram.results[r];
					int index = b + result.index;

					if (result.score >= cfg.first_tile_score_threshold) {
						uint32_t hit = (f_anchors[index].hit_offset >> 32);
						uint32_t offset = ((f_anchors[index].hit_offset << 32) >> 32);

						size_t chr_id = std::upper_bound(Index::chr_coord.cbegin(), Index::chr_coord.cend(), hit) - Index::chr_coord.cbegin() - 1;
						uint32_t chr_start = Index::chr_coord[chr_id];
						uint32_t chr_end = chr_start + Index::chr_len[chr_id];

						size_t read_num = std::upper_bound(data.fwAnchorBuckets.cbegin(), data.fwAnchorBuckets.cend(), index) - data.fwAnchorBuckets.cbegin() - 1;
						Read &read = reads[read_num];
						const size_t read_len = read.seq.size();

						uint32_t ref_tile_start = (hit + cfg.first_tile_size < chr_end) ? hit : ((chr_end > cfg.first_tile_size) ? chr_end - cfg.first_tile_size : 0);
						uint32_t query_tile_start = (offset + cfg.first_tile_size < read_len) ? offset : ((read_len > cfg.first_tile_size) ? read_len - cfg.first_tile_size : 0);
						uint32_t query_tile_size = std::min(size_t(cfg.first_tile_size), read_len);

						uint32_t ovl = (offset + (chr_end - hit));

						if (ovl > cfg.min_overlap / 2) {
							ExtendLocations loc;
							loc.read_num = read_num;
							loc.chr_id = chr_id;
							loc.score = result.score;
							loc.reference_pos = ref_tile_start + result.ref_max_pos;
							loc.query_pos = query_tile_start + result.query_max_pos;
                            loc.left_hit_offsets.assign(f_anchors[index].left_chained_hits.begin(), f_anchors[index].left_chained_hits.end());
							loc.right_hit_offsets.assign(f_anchors[index].right_chained_hits.begin(), f_anchors[index].right_chained_hits.end());
							read_extend_locations.push_back(loc);
						}
						num_extend_requests += 1;
					}
				}
			}
		}

		slopeFilter(read_extend_locations, output.fwLocations);

		// Reverse-complement reads
		auto& rc_anchors = data.rcAnchors;

		read_extend_locations.clear();

		for (size_t b = 0; b < rc_anchors.size(); b += max_requests)
		{
			size_t first_hit = b;
			size_t last_hit = std::min(b + max_requests, rc_anchors.size());

			size_t num_requests = last_hit - first_hit;

			if (num_requests > 0) {
				BatchAlignmentInputFieldsDRAM  request_batch_dram;
				BatchAlignmentResultDRAM result_batch_dram;
				request_batch_dram.requests.resize(num_requests);

				request_batch_dram.do_traceback = 0; //traceback;

				for (size_t c = first_hit, r = 0; c < last_hit; c++, r++)
				{
					uint32_t hit = (rc_anchors[c].hit_offset >> 32);
					uint32_t offset = ((rc_anchors[c].hit_offset << 32) >> 32);

					size_t chr_id = std::upper_bound(Index::chr_coord.cbegin(), Index::chr_coord.cend(), hit) - Index::chr_coord.cbegin() - 1;
					uint32_t chr_start = Index::chr_coord[chr_id];
					uint32_t chr_end = chr_start + Index::chr_len[chr_id];

					size_t read_num = std::upper_bound(data.rcAnchorBuckets.cbegin(), data.rcAnchorBuckets.cend(), c) - data.rcAnchorBuckets.cbegin() - 1;
					Read &read = reads[read_num];
					const size_t read_len = read.seq.size();
					char *read_char = (char *)read.seq.data();

					uint32_t ref_tile_start = (hit + cfg.first_tile_size < chr_end) ? hit : ((chr_end > cfg.first_tile_size) ? chr_end - cfg.first_tile_size : 0);
					uint32_t query_tile_start = (offset + cfg.first_tile_size < read_len) ? offset : ((read_len > cfg.first_tile_size) ? read_len - cfg.first_tile_size : 0);
					uint32_t ref_tile_size = std::min(uint32_t(cfg.first_tile_size), (chr_end - chr_start));
					uint32_t query_tile_size = std::min(size_t(cfg.first_tile_size), read_len);

					request_batch_dram.requests[r].index = c - b;
					request_batch_dram.requests[r].ref_size = ref_tile_size;
					request_batch_dram.requests[r].query_size = query_tile_size;

					request_batch_dram.requests[r].ref_bases_start_addr = ref_tile_start;
					request_batch_dram.requests[r].query_bases_start_addr = (read_char - g_DRAM->buffer) + read_len - (query_tile_start + query_tile_size);

					request_batch_dram.requests[r].max_tb_steps = 2 * cfg.first_tile_size;
					request_batch_dram.requests[r].score_threshold = 0;

					request_batch_dram.requests[r].align_fields = reverse_query + complement_query;
				}

                //fpga_writer_lock.lock_read();

				g_BatchAlignmentSIMD(token, g_DRAM->buffer, request_batch_dram, result_batch_dram);

                //fpga_writer_lock.unlock();

				num_filter_tiles += num_requests;

				for (int r = 0; r < num_requests; r++) {
					AlignmentResult result = result_batch_dram.results[r];
					int index = b + result.index;

					if (result.score >= cfg.first_tile_score_threshold) {
						uint32_t hit = (rc_anchors[index].hit_offset >> 32);
						uint32_t offset = ((rc_anchors[index].hit_offset << 32) >> 32);

						size_t chr_id = std::upper_bound(Index::chr_coord.cbegin(), Index::chr_coord.cend(), hit) - Index::chr_coord.cbegin() - 1;
						uint32_t chr_start = Index::chr_coord[chr_id];
						uint32_t chr_end = chr_start + Index::chr_len[chr_id];

						size_t read_num = std::upper_bound(data.rcAnchorBuckets.cbegin(), data.rcAnchorBuckets.cend(), index) - data.rcAnchorBuckets.cbegin() - 1;
						Read &read = reads[read_num];
						const size_t read_len = read.seq.size();

						uint32_t ref_tile_start = (hit + cfg.first_tile_size < chr_end) ? hit : ((chr_end > cfg.first_tile_size) ? chr_end - cfg.first_tile_size : 0);
						uint32_t query_tile_start = (offset + cfg.first_tile_size < read_len) ? offset : ((read_len > cfg.first_tile_size) ? read_len - cfg.first_tile_size : 0);
						uint32_t query_tile_size = std::min(size_t(cfg.first_tile_size), read_len);

						uint32_t ovl = (offset + (chr_end - hit));

						if (ovl > cfg.min_overlap / 2) {
							ExtendLocations loc;
							loc.read_num = read_num;
							loc.chr_id = chr_id;
							loc.score = result.score;
							loc.reference_pos = ref_tile_start + result.ref_max_pos;
							loc.query_pos = query_tile_start + result.query_max_pos;
                            loc.left_hit_offsets.assign(rc_anchors[index].left_chained_hits.begin(), rc_anchors[index].left_chained_hits.end());
							loc.right_hit_offsets.assign(rc_anchors[index].right_chained_hits.begin(), rc_anchors[index].right_chained_hits.end());
							read_extend_locations.push_back(loc);
						}
												//					rc_extend_locations.push(loc);
						num_extend_requests += 1;
					}
				}
			}
		}

		slopeFilter(read_extend_locations, output.rcLocations);

		return extender_input(extender_payload(reads, output), token);
	}

void filter_body::slopeFilter(std::deque<ExtendLocations> &read_extend_locations, std::vector<ExtendLocations> &extend_locations)
{
    extend_locations.clear();
    if (!read_extend_locations.empty())
    {
        std::sort(read_extend_locations.begin(), read_extend_locations.end(),
                [](const ExtendLocations &a, const ExtendLocations &b) {
                return ((a.read_num < b.read_num) || ((a.read_num == b.read_num) && (a.score > b.score)) || ((a.read_num == b.read_num) && (a.score == b.score) && (a.reference_pos < b.reference_pos)) || ((a.read_num == b.read_num) && (a.score == b.score) && (a.reference_pos == b.reference_pos) && (a.query_pos < b.query_pos)));
                });

        for (auto l1 = read_extend_locations.begin(); l1 != read_extend_locations.end(); l1++)
        {
            if (l1->read_num == -1)
                continue;

            extend_locations.push_back(*l1);

            for (auto l2 = l1 + 1; l2 != read_extend_locations.end(); l2++)
            {
                if (l2->read_num == -1)
                    continue;

                if (l2->read_num != l1->read_num)
                    break;

                /*uint32_t rs1 = (l1->left_hit_offsets.front() >> 32);
                uint32_t qs1 = ((l1->left_hit_offsets.front() << 32) >> 32);
                uint32_t re1 = (l1->right_hit_offsets.front() >> 32);
                uint32_t qe1 = ((l1->right_hit_offsets.front() << 32) >> 32);
                uint32_t r2 = l2->reference_pos;
                uint32_t q2 = l2->query_pos;*/




                /*if ((r2 >= rs1) && (r2 <= re1) && (q2 >= qs1) && (q2 <= qe1)) {
                    l2->read_num = -1;
                    num_slope_filtered += 1;
                }
				else {*/
					float r1 = l1->reference_pos;
					float q1 = l1->query_pos;
					float r2 = l2->reference_pos;
					float q2 = l2->query_pos;
					float slope = std::abs((r1 - r2) / (q1 - q2) - 1);
					if (std::abs((r1 - r2) / (q1 - q2) - 1) <= cfg.slope_threshold)
					{
						l2->read_num = -1;
						num_slope_filtered += 1;
					}

				//}



            }
        }
    }
}
