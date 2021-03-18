#ifndef TRIM_1_GFQ_H
#define TRIM_1_GFQ_H

#include "scc_common.h"
#include "util.h"
#include "wtime.h"
#include <omp.h>

inline void trim_1_first(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end
        )
{
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id] == 0)
        {
            scc_id[vert_id] = -1;
            continue;
        }
        else
        {
            if(bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id] == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
        }
    }
}

inline void trim_1_from_sub(
    index_t *scc_id,
    index_t *fw_beg_pos,
    index_t *bw_beg_pos,
    index_t *sub_fw_beg_pos,
    index_t *sub_bw_beg_pos,
    index_t vert_beg,
    index_t vert_end,
    index_t *mark,
    index_t *thread_bin,
    index_t *thread_queue,
    index_t *sub_graph,
    index_t tid,
    depth_t *fw_sa
)
{
    vertex_t thread_queue_frontier=0;
    index_t thread_bin_size=0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(sub_graph[vert_id]!=1)
        {
            continue;
        }

        if(fw_beg_pos[vert_id+1]-fw_beg_pos[vert_id]==0)
        {
            scc_id[vert_id]=-1;
            mark[vert_id]=2;
            fw_sa[vert_id]=0;
            thread_bin_size++;
            thread_queue[thread_queue_frontier]=vert_id;
            ++thread_queue_frontier;
            continue;
        }

        //in sub graph triming
        if(sub_fw_beg_pos[vert_id+1] - sub_fw_beg_pos[vert_id] == 0)
        {
            scc_id[vert_id] = -1;
            if(MYTEST)
            {
                std::cout<<"out degree=0::"<<vert_id<<std::endl;
            }
            continue;
        }
        else
        {
            if(sub_bw_beg_pos[vert_id+1] - sub_bw_beg_pos[vert_id] == 0)
            {
                scc_id[vert_id] = -1;
                if(MYTEST)
                {
                    std::cout<<"in degree=0::"<<vert_id<<std::endl;
                }
                continue;
            }
        }
    }
    thread_bin[tid]=thread_bin_size;
}

inline void marking_from(
    index_t *scc_id,
    index_t *bw_beg_pos,
    index_t *sub_bw_beg_pos,
    index_t vert_beg,
    index_t vert_end,
    index_t *sub_graph,
    index_t *mark
)
{

}

inline void trim_1_first_gfq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid
        )
{
    //step 1: get thread_bin_size during trim size-1
    index_t thread_bin_size = 0;
    
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id] == 0)
        {
            scc_id[vert_id] = -1;
            continue;
        }
        else
            if(bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id] == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            else
            {
                thread_bin_size ++;
            }
    }
    
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }
    
    #pragma omp barrier

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
}

//step 1.2: trim size_1
inline void trim_1_normal(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr
        )
{
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;

            }
        }
    }
}

inline void trim_1_from_fq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        index_t *small_queue
        )
{
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;

            }
        }
    }
}

// return the number of trimmed vertices
inline void trim_1_normal_only_size(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        const index_t thread_count,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid
        )
{
    //step 1: get thread_bin_size during trim size-1
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                thread_bin_size ++;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;
                thread_bin_size ++;
                continue;
            }
            
        }
    }
    
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

//    //step 3: write the vertices into fq
//    vertex_t start_pos = prefix_sum[tid];
//    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
//    {
//        if(scc_id[vert_id] == 0)
//        {
//            frontier_queue[start_pos++] = vert_id;
//        }
//    }
}

//step 1.2: trim size_1
inline void trim_1_normal_gfq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid
        )
{
    //step 1: get thread_bin_size during trim size-1
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            
            thread_bin_size ++;
        }
    }
    
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
}

inline void trim_1_from_fq_gfq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            thread_bin_size++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
    #pragma omp barrier
    if(DEBUG)
    {
        if(tid == 0)
        {
            printf("In normal trim, thread bin size, %d\n", prefix_sum[55]);
        }
    }
    //step 4: write back to small_queue
    for(index_t i=prefix_sum[tid]; i<prefix_sum[tid]+thread_bin[tid]; ++i)
    {
        frontier_queue[i] = temp_queue[i];
    }
}

inline static void get_queue(
        vertex_t *thread_queue,
        vertex_t *thread_bin,
        index_t *prefix_sum,
        index_t tid,
        vertex_t *temp_queue,
        vertex_t *fw_sa
        )
{
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }
    #pragma omp barrier
    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = start_pos; vert_id < start_pos + thread_bin[tid]; ++vert_id)
    {
        temp_queue[vert_id] = thread_queue[vert_id-start_pos];
    }
}

inline void generate_temp_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *mark,
        index_t *temp_queue,
        index_t *fw_sa
){
    index_t thread_bin_size=0;
    for(vertex_t vert_id=vert_beg;vert_id<vert_end;++vert_id)
    {
        if(mark[vert_id]==1)
        {
            thread_bin_size++;
        }
    }
    thread_bin[tid]=thread_bin_size;

    prefix_sum[tid]=0;
    #pragma omp barrier
    for(index_t i=0;i<tid;++i)
    {
        prefix_sum[tid]+=thread_bin[i];
    }

    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(mark[vert_id] == 1)
        {
            temp_queue[start_pos++] = vert_id; //此时的scc=0的点
            fw_sa[vert_id]=0;
            mark[vert_id]=0;
        }
    }
    
}


// Using prefix sum to generate frontier queue 
inline static void generate_frontier_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *sub_graph
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(sub_graph[vert_id]==1 && scc_id[vert_id] == 0)
        {
            thread_bin_size++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }


    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(sub_graph[vert_id]==1 && scc_id[vert_id] == 0)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
//    #pragma omp barrier
//    return prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

}

// Using prefix sum to generate frontier queue 
inline static void generate_not_marked(
        const index_t vert_count,
        index_t *mark,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *sub_graph,
        index_t *not_marked
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(sub_graph[vert_id]==1 && mark[vert_id] == 0)
        {
            thread_bin_size++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }


    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(sub_graph[vert_id]==1 && mark[vert_id] == 0)
        {
            not_marked[start_pos++] = vert_id;
        }
    }
}

// Using prefix sum to generate frontier queue 
inline static void generate_sub_not_marked(
        const index_t vert_count,
        index_t *mark,
        const index_t thread_count,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *sub_graph,
        index_t *thread_queue,
        index_t *not_marked
        )
{
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        vertex_t v=not_marked[vert_id];
        if(mark[v] == 0)
        {
            thread_queue[thread_bin_size++]=v;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    vertex_t idx=0;
    for(vertex_t vert_id = 0; vert_id < thread_bin_size; ++vert_id)
    {        
        not_marked[start_pos++] = thread_queue[vert_id];
    }
}


// Using prefix sum to generate frontier queue 
inline static void generate_sub_frontier_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *sub_graph,
        index_t *thread_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        vertex_t v=frontier_queue[vert_id];
        if(sub_graph[v]==1 && scc_id[v] == 0)
        {
            thread_queue[thread_bin_size++]=v;
            if(MYTEST)
            {
                std::cout<<"thread_bin: "<<vert_id<<std::endl;
            }
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    vertex_t idx=0;
    for(vertex_t vert_id = 0; vert_id < thread_bin_size; ++vert_id)
    {        
        frontier_queue[start_pos++] = thread_queue[vert_id];
    }
//    #pragma omp barrier
//    return prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

}

// Using prefix sum to generate frontier queue 
inline static void generate_sub_marking_frontier_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *sub_graph,
        index_t *mark
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(sub_graph[vert_id]==1 && mark[vert_id] == 0)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(sub_graph[vert_id]==1 && mark[vert_id] == 0)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
//    #pragma omp barrier
//    return prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

}



// Using prefix sum to generate frontier queue 
inline static void generate_marking_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *sub_graph,
        index_t *mark
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(mark[vert_id]==1)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(mark[vert_id]==1)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
//    #pragma omp barrier
//    return prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

}





inline static void gfq_from_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
    #pragma omp barrier
    //step 4: write back to small_queue
    for(index_t i=prefix_sum[tid]; i<prefix_sum[tid]+thread_bin[tid]; ++i)
    {
        small_queue[i] = temp_queue[i];
    }

}

inline static void bw_gfq_from_fw(
        index_t *fw_sa,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(fw_sa[vert_id] != -1)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    #pragma omp barrier
    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(fw_sa[vert_id] != -1)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
//    #pragma omp barrier
//    //step 4: write back to small_queue
//    for(index_t i=prefix_sum[tid]; i<prefix_sum[tid]+thread_bin[tid]; ++i)
//    {
//        small_queue[i] = temp_queue[i];
//    }

}

inline static void gfq_fw_bw_from_queue(
        index_t *sa,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(sa[vert_id] == -1)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(sa[vert_id] == -1)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
}

inline void makring(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *fw_sa,//初始都是-1，并行前创建的
        index_t *vertex_cur,
        index_t *vertex_front,
        index_t tid,
        index_t thread_count,
        const int alpha,
        const int beta,
        const int gamma,
        vertex_t *frontier_queue,
        vertex_t fq_size,
        const double avg_degree,
        vertex_t vertex_visited,
        vertex_t *temp_queue,
        index_t *prefix_sum,
        vertex_t upper_bound,
        vertex_t *thread_queue,
        vertex_t queue_size,
        index_t *mark
)
{
    depth_t level = 0;
    bool is_top_down = true;
    bool is_top_down_queue = false;

    #pragma omp barrier
    while(true)
    {
        double ltm= wtime();
        vertex_t vertex_frontier = 0;
        
//        printf("tid, %d, %d, %d, %d\n", tid, step, queue_beg, queue_end);
        #pragma omp barrier
        if(is_top_down)
        {
            vertex_t step = queue_size / thread_count;
            vertex_t queue_beg = tid * step;
            vertex_t queue_end = (tid == thread_count - 1 ? queue_size: queue_beg + step);
            for(vertex_t q_vert_id=queue_beg; q_vert_id<queue_end; q_vert_id++)
            {
                vertex_t vert_id = temp_queue[q_vert_id];//tmp就是BFS队列
                //in fq, scc_id[vert_id] is always not 0
                if(mark[vert_id]!=0 && fw_sa[vert_id] != -1)
                {
                    index_t my_beg = fw_beg_pos[vert_id];
                    index_t my_end = fw_beg_pos[vert_id+1];
                    for(; my_beg<my_end; my_beg++)//邻居
                    {
                        vertex_t nebr=fw_csr[my_beg];
                        if(mark[nebr] == 0 && fw_sa[nebr] == -1)
                        {
                            if(MYTEST)
                            {
                                std::cout<<"up down from"<<vert_id<<"to "<<nebr<<std::endl;
                            }
                            fw_sa[nebr] = level+1;
                            mark[nebr]=2;
                            scc_id[nebr]=-1;
                            thread_queue[vertex_frontier] = nebr;//邻居入的是thread_queue。thread_queue拼成temp_queue
                            vertex_frontier++;
                        }
                    }
                }
            }
        }
        else
            if(!is_top_down_queue)
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(mark[vert_id] == 0 && fw_sa[vert_id] == -1)
                    {
                        index_t my_beg = bw_beg_pos[vert_id];
                        index_t my_end = bw_beg_pos[vert_id+1];
//                        my_work_curr+=my_end-my_beg;

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=bw_csr[my_beg];
                            if(mark[nebr] != 0 && fw_sa[nebr] != -1)
//                            if(scc_id[vert_id] == 0 && fw_sa[nebr] == level)
                            {
                                if(MYTEST)
                                {
                                std::cout<<"bottom up from"<<vert_id<<"to "<<nebr<<std::endl;
                                }
                                fw_sa[vert_id] = level+1;
                                mark[vert_id]=2;
                                scc_id[vert_id]=-1;
                                vertex_frontier++;
//                                front_count++;
                                break;
                            }
                        }
                    }
                }
            }
            else
            {
                vertex_t end_queue = upper_bound;
                index_t head = 0;
                index_t tail = 0;
                //std::queue<index_t> q;
                vertex_t step = queue_size / thread_count;
                vertex_t queue_beg = tid * step;
                vertex_t queue_end = (tid == thread_count - 1 ? queue_size: queue_beg + step);

                //Option 1: put current level vertices into fq
                for(vertex_t q_vert_id=queue_beg; q_vert_id<queue_end; q_vert_id++)
                {
                    thread_queue[tail] = temp_queue[q_vert_id]; //temp是level队列
                    tail ++;
                }
                while(head != tail)//对thread——queue中的元素进行BFS
                {
                    vertex_t temp_v = thread_queue[head++];
//                    front_count ++;
                    if(head == end_queue)
                        head = 0;
                    index_t my_beg = fw_beg_pos[temp_v];
                    index_t my_end = fw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        
                        if(mark[w] == 0 && fw_sa[w] == -1)
                        {
                            std::cout<<"last BFS from"<<temp_v<<"to"<<w<<std::endl;
                            thread_queue[tail++] = w;
                            if(tail == end_queue)
                                tail = 0;
                            fw_sa[w] = level + 1;
                            mark[w]=2;
                            scc_id[w]=-1;
                        }
                    }
                }
            }

        
        vertex_front[tid] = vertex_frontier;//vertexfrontier指向的是thread queue
        #pragma omp barrier
        vertex_frontier = 0;

        for(index_t i=0; i<thread_count; ++i)
        {
            vertex_frontier += vertex_front[i];
        }
        vertex_visited += vertex_frontier;//所有线程thread queue的内元素数量和
//        #pragma omp barrier
        if(VERBOSE)
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
			if(tid==0 && level < 50) 
                std::cout<<"Level-"<<(int)level<<" "
//				<<"-frontier-time-visited:"
				<<vertex_frontier<<" "
                <<fq_size<<" "
                <<(double)(fq_size)/vertex_frontier<<" "
				<<(wtime() - ltm) * 1000<<"ms "
				<<vertex_visited<<" "
                <<edge_frontier<<" "
                <<edge_remaider<<" "
                <<edge_remaider/edge_frontier<<"\n";
        }
        
        if(vertex_frontier == 0) break;
        
        if(is_top_down) 
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;
//            printf("edge_remainder/alpha = %g, edge_froniter = %g\n", edge_remainder / alpha, edge_frontier);
            if(!is_top_down_queue && (edge_remainder / alpha) < edge_frontier)//转bottom up
            {
                is_top_down = false;
                if(VERBOSE)
                {
                    if(tid==0)
                    {
//                        double Nf = vertex_frontier;
//                        double Nu = fq_size - vertex_visited;
//                        double Mf = Nf * Nf / Nu + avg_degree * (Nu - Nf);
//                        printf("mf=%.0lf, mu=%.0lf, alpha=%d, Mf=%.0lf\n", edge_frontier, edge_remainder, ALPHA, Mf);
                        std::cout<<"--->Switch to bottom up\n";
                    }
                }
            }

        }
        else
            if((!is_top_down && !is_top_down_queue && (fq_size*1.0/beta) > vertex_frontier) || (!is_top_down && !is_top_down_queue && level > gamma))
//            if(level > 10)
            {
                //if(!is_top_down_queue)
                
                vertex_frontier = 0;
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(mark[vert_id] != 0 && fw_sa[vert_id] == level + 1)
                    {
                        thread_queue[vertex_frontier] = vert_id;
                        vertex_frontier++;
                    }
                }
                vertex_front[tid] = vertex_frontier;
                
                is_top_down = false;
                is_top_down_queue = true;

                if(VERBOSE)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to top down queue\n";
                }
            }
        
        #pragma omp barrier

        if(is_top_down || is_top_down_queue)
        {
            get_queue(thread_queue,//合并thread-queue放到temp里
                    vertex_front,
                    prefix_sum,
                    tid,
                    temp_queue,
                    fw_sa);
            queue_size = prefix_sum[thread_count-1] + vertex_front[thread_count-1];
//            if(VERBOSE)
//            {
//                if(tid == 0)
//                    printf("queue_size, %d\n", queue_size);
//            }
        }
        #pragma omp barrier
        
        level ++;
    }
}

#endif


