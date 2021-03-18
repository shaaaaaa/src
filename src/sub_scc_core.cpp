#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <omp.h>
#include <algorithm>
#include "wtime.h"
#include "graph.h"
//#include "frontier_queue.h"
#include "scc_common.h"
//#include "trim_1.h"
#include "trim_1_gfq.h"
#include "trim_2_3.h"
#include "color_propagation.h"
#include "fw_bw.h"
#include "openmp_wcc.hpp"
#define INF -1

//0 trim, 1 largest SCC, 2 small SCC, 3 total time
//4 trim_size_1, 5 trim_size_2, 6 pivot_selection, 7 fw_bfs, 8 bw_bfs, 9 color propagation, 10 color identify, 11 color_init 

void sub_scc_detection(
        const graph *g,
        const int alpha, 
        const int beta,
        const int gamma,
        const double theta,
        const index_t thread_count,
        double *avg_time
        )
{
    const index_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;
    if(DEBUG)
        printf("vert_count = %d, edge_count = %d, avg_degree = %.3lf\n", vert_count, edge_count, avg_degree);
    //step 0: initialization
//    return;

    index_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    index_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    index_t *sub_fw_beg_pos = g->sub_fw_beg_pos;
    vertex_t *sub_fw_csr = g->sub_fw_csr;
    index_t *sub_bw_beg_pos = g->sub_bw_beg_pos;
    vertex_t *sub_bw_csr = g->sub_bw_csr;
    index_t *sub_graph=g->sub_graph;

    if(MYTEST)
    {
        for(int i=0;i<vert_count;i++)
        {
            if(sub_graph[i]==1)
            {
                std::cout<<"i:"<<i<<"exsit"<<std::endl;
            }

            index_t my_beg = sub_fw_beg_pos[i];
            index_t my_end = sub_fw_beg_pos[i+1];
            for(; my_beg<my_end; my_beg++)//邻居
            {
                vertex_t nebr=sub_fw_csr[my_beg];
                std::cout<<i<<" to "<<nebr<<std::endl;
            }

            my_beg = sub_bw_beg_pos[i];
            my_end = sub_bw_beg_pos[i+1];
            for(; my_beg<my_end; my_beg++)//邻居
            {
                vertex_t nebr=sub_bw_csr[my_beg];
                std::cout<<i<<" from "<<nebr<<std::endl;
            }

        }
    }

    if(VERBOSE)
    {
        for(int i=fw_beg_pos[vert_count]; i<fw_beg_pos[vert_count+1]; ++i)
            printf("%d\n", fw_csr[i]);
    }

    index_t *scc_id = new index_t[vert_count + 1];
   
    index_t *color = new index_t[vert_count + 1];
//    index_t *color_times = new index_t[vert_count + 1];

    index_t *max_pivot_list = new index_t[thread_count];
    index_t *max_degree_list = new index_t[thread_count];

    index_t *not_marked=new index_t[vert_count+1];
    

//    index_t *mul_degree = new index_t[vert_count + 1];
//    index_t *degree_prop = new index_t[vert_count + 1];
	
    depth_t *fw_sa;
	depth_t *bw_sa;

//    fw_sa = new depth_t[vert_count+1];
//    bw_sa = new depth_t[vert_count+1];

	if(posix_memalign((void **)&fw_sa,getpagesize(),
		sizeof(depth_t)*(vert_count+1)))
		perror("posix_memalign");
	
    if(posix_memalign((void **)&bw_sa,getpagesize(),
		sizeof(depth_t)*(vert_count+1)))
		perror("posix_memalign");
    
    index_t *small_queue = new index_t[vert_count + 1];
    index_t *temp_queue = new index_t[vert_count + 1];
    index_t *inter_queue = new index_t[vert_count + 1];
    index_t *wcc_fq= new index_t[vert_count + 1];
    
    index_t *thread_bin = new index_t[thread_count];
    index_t *prefix_sum = new index_t[thread_count];
	
    index_t *front_comm=new index_t[thread_count];	
	index_t *work_comm=new index_t[thread_count];
    bool *color_change = new bool[thread_count];
    memset(color_change, 0, sizeof(bool) * thread_count);

    index_t *mark=new index_t[vert_count+1];

    //WCC + FW-BW
//    index_t *wcc_color = new index_t[vert_count + 1];
//    index_t *color_redirect = new index_t[vert_count + 1];
//    bool *is_redirect = new bool[thread_count];
//    color_t *global_color = new color_t[1];
//    global_color[0] = 0;

    vertex_t wcc_fq_size = 0;

//    memset(fw_sa, -1, sizeof(depth_t)*vert_count);
//	memset(bw_sa, -1, sizeof(depth_t)*vert_count);
//    memset(scc_id, 0, sizeof(index_t) * (vert_count + 1));

//Initialization
    
    if(DEBUG)
    {
        printf("Initialization\n");
    }

//    #pragma omp parallel for
    for(index_t i=0; i<vert_count + 1; ++i)
    {
        color[i] = i;
//            color_times[i] = 0;
        fw_sa[i] = -1;
        bw_sa[i] = -1;
        scc_id[i] = 0;
//            wcc_color[i] = -1;
//            color_redirect[i] = i;
        mark[i]=0;
    }

    //step 1: trim size-1
    //step 2: largest_scc, asynchronize, direction optimized
    //step 3: trim size-1 && size-2 && size-3
    /// generate new FQ
    //step 4: small sccs, graph coloring
    
    if(DEBUG)
    {
        printf("Parallel starts\n");
    }
    vertex_t vertex_fw = 0;
    vertex_t vertex_bw = 0;
    index_t size_3_1 = 0;
    index_t size_3_2 = 0;
    bool changed = false;
    double end_time;
    double start_time = wtime();
    #pragma omp parallel \
    num_threads(thread_count) 
//    \shared(global_color)
    {
        const index_t tid = omp_get_thread_num();
        index_t step = vert_count / thread_count;
        index_t vert_beg = tid * step;
        index_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
        double time_size_1_first;
        double time_size_1;
        double time_fw;
        double time_bw; 
        double time_size_2;
        double time_size_3;
        double time_gfq;
        double time_color_1;
        double time_color_2;
        double time_color;
        double pivot_time;
        double time_color_init;
        double time_wcc;
        double time_mice_fw_bw; 
        const vertex_t upper_bound = vert_count / thread_count * 5;
        vertex_t *thread_queue = new vertex_t[vert_count+1];
        
        
        double time = wtime();
// not using frontier queue for trim
//
        index_t trim_times = 1;

// change to control by trimmed vertices
        
        index_t queue_size=0;

        trim_1_from_sub(scc_id,
                fw_beg_pos, 
                bw_beg_pos,
                sub_fw_beg_pos,
                sub_bw_beg_pos,
                vert_beg,
                vert_end,
                mark,
                thread_bin,//数量存在thread-bin里了
                thread_queue,//结果存在thread-queue里了
                sub_graph,
                tid,
                fw_sa);
        trim_times ++;
        #pragma omp barrier
        if(tid == 0)
        {
            time_size_1_first = wtime() - time;
        }

//trim1后,初始化mark起点
        get_queue(
            thread_queue,//把这个
            thread_bin,
            prefix_sum,
            tid,
            temp_queue,//合到这里了
            fw_sa
        );
        #pragma omp barrier

        index_t temp_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

        if(MYTEST && tid==0)
        {
            for(index_t idx=0;idx<vert_count;idx++)
            {
                std::cout<<"small queue before gen::"<<small_queue[idx]<<std::endl;
            }
        }

        generate_not_marked(
            vert_count,
            mark,
            thread_count,
            small_queue,//目前所有没有确定的点
            thread_bin,
            prefix_sum,
            vert_beg,
            vert_end,
            tid,
            sub_graph,
            not_marked
        );
        #pragma omp barrier

        index_t mark_fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
        index_t mark_step = mark_fq_size / thread_count;
        index_t mark_vert_beg = tid * mark_step;
        index_t mark_vert_end = (tid == thread_count - 1 ? mark_fq_size : mark_vert_beg + mark_step);

        if(MYTEST && tid==0)
        {
            std::cout<<"mark_fq_size "<<mark_fq_size<<std::endl;
        }

        time = wtime();

        makring(
            scc_id,
            sub_bw_beg_pos,
            sub_fw_beg_pos,
            mark_vert_beg,
            mark_vert_end,
            sub_bw_csr,
            sub_fw_csr,
            fw_sa,//初始都是-1，并行前创建的
            front_comm,
            work_comm,
            tid,
            thread_count,
            alpha,
            beta,
            gamma,
            not_marked,
            mark_fq_size,
            avg_degree,
            vertex_fw,
            temp_queue,//要marking的起点
            prefix_sum,
            upper_bound,
            thread_queue,
            temp_size,
            mark
        );
        #pragma omp barrier

        time_size_1_first = wtime() - time;

        if(MYTEST && tid==0)
        {
            for(index_t idx=0;idx<vert_count;idx++)
            {
                std::cout<<"small queue before::"<<small_queue[idx]<<std::endl;
            }
        }

        generate_frontier_queue(vert_count,
            scc_id,
            thread_count,
            small_queue,//目前所有没有确定的点
            thread_bin,
            prefix_sum,
            vert_beg,
            vert_end,
            tid,
            sub_graph
        );
        #pragma omp barrier

        index_t fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
        step = fq_size / thread_count;
        vert_beg = tid * step;//这是smalll_queue的size算的
        vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);

        if(tid == 0)
        {
            time_size_1 = wtime() - time;
            printf("time_size_1, %.3lf, fq_size, %d\n", time_size_1 * 1000, fq_size);
        }

        if(MYTEST && tid==0)
        {
            for(index_t idx=0;idx<fq_size;idx++)
            {
                std::cout<<"small queue::"<<small_queue[idx]<<std::endl;
            }
        }
        
        if(DEBUG)
        {
            if(tid == 0)
            {
                printf("fq_size, %d\n", fq_size);
            }
        }
        
        time = wtime();

        

        if(fq_size > 0)
        {
     // not using trim-2
            if(DEBUG)
            {
                if(tid == 0)
                    printf("fq_size, %d\n", fq_size);
            }
            time = wtime();

            //init_fw_sa_only(vert_beg,
            //        vert_end,
            //        fw_sa);
            //#pragma omp barrier


            //step 3.2: trim size_2
    //        trim_2_from_graph(scc_id,
    //                fw_beg_pos,
    //                bw_beg_pos,
    //                vert_beg,
    //                vert_end,
    //                fw_csr,
    //                bw_csr);
            if(MYTEST && tid==0)
            {
                for(index_t idx=0;idx<vert_count;idx++)
                {
                    std::cout<<"scc_id::--->"<<scc_id[idx]<<std::endl;
                }
            }

            trim_2_from_fq(scc_id,
                    sub_fw_beg_pos,
                    sub_bw_beg_pos,
                    vert_beg,
                    vert_end,
                    sub_fw_csr,
                    sub_bw_csr,
                    small_queue,
                    thread_queue,
                    thread_bin,
                    tid,
                    fw_sa,
                    mark);
            #pragma omp barrier
            if(tid == 0)
                time_size_2 = wtime() - time;
            
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("time_size_2, %.3lf, fq_size, %d\n", time_size_2 * 1000, fq_size);
                }
            }
//note! just for record            
            //time_size_1 = 0;
            trim_times = 0;
            #pragma omp barrier
            
            
// not using trim-3
            //step 3.3: trim size-3
            //size-3: type_1, A --> B --> C --> A
            time = wtime();
    //        trim_3_1_from_graph(scc_id, 
    //                fw_beg_pos,
    //                bw_beg_pos,
    //                vert_beg,
    //                vert_end,
    //                fw_csr,
    //                bw_csr);

            trim_3_1_from_fq(scc_id,
                    sub_fw_beg_pos,
                    sub_bw_beg_pos,
                    vert_beg,
                    vert_end,
                    sub_fw_csr,
                    sub_bw_csr,
                    small_queue,
                    thread_queue,
                    thread_bin,
                    tid,
                    fw_sa,
                    mark);
            #pragma omp barrier
            double time_size_3_1 = wtime() - time;

            if(DEBUG)
            {
                if(tid == 0)
                    printf("time_size_3_1, %.3lf\n", time_size_3_1 * 1000);
            }

            //step 3.3: trim size-3
            //size-3: type_1, A --> B --> A --> C --> A 
            //starting from hub vertex A
            
            time = wtime();
    //        trim_3_2_from_graph(scc_id,
    //                fw_beg_pos,
    //                bw_beg_pos,
    //                vert_beg,
    //                vert_end,
    //                fw_csr,
    //                bw_csr);
            trim_3_2_from_fq(scc_id,
                    sub_fw_beg_pos,
                    sub_bw_beg_pos,
                    vert_beg,
                    vert_end,
                    sub_fw_csr,
                    sub_bw_csr,
                    small_queue,
                    thread_queue,
                    thread_bin,
                    tid,
                    fw_sa,
                    mark);
            
            trim_times = 0;
            
/// note!
//            time_size_1 = 0;
            #pragma omp barrier
            double time_size_3_2 = wtime() - time;

            if(DEBUG)
            {
                if(tid == 0)
                    printf("time_size_3_2, %.3lf\n", time_size_3_2 * 1000);
            }
            time_size_3 = time_size_3_1 + time_size_3_2;
            
            
            //step 3.4: trim size_1
                
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("trim_1 times before Color, %d, fq_size, %d\n", trim_times, fq_size);
                }
            }
          
            #pragma omp barrier
            

            get_queue(
                thread_queue,//前几次trim里产生的起点都在这里，合并到temp-queue
                thread_bin,
                prefix_sum,
                tid,
                temp_queue,
                fw_sa
            );
            #pragma omp barrier

            temp_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

            if(MYTEST && tid==0)
            {
                std::cout<<"temp_size"<<temp_size<<std::endl;
                for(index_t idx=0;idx<vert_count;idx++)
                {
                    std::cout<<"sssmall queue"<<small_queue[idx]<<std::endl; 
                }
            }


            generate_sub_not_marked(vert_count,//重组small queue
                mark,
                thread_count,
                thread_bin,//尚未访问点个数
                prefix_sum,
                mark_vert_beg,
                mark_vert_end,
                tid,
                sub_graph,
                thread_queue,
                not_marked
            );
            #pragma omp barrier

            mark_fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            mark_step = mark_fq_size / thread_count;
            mark_vert_beg = tid * mark_step;
            mark_vert_end = (tid == thread_count - 1 ? mark_fq_size : mark_vert_beg + mark_step);

            if(MYTEST && tid==0)
            {
                std::cout<<"fqqqq"<<fq_size<<std::endl;
                for(index_t idx=0;idx<fq_size;idx++)
                {
                    std::cout<<"temp_queue::--->"<<small_queue[idx]<<std::endl;
                }
            }

            #pragma omp barrier
            time=wtime();
            makring(
                scc_id,
                sub_bw_beg_pos,
                sub_fw_beg_pos,
                mark_vert_beg,
                mark_vert_end,
                sub_bw_csr,
                sub_fw_csr,
                fw_sa,//初始都是-1，并行前创建的
                front_comm,
                work_comm,
                tid,
                thread_count,
                alpha,
                beta,
                gamma,
                not_marked,
                mark_fq_size,
                avg_degree,
                vertex_fw,
                temp_queue,
                prefix_sum,
                upper_bound,
                thread_queue,
                temp_size,
                mark
            );
            #pragma omp barrier
            time_size_3+=wtime()-time;
            

            generate_sub_frontier_queue(vert_count,
                scc_id,
                thread_count,
                small_queue,//目前所有没有确定的点
                thread_bin,
                prefix_sum,
                vert_beg,
                vert_end,
                tid,
                sub_graph,
                thread_queue
            );
            #pragma omp barrier

            fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            step = fq_size / thread_count;
            vert_beg = tid * step;//这是smalll_queue的size算的
            vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);


            time = wtime();
// color WCC, with both FW & BW
            coloring_wcc(fq_size,
                    scc_id,
                    thread_count,
                    small_queue,
                    vert_beg,
                    vert_end,
                    tid,
                    color,
                    color_change,
                    sub_fw_beg_pos,
                    sub_fw_csr,
                    sub_bw_beg_pos,
                    sub_bw_csr);

            #pragma omp barrier
            if(tid == 0)
            {
                time_wcc = wtime() - time;
                //printf("WCC time (ms), %.3lf\n", time_wcc * 1000);
//                printf("global color, %d\n", global_color[0]);
//                for(int i=0; i<fq_size; ++i)
//                {
//                    printf("%d, %d\n", small_queue[i], wcc_color[small_queue[i]]);
//                }
            }
            if(MYTEST && tid==0)
            {
                for(int idx=0;idx<vert_count;idx++)
                {
                    std::cout<<"color:"<<idx<<"::"<<color[idx]<<std::endl;
                }
            }

// using fw-bw to detect mice SCCs
            init_fw_sa(vert_beg,
                    vert_end,
                    fw_sa,
                    small_queue,
                    fq_size,
                    wcc_fq,
                    tid,
                    color,
                    wcc_fq_size);
            #pragma omp barrier
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("wcc_time, %.3lf\n", time_wcc * 1000);
                    printf("wcc_fq, %d\n", wcc_fq_size);
                }
            }
//            if(tid == 0)
//            {

            time = wtime();
            mice_fw_bw(color,
                    scc_id,
                    sub_fw_beg_pos,
                    sub_bw_beg_pos,
                    sub_fw_csr,
                    sub_bw_csr,
                    fw_sa,
                    tid,
                    thread_count,
                    small_queue,
                    fq_size,
                    wcc_fq,
                    wcc_fq_size,
                    thread_queue,
                    thread_bin,
                    mark
            ); 
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("mice fw bw success");
                }
            }
//            }
            #pragma omp barrier
            if(tid == 0)
            {
                time_mice_fw_bw = wtime() - time;
//                printf("FW-BW mice (ms), %.3lf\n", time_mice_fw_bw * 1000);
            }

            get_queue(
                thread_queue,
                thread_bin,
                prefix_sum,
                tid,
                temp_queue,
                fw_sa
            );
            #pragma omp barrier

            index_t temp_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

            generate_sub_not_marked(vert_count,
                mark,
                thread_count,
                thread_bin,
                prefix_sum,
                mark_vert_beg,
                mark_vert_end,
                tid,
                sub_graph,
                thread_queue,
                not_marked
            );
            #pragma omp barrier

            if(tid == 0)
            {
                time_size_1 = wtime() - time;
            }
            mark_fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            mark_step = mark_fq_size / thread_count;
            mark_vert_beg = tid * step;
            mark_vert_end = (tid == thread_count - 1 ? mark_fq_size : mark_vert_beg + mark_step);

            time = wtime();

            makring(
                scc_id,
                sub_bw_beg_pos,
                sub_fw_beg_pos,
                mark_vert_beg,
                mark_vert_end,
                sub_bw_csr,
                sub_fw_csr,
                fw_sa,//初始都是-1，并行前创建的
                front_comm,
                work_comm,
                tid,
                thread_count,
                alpha,
                beta,
                gamma,
                small_queue,
                mark_fq_size,
                avg_degree,
                vertex_fw,
                temp_queue,
                prefix_sum,
                upper_bound,
                thread_queue,
                temp_size,
                mark
            );
            #pragma omp barrier

            time_mice_fw_bw+=wtime()-time;
        }
        if(tid == 0)
        {
            avg_time[0] += time_size_1_first + time_size_1 + time_size_2 + time_size_3;
//            avg_time[0] += time_size_1;
            avg_time[1] += time_fw + time_bw + pivot_time;
            avg_time[2] += time_wcc + time_mice_fw_bw;
            avg_time[4] += time_size_1_first + time_size_1;
//            avg_time[4] += time_size_1;
            avg_time[5] += time_size_2;
            avg_time[6] += pivot_time;
            avg_time[7] += time_fw;
            avg_time[8] += time_bw;
            avg_time[9] += time_wcc;
            avg_time[10] += time_mice_fw_bw;
//            avg_time[11] += time_color_init;
            
            avg_time[13] += time_size_3;
            avg_time[14] += time_gfq;
            
//            if(DEBUG)
//            {
//                for(int i=0; i< 15; ++i)
//                    printf("%.3lf\n", avg_time[i]);
//                printf("algorithm run finishes\n");
//            }
        }
        if(OUTPUT_TIME)
        {
            if(tid == 0)
            {
                printf("\ntime size_1_first, %.3lf\ntime size_1, %.3lf\ntime pivot, %.3lf\nlargest fw, %.3lf\nlargest bw, %.3lf\nlargest fw/bw, %.3lf\ntrim size_2, %.3lf\ntrim size_3, %.3lf\nwcc time, %.3lf\nmice fw-bw time, %.3lf\nmice scc time, %.3lf\ntotal time, %.3lf\n", time_size_1_first * 1000, time_size_1 * 1000, pivot_time * 1000, time_fw * 1000, time_bw * 1000, (pivot_time + time_fw + time_bw) * 1000, time_size_2 * 1000, time_size_3 * 1000, time_wcc * 1000, time_mice_fw_bw * 1000, (time_wcc + time_mice_fw_bw) * 1000, (time_size_1_first + time_size_1 + pivot_time + time_fw + time_bw + time_size_2 + time_size_3 + time_wcc + time_mice_fw_bw) * 1000);
            }
        }
        #pragma omp barrier
    }
    end_time = wtime() - start_time;
    avg_time[3] += end_time;
    if(DEBUG)
        printf("total time, %.3lf\n", end_time * 1000);

    get_scc_result(scc_id,
            vert_count);

    delete[] mark;
    delete[] not_marked;
    
    delete[] scc_id;
    delete[] color;
    delete[] max_pivot_list;
    delete[] max_degree_list;
    delete[] small_queue;
    delete[] temp_queue;
    delete[] thread_bin;
    delete[] prefix_sum;
    delete[] front_comm;
	delete[] work_comm;
    delete[] color_change;
//    delete[] wcc_color;
//    delete[] color_redirect;
//    delete[] is_redirect;

}


