#include "graph.h"
graph::graph(
		const char *fw_beg_file,
		const char *fw_csr_file,
		const char *bw_beg_file,
		const char *bw_csr_file)
{
	double tm=wtime();
	
	typedef int index_tt;
	typedef int vertex_tt;

	vert_count=fsize(fw_beg_file)/sizeof(index_tt) - 1;
	edge_count=fsize(fw_csr_file)/sizeof(vertex_tt);
    
    //fw
	FILE *file=fopen(fw_beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<fw_beg_file<<" cannot open\n";
		exit(-1);
	}

	index_tt *tmp_beg_pos = new index_tt[vert_count+1];
	index_tt ret=fread(tmp_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(fw_csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<fw_csr_file<<" cannot open\n";
		exit(-1);
	}

	vertex_tt *tmp_csr = new vertex_tt[edge_count];
	ret=fread(tmp_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);

    //converting to uint32_t
	fw_beg_pos = new index_t[vert_count+1];
	fw_csr = new vertex_t[edge_count];
	
	for(index_t i=0;i<vert_count+1;++i)
		fw_beg_pos[i]=(index_t)tmp_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		fw_csr[i]=(vertex_t)tmp_csr[i];

    // bw

	file=fopen(bw_beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<bw_beg_file<<" cannot open\n";
		exit(-1);
	}

	tmp_beg_pos = new index_tt[vert_count+1];
	ret=fread(tmp_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(bw_csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<bw_csr_file<<" cannot open\n";
		exit(-1);
	}

	tmp_csr = new vertex_tt[edge_count];
	ret=fread(tmp_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);

    //converting to uint32_t
	bw_beg_pos = new index_t[vert_count+1];
	bw_csr = new vertex_t[edge_count];
	
	for(index_t i=0;i<vert_count+1;++i)
		bw_beg_pos[i]=(index_t)tmp_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		bw_csr[i]=(vertex_t)tmp_csr[i];

	delete[] tmp_beg_pos;
	delete[] tmp_csr;

	std::cout<<"Graph load (success): "<<vert_count<<" verts, "
		<<edge_count<<" edges "<<wtime()-tm<<" second(s)\n";
}

graph::graph(
		const char *fw_beg_file,
		const char *fw_csr_file,
		const char *bw_beg_file,
		const char *bw_csr_file,
		const char *sub_fw_beg_file,
		const char *sub_fw_csr_file,
		const char *sub_bw_beg_file,
		const char *sub_bw_csr_file,
		const char *sub_graph_file)
{
	double tm=wtime();
	
	typedef int index_tt;
	typedef int vertex_tt;

	vert_count=fsize(fw_beg_file)/sizeof(index_tt) - 1;
	edge_count=fsize(fw_csr_file)/sizeof(vertex_tt);
    
    //fw
	FILE *file=fopen(fw_beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<fw_beg_file<<" cannot open\n";
		exit(-1);
	}

	index_tt *tmp_beg_pos = new index_tt[vert_count+1];
	index_tt ret=fread(tmp_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(fw_csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<fw_csr_file<<" cannot open\n";
		exit(-1);
	}

	vertex_tt *tmp_csr = new vertex_tt[edge_count];
	ret=fread(tmp_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);

    //converting to uint32_t
	fw_beg_pos = new index_t[vert_count+1];
	fw_csr = new vertex_t[edge_count];
	
	for(index_t i=0;i<vert_count+1;++i)
		fw_beg_pos[i]=(index_t)tmp_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		fw_csr[i]=(vertex_t)tmp_csr[i];

    // bw

	file=fopen(bw_beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<bw_beg_file<<" cannot open\n";
		exit(-1);
	}

	tmp_beg_pos = new index_tt[vert_count+1];
	ret=fread(tmp_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(bw_csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<bw_csr_file<<" cannot open\n";
		exit(-1);
	}

	tmp_csr = new vertex_tt[edge_count];
	ret=fread(tmp_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);

    //converting to uint32_t
	bw_beg_pos = new index_t[vert_count+1];
	bw_csr = new vertex_t[edge_count];
	
	for(index_t i=0;i<vert_count+1;++i)
		bw_beg_pos[i]=(index_t)tmp_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		bw_csr[i]=(vertex_t)tmp_csr[i];

	delete[] tmp_beg_pos;
	delete[] tmp_csr;

//----------------------------------------------
	file=fopen(sub_fw_beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<sub_fw_beg_file<<" cannot open\n";
		exit(-1);
	}

	index_tt *tmp_sub_fw_beg_pos = new index_tt[vert_count+1];
	ret=fread(tmp_sub_fw_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(sub_fw_csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<sub_fw_csr_file<<" cannot open\n";
		exit(-1);
	}


	vertex_tt *tmp_sub_fw_csr = new vertex_tt[edge_count];
	ret=fread(tmp_sub_fw_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);

	sub_fw_beg_pos = new index_t[vert_count+1];
	sub_fw_csr = new vertex_t[edge_count];
	
	for(index_t i=0;i<vert_count+1;++i)
		sub_fw_beg_pos[i]=(index_t)tmp_sub_fw_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		sub_fw_csr[i]=(vertex_t)tmp_sub_fw_csr[i];

	delete[] tmp_sub_fw_beg_pos;
	delete[] tmp_sub_fw_csr;

//------------------------------------------------
	file=fopen(sub_bw_beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<sub_bw_beg_file<<" cannot open\n";
		exit(-1);
	}

	index_tt *tmp_sub_bw_beg_pos = new index_tt[vert_count+1];
	ret=fread(tmp_sub_bw_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(sub_bw_csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<sub_bw_csr_file<<" cannot open\n";
		exit(-1);
	}


	vertex_tt *tmp_sub_bw_csr = new vertex_tt[edge_count];
	ret=fread(tmp_sub_bw_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);

	sub_bw_beg_pos = new index_t[vert_count+1];
	sub_bw_csr = new vertex_t[edge_count];

	
	for(index_t i=0;i<vert_count+1;++i)
		sub_bw_beg_pos[i]=(index_t)tmp_sub_bw_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		sub_bw_csr[i]=(vertex_t)tmp_sub_bw_csr[i];

	delete[] tmp_sub_bw_csr;
	delete[] tmp_sub_bw_beg_pos;
	//------------------------------------------------------
	file=fopen(sub_graph_file, "rb");
	if(file==NULL)
	{
		std::cout<<sub_graph_file<<" cannot open\n";
		exit(-1);
	}

	index_tt *tmp_sub_graph = new index_tt[vert_count+1];
	ret=fread(tmp_sub_graph, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	sub_graph = new index_t[vert_count+1];
	
	for(index_t i=0;i<vert_count+1;++i)
		sub_graph[i]=(index_t)tmp_sub_graph[i];


	delete[] tmp_sub_graph;

	std::cout<<"Graph load (success): "<<vert_count<<" verts, "
		<<edge_count<<" edges "<<wtime()-tm<<" second(s)\n";
}
