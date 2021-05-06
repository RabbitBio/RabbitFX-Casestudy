/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date:   Feb, 2020
 * This program is part of the Ktrim package
**/

#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <zlib.h>
#include "common.h"
#include "util.h"
#include "io/FastxIO.h"
#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include "io/DataQueue.h"
#include "io/Reference.h"
//using namespace std;
#include <sys/time.h>

double get_time(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}


void inline CSEREAD_resize( CSEREAD * cr, int n ) {
	cr->seq[ n] = 0;
	cr->qual[n] = 0;
	cr->size    = n;
}

void find_seed( vector<unsigned int> &seed, CSEREAD *read, const ktrim_param & kp ) {
	seed.clear();
	register char *poffset  = read->seq;
	register char *indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index1 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset );
		indexloc ++;
	}
	poffset  = read->seq + OFFSET_INDEX3;
	indexloc = poffset;
	while( true ) {
		indexloc = strstr( indexloc, kp.adapter_index3 );
		if( indexloc == NULL )
			break;
		seed.push_back( indexloc - poffset );
		indexloc ++;
	}
	sort( seed.begin(), seed.end() );
}

void workingThread_SE_C( unsigned int tn, unsigned int start, unsigned int end, CSEREAD *workingReads,
							ktrim_stat * kstat, writeBuffer * writebuffer, const ktrim_param & kp ) {

//	fprintf( stderr, "=== working thread %d: %d - %d\n", tn, start, end ), "\n";

	writebuffer->b1stored[tn] = 0;

	register int i, j;
	register unsigned int last_seed;
	vector<unsigned int> seed;
	vector<unsigned int> :: iterator it;
	const char *p, *q;

	register CSEREAD * wkr = workingReads + start;
	for( unsigned int ii=start; ii!=end; ++ii, ++wkr ) {
		//fprintf( stderr, "working: %d, %s\n", ii, wkr->id );
		//std::cout << "working : "<< ii << ", " << wkr->id << std::endl;
		// quality control
		p = wkr->qual;
		j = wkr->size;
		//std::cout << "working info:\n " << std::string(wkr->seq, wkr->size) << "\n"
		//	<< std::string(wkr->qual, wkr->size) << std::endl ;
			
		//std::cout << "p is " << p << ",  j is " << j << std::endl;
		// update in v1.2: support window check
		i = get_quality_trim_cycle_se( p, j, kp );
		//std::cout<< "before trim :" << wkr->id << " i  is " << i << "\n";

		if( i == 0 ) { // not long enough
			++ kstat->dropped[ tn ];
			continue;
		}
		if( i != j ) {  // quality-trim occurs
			CSEREAD_resize( wkr, i);
		}

		// looking for seed target, 1 mismatch is allowed for these 2 seeds
		// which means seq1 and seq2 at least should take 1 perfect seed match
		find_seed( seed, wkr, kp );

		last_seed = impossible_seed;	// a position which cannot be in seed
		for( it=seed.begin(); it!=seed.end(); ++it ) {
			if( *it != last_seed ) {
//				fprintf( stderr, " check seed: %d\n", *it );
			// as there maybe the same value in seq1_seed and seq2_seed,
			// use this to avoid re-calculate that pos
				if( check_mismatch_dynamic_SE_C( wkr, *it, kp ) )
					break;
		
				last_seed = *it;
			}
		}
		if( it != seed.end() ) {	// adapter found
			++ kstat->real_adapter[tn];
			if( *it >= kp.min_length )	{
				CSEREAD_resize( wkr, *it );
			} else {	// drop this read as its length is not enough
				++ kstat->dropped[tn];
				continue;
			}
		} else {	// seed not found, now check the tail 2, if perfect match, drop these 2; Single-end reads do not check tail 1
			i = wkr->size - 2;
			p = wkr->seq;
			if( p[i]==kp.adapter_r1[0] && p[i+1]==kp.adapter_r1[1] ) {
				++ kstat->tail_adapter[tn];
				if( i < kp.min_length ) {
					++ kstat->dropped[tn];
					continue;
				}
				CSEREAD_resize( wkr, i );
			}
		}
		writebuffer->b1stored[tn] += sprintf( writebuffer->buffer1[tn]+writebuffer->b1stored[tn],
												"%s\n%s\n+\n%s\n", wkr->id, wkr->seq, wkr->qual);
		//std::cout<< "after trim :" << wkr->id << "\n";
	}
}

int producer_fastq_task(std::string file, mash::fq::FastqDataPool* fastqPool, FqChunkQueue &dq){
    mash::fq::FastqFileReader *fqFileReader;
    fqFileReader = new mash::fq::FastqFileReader(file, *fastqPool);
    int n_chunks = 0;
    int line_sum = 0;
		double pstart = get_time();
    while(true){
        mash::fq::FastqChunk *fqchunk = new mash::fq::FastqChunk;
        fqchunk->chunk = fqFileReader->readNextChunk();
        if (fqchunk->chunk == NULL) break;
        n_chunks++;
        //std::cout << "readed chunk: " << n_chunks << std::endl;
        dq.Push(n_chunks, fqchunk->chunk);
    }
    dq.SetCompleted();
		double pend = get_time();
    //std::cout << "file " << file << " has " << n_chunks << " chunks" << 
	//					 "use time: " << pend - pstart << std::endl;
    return 0;
}
int myChunkFormat(mash::fq::FastqChunk* &fqChunk, std::vector<CSEREAD> &data, bool mHasQuality = true){
	mash::fq::FastqDataChunk * chunk = fqChunk->chunk;
	uint64_t seq_count = 0;
	uint64_t pos_ = 0;
	neoReference ref;
	CSEREAD read;
	while(true){
		ref.base = chunk->data.Pointer();
		ref.pname = pos_;
		if(mash::fq::neoGetLine(chunk, pos_, ref.lname)){
			ref.pseq = pos_; 
		} 
		else{ break;}
		mash::fq::neoGetLine(chunk, pos_, ref.lseq); 
		ref.pstrand = pos_; 
		mash::fq::neoGetLine(chunk, pos_, ref.lstrand); 
		ref.pqual = pos_;  
		mash::fq::neoGetLine(chunk, pos_, ref.lqual);
		seq_count++;
		//std::cout << "info: " << std::string((char*)ref.base + ref.pname, ref.lname) << " \n"
		//		<< std::string((char*)ref.base + ref.pseq , ref.lseq)  << "\n"
		//		<< std::string((char*)ref.base + ref.pqual, ref.lqual) << std::endl;
		//print_read(ref);
		//memcpy(read.id, (char*)ref.base + ref.pname, ref.lname);
		//memcpy(read.seq, (char*)ref.base + ref.pseq, ref.lseq);
		//memcpy(read.qual, (char*)ref.base + ref.pqual, ref.lqual);
		read.id   = (char*)ref.base + ref.pname;
		read.seq  = (char*)ref.base + ref.pseq;
		read.qual = (char*)ref.base + ref.pqual;

		register unsigned int i = ref.lseq;
		read.size = i;
		read.id[ref.lname] = 0;
		read.seq[i] = 0;
		read.qual[i] = 0;
		data.emplace_back(read);
		//std::cout<< read.id <<", lseq is "<<ref.lseq<<" size: "<<read.size<<std::endl;
	}

	return seq_count;
}

void consumer_fastq_task(mash::fq::FastqDataPool* fastqPool, FqChunkQueue &dq, writeBufferQueue &dq2, ktrim_param* kp, ktrim_stat &kstat){
	//format 
	long line_sum = 0;
    mash::int64 id = 0;
    mash::fq::FastqChunk *fqchunk = new mash::fq::FastqChunk;



  	mash::int64 wb_id = 0;
  	while(dq.Pop(id, fqchunk->chunk)){
		std::vector<CSEREAD> data;
		data.reserve(10000);
		int loaded = myChunkFormat(fqchunk, data, true);
		kstat.reads[0] += loaded;

		writeBuffer *wb = new writeBuffer;
		wb->buffer1  = new char * [ 1 ];
		wb->b1stored = new unsigned int [ 1 ];	
		wb->buffer1[0] = new char[ BUFFER_SIZE_PER_BATCH_READ_ST ];

		register unsigned int last_seed;
		vector<unsigned int> seed;
		vector<unsigned int> :: iterator it;
		//for(CSEREAD& d:data){
		//	std::cout << d.id <<", size : "<< d.size << ":\n" << std::endl;
		//	std::cout << "......................." << std::endl;
		//}
		workingThread_SE_C(0, 0, loaded, data.data(), &kstat, wb, *kp);
		//std::cout << "-----------------------------------" << std::endl;
		//std::cout << wb->buffer1[0] << std::endl;
		dq2.Push(wb_id, wb);
		wb_id++;

		fastqPool->Release(fqchunk->chunk);
	}
	dq2.SetCompleted();
	//std::cout << "----one consumer finished------" << std::endl;
}

void writer_fastq_task(writeBufferQueue& dq, ktrim_param *kp, ktrim_stat* kstats){
	string fileName = kp->outpre;
	fileName += ".read1.fq";
	FILE* fout1 = fopen(fileName.c_str(), "wt");
	writeBuffer *writebuffer;
	mash::int64 id = 0;
	double wstart = get_time();
	while(dq.Pop(id, writebuffer)){
		fwrite(writebuffer->buffer1[0], sizeof(char), writebuffer->b1stored[0], fout1);

		delete writebuffer->buffer1[0];
		delete writebuffer->b1stored;
		delete writebuffer;
	}
	//double before_close = get_time();
	fclose(fout1);
	//double wend = get_time();
	//std::cout << "write time before close: " << before_close - wstart << std::endl;
	//std::cout << "write time all: " << wend - wstart << std::endl;
	// write trim.log
	fileName = kp->outpre;
	fileName += ".trim.log";
	ofstream fout( fileName.c_str() );
	if( fout.fail() ) { 
		fprintf( stderr, "\033[1;34mError: cannot write log file!\033[0m\n" );
		//return 105;
	}
	unsigned int dropped_all=0, real_all=0, tail_all=0;	
	unsigned int line = 0;
	for(int i = 0; i < kp->thread; i++){
		dropped_all += kstats[i].dropped[0];
		real_all += kstats[i].real_adapter[0];
		tail_all += kstats[i].tail_adapter[0];
		line += kstats[i].reads[0];
	}

	fout << "Total\t"    << line			<< '\n'
		 << "Dropped\t"  << dropped_all		<< '\n'
		 << "Aadaptor\t" << real_all	<< '\n'
		 << "TailHit\t"  << tail_all << '\n';
	fout.close();

}

int process_SE_C(ktrim_param *kp){
	mash::fq::FastqDataPool *fastqPool = new mash::fq::FastqDataPool(32, 1<<22);
	FqChunkQueue queue1(64, 1);  //  because 1 producer;
	writeBufferQueue queue2(128, kp->thread);

	std::thread producer(producer_fastq_task, kp->FASTQ1, fastqPool, std::ref(queue1));
	std::thread** threads = new std::thread*[kp->thread];
	ktrim_stat *kstats = new ktrim_stat[kp->thread]; 
	for(int t = 0; t < kp->thread; t++){
		kstats[t].dropped	   = new unsigned int [ 1 ];
		kstats[t].real_adapter = new unsigned int [ 1 ];
		kstats[t].tail_adapter = new unsigned int [ 1 ];
		kstats[t].reads = new unsigned int [ 1 ];
		kstats[t].dropped[0] = 0;
		kstats[t].real_adapter[0] = 0;
		kstats[t].tail_adapter[0] = 0;
		kstats[t].reads[0] = 0;

	}
	for(int t = 0; t < kp->thread; t++){
		threads[t] = new std::thread(std::bind(consumer_fastq_task, fastqPool, std::ref(queue1), std::ref(queue2), kp, kstats[t]));
	}
	std::thread writer(writer_fastq_task, std::ref(queue2), kp, kstats);

	producer.join();
	for(int t = 0; t < kp->thread; t++){
		threads[t]->join();
	}
	writer.join();
	return 0;
}
