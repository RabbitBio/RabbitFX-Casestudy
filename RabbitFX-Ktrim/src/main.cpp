//#include "io/FastxIO.h"
#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include <string>
#include <iostream>
#include "cmdparser/CLI11.hpp"
#include "io/DataQueue.h"
#include <thread>
#include "param_handler.h"
#include "util.h"
#include "pe_handler.h"
#include "se_handler.h"

typedef  rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> FqChunkQueue;

int count_line(rabbit::fq::FastqChunk* fqchunk){
    return 1000;
}

/*
int producer_fastq_task(std::string file, mash::fq::FastqDataPool* fastqPool, FqChunkQueue &dq){
    rabbit::fq::FastqFileReader *fqFileReader;
    //mash::fq::FastqReader *fastqReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool, false);
    //fastqReader = new mash::fq::FastqReader(*fqFileReader, *fastqPool);  //没有必要再分fastqreader和fastareader了，只要上面的filereader是不同的类型就可以了。函数重载readnextchunk和
    int n_chunks = 0;
    int line_sum = 0;
    while(true){
        mash::fq::FastqChunk *fqchunk = new mash::fq::FastqChunk;
        fqchunk->chunk = fqFileReader->readNextChunk();
        if (fqchunk->chunk == NULL) break;
        n_chunks++;
        std::cout << "readed chunk: " << n_chunks << std::endl;
        //dq.Push(n_chunks, fqchunk->chunk);
        //line_sum += count_line(fqchunk);
        //std::vector<neoReference> data;
        //data.resize(10000);
        //line_sum += mash::fq::chunkFormat(fqchunk, data, true);
        fastqPool->Release(fqchunk->chunk);
    }
    dq.SetCompleted();
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
    return 0;
}
void consumer_fastq_task(mash::fq::FastqDataPool* fastqPool, FqChunkQueue &dq){
    long line_sum = 0;
    mash::int64 id = 0;
    std::vector<Reference> data;
    mash::fq::FastqChunk *fqchunk = new mash::fq::FastqChunk;
    data.resize(10000);
    while(dq.Pop(id, fqchunk->chunk)){
        line_sum += mash::fq::chunkFormat(fqchunk, data, true);
        fastqPool->Release(fqchunk->chunk);
    }
    std::cout << "line_sum: " << line_sum << std::endl;
}

void print_chunk(mash::fa::FastaDataChunk *chunk){
    std::cout << "chunk size: " << chunk->size << std::endl;
    std::cout << "chunk head: " << std::string((char*)chunk->data.Pointer(), 100) << std::endl;
}
void print_fachunkpart_info(mash::fa::FastaChunk *fachunk){
    std::cout << "------------------chunk info-----------------" << std::endl;
    print_chunk(fachunk->chunk);
    mash::fa::FastaDataChunk *current_chunk = fachunk->chunk;
    while(current_chunk->next != NULL){
        std::cout << "next" << std::endl;
        current_chunk = current_chunk->next;
        print_chunk(current_chunk);
    }
}

int producer_fasta_task(std::string file){
    mash::fa::FastaDataPool *fastaPool = new mash::fa::FastaDataPool();
    mash::fa::FastaFileReader *faFileReader;
    //mash::fq::FastqReader *fastqReader;
    faFileReader = new mash::fa::FastaFileReader(file, *fastaPool, false);
    //fastqReader = new mash::fq::FastqReader(*fqFileReader, *fastqPool);  //没有必要再分fastqreader和fastareader了，只要上面的filereader是不同的类型就可以了。函数重载readnextchunk和
    int n_chunks = 0;
    int line_sum = 0;
    while(true){
        mash::fa::FastaChunk *fachunk = new mash::fa::FastaChunk;
        fachunk = faFileReader->readNextChunkList();
        //fachunk = faFileReader->readNextChunk();
        if (fachunk == NULL) break;
        n_chunks++;
        //line_sum += count_line(fqchunk);
        std::vector<Reference> data;
        //print_fachunkpart_info(fachunk);
        //-----relaease
        mash::fa::FastaDataChunk * tmp = fachunk->chunk;
        do{
            fastaPool->Release(tmp);
            tmp = tmp->next;
        }while(tmp != NULL);
        //------release
        //line_sum += mash::fa::chunkFormat(*fachunk, data);
    }
    std::cout << "file " << file << " has " << line_sum << " lines" << std::endl;
    return 0;

    //result record: readnextchunklist: 2.85//3.25
    //               readnextchunk:     2.76
}
*/
int main(int argc, char** argv){

    static ktrim_param kp;
    init_param( kp );
    int retValue = process_cmd_param(argc, argv, kp);

    if(retValue == 100){
        return 0;
    }else if(retValue != 0){
        return retValue;
    }
    //----------------------------------------
    if(kp.FASTQ2 == NULL) { // se process
        retValue = process_SE_C(&kp);
    } else { //pe process
        retValue = process_PE_C(&kp);
    }  


    /*
    mash::fq::FastqDataPool *fastqPool = new mash::fq::FastqDataPool(32, 1<<22);
    FqChunkQueue queue1(64, 1);  //  because 1 producer;
    std::thread producer(producer_fastq_task, filename, fastqPool, std::ref(queue1));
    std::thread** threads = new std::thread*[th];
    for(int t = 0; t < th; t++){
        threads[t] = new std::thread(std::bind(consumer_fastq_task, fastqPool, std::ref(queue1)));
    }
    producer.join();
    for(int t = 0; t < th; t++){
        threads[t]->join();
    }
    */
    return 0;
}
