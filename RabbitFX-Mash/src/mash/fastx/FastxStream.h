/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00

  last modified by Zekun Yin 2020/5/18
*/
#ifndef H_FASTQSTREAM
#define H_FASTQSTREAM

#include "../Sketch.h"

#include "Globals.h"

#include "Buffer.h"
#include "FastxChunk.h"
#include "utils.h"
#include <iostream>
#include <string>
#include <zlib.h> //support gziped files, functional but inefficient
#include "FileReader.h"


#if defined (_WIN32)
#   define _CRT_SECURE_NO_WARNINGS
#   pragma warning(disable : 4996) // D_SCL_SECURE
#   pragma warning(disable : 4244) // conversion uint64 to uint32
//# pragma warning(disable : 4267)
#   define FOPEN    fopen
#   define FDOPEN    fdopen
#   define FSEEK    _fseeki64
#   define FTELL    _ftelli64
#   define FCLOSE   fclose
#elif __APPLE__ // Apple by default suport 64 bit file operations (Darwin 10.5+)
#   define FOPEN    fopen
#   define FSEEK    fseek
#   define FTELL    ftell
#   define FCLOSE   fclose
#else
#   if !defined(_LARGEFILE_SOURCE)
#       define _LARGEFILE_SOURCE
#       if !defined(_LARGEFILE64_SOURCE)
#           define _LARGEFILE64_SOURCE
#       endif
#   endif
#   if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS != 64)
#       undef _FILE_OFFSET_BITS
#   endif
#   if !defined(_FILE_OFFSET_BITS)
#       define _FILE_OFFSET_BITS 64
#   endif
#   define FOPEN    fopen64
#   define FDOPEN   fdopen
#   define FSEEK    fseeko64
#   define FTELL    ftello64
#   define FCLOSE   fclose
#endif

namespace mash
{

namespace fa
{

class FastaFileReader
{
private:
	//static const uint32 SwapBufferSize = 1 << 20;
	static const uint32 SwapBufferSize = 1 << 26;//16MB
	//static const uint32 SwapBufferSize = 1 << 13;
  FileReader* mFaReader;

public:
	FastaFileReader(const std::string& fileName_, uint64 halo = 21, bool isZippedNew = false)
		:	swapBuffer(SwapBufferSize)
		,	bufferSize(0)
		,	eof(false)
		,	usesCrlf(false)
		,	totalSeqs(0)
		,	mHalo(halo)
		,	isZipped(isZippedNew)
	{	
		//if(ends_with(fileName_,".gz"))
		if(isZipped || ends_with(fileName_, ".gz"))
		{
      isZipped = true;
		}
    mFaReader = new FileReader(fileName_, isZipped);
			
	}

  FastaFileReader(const int fd, uint64 halo = 21, bool isZippedNew = false)
		:	swapBuffer(SwapBufferSize)
		,	bufferSize(0)
		,	eof(false)
		,	usesCrlf(false)
		,	totalSeqs(0)
		,	mHalo(halo)
		,	isZipped(isZippedNew)
	{	
		//if(ends_with(fileName_,".gz"))
    mFaReader = new FileReader(fd, isZipped);
			
	}

	~FastaFileReader()
	{
		//std::cerr << "totalSeqs: " << this->totalSeqs << std::endl;
		if(mFile != NULL || mZipFile != NULL)
			delete mFaReader;
	}

	bool Eof() const
	{
		return eof;
	}

	bool ReadNextChunk(FastaChunk* chunk_, SeqInfos& seqInfos);

	void Close()
	{
		if(mFile != NULL){
			FCLOSE(mFile);
			mFile = NULL;
		}

		if(mZipFile !=NULL && isZipped){
			gzclose(mZipFile);
			mZipFile=NULL;
		}

	}

	int64 Read_back(byte* memory_, uint64 size_)
	{	
		if(isZipped){
			int64 n = gzread(mZipFile,memory_,size_);
			if(n == -1)
				std::cerr << "Error to read gzip file" << std::endl;
			return n;
		}
		else{
			int64 n = fread(memory_, 1, size_, mFile) ;
			return n;
		}
	}

  int64 Read(byte* memory_, uint64 size_){
    return mFaReader->Read(memory_, size_);
  }

private:
	core::Buffer	swapBuffer;
	uint64			bufferSize;
	bool			eof;
	bool			usesCrlf;
	bool			isZipped;

	FILE*           mFile = NULL;
	gzFile          mZipFile = NULL;

	uint64 			mHalo;

public:
	uint64          totalSeqs;
	uint64			gid = 0;

	//uint64 lastOneReadPos;
	//uint64 lastTwoReadPos;

	//uint64 GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_);
	//uint64 GetPreviousRecordPos(uchar* data_, uint64 pos_, const uint64 size_);
private:
	uint64 FindCutPos(FastaChunk* dataChunk_, uchar* data_, const uint64 size_, const uint64 halo_, SeqInfos& seqInfos);

	void FastaSkipToEol(uchar* data_, uint64& pos_, const uint64 size_)
	{
		//cout << "SkipToEol " << pos_ << " " << size_ << endl;
		ASSERT(pos_ <= size_);

		while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_)
			++pos_;

		if (data_[pos_] == '\r' && pos_ < size_)
		{
			if (data_[pos_ + 1] == '\n')
			{
				usesCrlf = true;
				++pos_;
			}
		}
	}

	void SkipToEol(uchar* data_, uint64& pos_, const uint64 size_)
	{
		//cout << "SkipToEol " << pos_ << " " << size_ << endl;
		ASSERT(pos_ < size_);

		while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_)
			++pos_;

		if (data_[pos_] == '\r' && pos_ < size_)
		{
			if (data_[pos_ + 1] == '\n')
			{
				usesCrlf = true;
				++pos_;
			}
		}
	}
	//跳转到行首
	void SkipToSol(uchar* data_, uint64& pos_, const uint64 size_){
		//std::cout<<"SkipToSol:"<<data_[pos_]<<std::endl;
		ASSERT(pos_ < size_);
		if(data_[pos_] =='\n'){
			--pos_;
		}
		if(data_[pos_] =='\r'){
			usesCrlf = true;
			pos_--;
		}
		//找到换行符
		while(data_[pos_] != '\n' &&  data_[pos_] != '\r'){
			//std::cout<<"pos_;"<<pos_<<std::endl;
			--pos_;
		}
		if(data_[pos_] =='\n'){
			--pos_;
		}
		if(data_[pos_] =='\r'){
			usesCrlf = true;
			pos_--;
		}

	}
	

	bool FindEol(uchar* data_, uint64& pos_, const uint64 size_){
		bool found = false;
		uint64 pos0 = pos_;
		//TODO follow SkipToEol for both \n and \n\r
		while(pos0 < size_){
			if(data_[pos0] == '\n') 
			{
				//std::cout << "enter at: " << pos0 << std::endl;
				//std::string name((char*)(data_ + pos_), pos0-pos_);
				//std::cerr << name << std::endl;
				pos_ = pos0;
				found = true;
				break;
			}else{
				pos0++;
			}
		}
		return found;
	}

};


} // namespace fa

namespace fq
{

class FastqFileReader
{
private:
	static const uint32 SwapBufferSize = 1 << 20; // the longest FASTQ sequence todate is no longer than 1Mbp. 
  FileReader* mFqReader;

public:
	FastqFileReader(const std::string& fileName_, bool isZippedNew = false)
		:	swapBuffer(SwapBufferSize)
		,	bufferSize(0)
		,	eof(false)
		,	usesCrlf(false)
		,	isZipped(isZippedNew)
	{	
		//if(ends_with(fileName_,".gz"))
		if(isZipped){
      isZipped = true;
		}
    mFqReader = new FileReader(fileName_, isZipped);
	}

  FastqFileReader(int fd, bool isZippedNew = false)
      : swapBuffer(SwapBufferSize)
      , bufferSize(0)
      , eof(false)
      , usesCrlf(false)
      , isZipped(isZippedNew)
  {
    mFqReader = new FileReader(fd, isZipped);
  }

  ~FastqFileReader()
	{
		//if( mFile != NULL )
		if(mFile != NULL || mZipFile !=NULL)
			delete mFqReader;
	}

	bool Eof() const
	{
		return eof;
	}

	bool ReadNextChunk(FastqDataChunk* chunk_);
	bool ReadNextPairedChunk(FastqDataChunk* chunk_);
	void Close()
	{
		if(mFile != NULL){
			FCLOSE(mFile);
			mFile = NULL;
		}
		if(mZipFile !=NULL && isZipped){
			gzclose(mZipFile);
			mZipFile=NULL;
		}

	}

	int64 Read_back(byte* memory_, uint64 size_)
	{	
		if(isZipped){
			int64 n = gzread(mZipFile,memory_,size_);
			if(n == -1)
				std::cerr<<"Error to read gzip file" <<std::endl;
			return n;
		}
		else{
		  int64 n = fread(memory_, 1, size_, mFile) ;
		  return n;
		}
  }
  int64 Read(byte *memory_, uint64 size_)
  {
    return mFqReader->Read(memory_, size_);
  }

private:
	core::Buffer	swapBuffer;
	uint64			bufferSize;
	bool			eof;
	bool			usesCrlf;
	bool			isZipped;
	FILE*           mFile = NULL;
	gzFile          mZipFile = NULL;

	uint64 lastOneReadPos;
	uint64 lastTwoReadPos;

	uint64 GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_);
	uint64 GetPreviousRecordPos(uchar* data_, uint64 pos_, const uint64 size_);

	void SkipToEol(uchar* data_, uint64& pos_, const uint64 size_)
	{
		//cout << "SkipToEol " << pos_ << " " << size_ << endl;
		ASSERT(pos_ < size_);

		while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_)
			++pos_;

		if (data_[pos_] == '\r' && pos_ < size_)
		{
			if (data_[pos_ + 1] == '\n')
			{
				usesCrlf = true;
				++pos_;
			}
		}
	}
	//just to the start of a line
	void SkipToSol(uchar* data_, uint64& pos_, const uint64 size_){
		//std::cout<<"SkipToSol:"<<data_[pos_]<<std::endl;
		ASSERT(pos_ < size_);
		if(data_[pos_] =='\n'){
			--pos_;
		}
		if(data_[pos_] =='\r'){
			usesCrlf = true;
			pos_--;
		}
		//find next '\n' or \n\r
		while(data_[pos_] != '\n' &&  data_[pos_] != '\r'){
			//std::cout<<"pos_;"<<pos_<<std::endl;
			--pos_;
		}
		if(data_[pos_] =='\n'){
			--pos_;
		}
		if(data_[pos_] =='\r'){
			usesCrlf = true;
			pos_--;
		}

	}
	


};


} // namespace fq

} // namespace mash

#endif // H_FASTQSTREAM
