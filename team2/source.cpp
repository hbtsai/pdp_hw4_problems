#ifndef CMR_H_
#define CMR_H_
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
using namespace std;

class Camera {
public:
	double K[9];//camera matrix, 3Dcspace to 2Dcspace
	double R[9];//rotation matrix, 3Dgspace to 3D cspace
	double T[3];//camera position,
	void print(){
		cout<<std::fixed;
		cout<<std::setprecision(5);
		for(int i = 0 ; i < 9 ; ++i){
			cout<<setw(11)<<K[i]<<" ";
		}cout<<endl;
		for(int i = 0 ; i < 9 ; ++i){
			cout<<setw(11)<<R[i]<<" ";
		}cout<<endl;
		for(int i = 0 ; i < 3 ; ++i ){
			cout<<setw(11)<<T[i]<<" ";
		}cout<<endl;
	}
	Camera(){}
	Camera(const Camera& cmr){
		memcpy(K,cmr.K,sizeof(double)*9);
		memcpy(R,cmr.R,sizeof(double)*9);
		memcpy(T,cmr.T,sizeof(double)*3);
	}
	Camera(const char* fileName,const char* cmrName,int format){
		if(format==0)readParametersMDB(fileName,cmrName);
	}
	void readParametersMDB(const char* fileName, const char* cmrName){// middle burry format
		ifstream ifs;
		ifs.open(fileName);
		if(!ifs.is_open()){
			cout<<"Open "<<string(fileName)<<" error!!"<<endl;
			return;
		}else{
			int totalcmrNum;
			ifs>>totalcmrNum;
			string buffer;
			string target=cmrName;
			target+=string(".png");
			for(int i = 0 ;i < totalcmrNum ; ++i){
				ifs>>buffer;
				if(target==buffer){
					for(int j = 0 ; j < 9 ;++j)	ifs>>K[j];
					for(int j = 0 ; j < 9 ;++j) ifs>>R[j];
					for(int j = 0 ; j < 3 ;++j) ifs>>T[j];
					break;
				}else {
					for(int j = 0 ; j < 21 ; ++j)ifs>>buffer;
				}
			}
			ifs.close();
			return;
		}
	}
};

#endif


#ifndef VMATH_H_
#define VMATH_H_
//#include "cmr.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
using namespace std;

struct Vec{
	double x;
	double y;
	double z;
};

class Patch{
	Vec* vec;// posision of the points
public:
	size_t _size;
	Patch(int size);
	Vec& operator()(int ii,int jj){
		return vec[ii*_size+jj];
	}
	void initializePos(int ii , int jj);
	~Patch();
};

inline Vec cvt2DCto3DC(Camera& cmr, double photox, double photoy, double targetDepth);
inline Vec cvt3DCto2DC(Camera& cmr,  Vec& posC);
inline Vec cvt3DCto3DG(Camera& cmr, Vec& posC);
inline Vec cvt3DGto3DC(Camera& cmr, Vec& posG);

void setPatch2DCto3DC(Patch& patch, Camera& cmr, Patch& mainPatch, double targetDepth);
void setPatch3DCto3DG(Patch& refPatchG ,Camera&  cmr, Patch& refPatchC);
void setPatch3DGto2DC(Patch& refPatchC ,Camera&  cmr, Patch& refPatchG);

#endif



#ifndef BMP_H_
#define BMP_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct _BMP_header
{
    unsigned short magicno;
    unsigned int filesize;
    unsigned int reserved;
    unsigned int offset;
    unsigned int dib_header_size;
    unsigned int width;
    unsigned int height;
    unsigned short nColorPlanes;
    unsigned short nBits;
    unsigned int BI_RGB;
    unsigned int raw_size;
    unsigned int resolutionH;
    unsigned int resolutionV;
    unsigned int nColors;
    unsigned int imporColor;
}BMP_header;

typedef struct _Pixel
{
    unsigned char R;
    unsigned char G;
    unsigned char B;
}Pixel;

Pixel** readBMP(FILE*,BMP_header*);
void writeBMP(FILE*,FILE*,BMP_header,Pixel*);
void writeBMP2(FILE* fp,BMP_header* bmpheader, Pixel** pixel);

class SuperBMP{
public:
	BMP_header bmpHeader;
	int widthStep;
	Pixel*	   data;
public:
	SuperBMP(const SuperBMP& bmp){
		bmpHeader=bmp.bmpHeader;
		widthStep=get_widthStep();
		data=new Pixel[widthStep*bmpHeader.height];
		memcpy(data,bmp.data,sizeof(Pixel)*widthStep*bmpHeader.height);
	}
	SuperBMP(const BMP_header& header){
		bmpHeader=header;
		widthStep=get_widthStep();
		data=new Pixel[widthStep*bmpHeader.height];
	}
	int get_widthStep(){
		int nPadding = (bmpHeader.width*3/4+1)*4 - bmpHeader.width*3;
		widthStep = nPadding + bmpHeader.width;
		return widthStep;
	}
	bool readImage(const char* fileName){
		FILE* ifile=fopen(fileName,"rb");
		if(ifile==NULL){
			printf("Read %s error!!\n",fileName);fflush(stdout);
			return false;
		}
		Pixel**    rowData;
		rowData = readBMP(ifile, &bmpHeader);
		fclose(ifile);
	    int nPadding = (bmpHeader.width*3/4+1)*4 - bmpHeader.width*3;
		widthStep = nPadding + bmpHeader.width;
		data = new Pixel[widthStep* bmpHeader.height  ];
		for(int i = 0 ; i < bmpHeader.height ; ++i){
			memcpy(data+i*widthStep,rowData[i],bmpHeader.width*sizeof(Pixel));
		}
		for(int i = 0 ; i < bmpHeader.height ; ++i){
			delete rowData[i];
		}
		delete rowData;
		return true;
	}
	void writeImage(const char* fileName){
		FILE* ofile=fopen(fileName,"wb");
		Pixel**    rowData;
		rowData = new Pixel*[bmpHeader.height];
		for(int i = 0 ; i < bmpHeader.height ; ++i){
			rowData[i]=new Pixel[bmpHeader.width];
		}
		for(int i = 0 ; i < bmpHeader.height ; ++i){
			memcpy(rowData[i],data+i*widthStep,bmpHeader.width*sizeof(Pixel));
		}
		writeBMP2(ofile,&bmpHeader,rowData);
		for(int i = 0 ; i < bmpHeader.height ; ++i){
			delete rowData[i];
		}
		delete rowData;
		fclose(ofile);
	}
	void initialize(SuperBMP& bmp){
		bmpHeader=bmp.bmpHeader;
		widthStep =bmp.widthStep;
		data = new Pixel [bmpHeader.height * bmp.widthStep];
		memset(data, 0, bmpHeader.height* bmp.widthStep*sizeof(Pixel));
	}
	inline Pixel& operator()(size_t height, size_t width){
		return data[(bmpHeader.height-height-1)*widthStep+width];
	}
	size_t geth(){return bmpHeader.height;}
	size_t getw(){return bmpHeader.width;}
	SuperBMP(){
		data=NULL;
	}
	SuperBMP(const char* fileName){
		readImage(fileName);
	}
	~SuperBMP(){
		delete data;
	}

};


#endif //BMP_H_
///////////////////////////////BMP.CPP//////////////////////////////////////
#include <cstdlib>
//#include "bmp.h"
#include "string.h"
void writeOn(void* ptr,void* data, int size) {//2bytes
	memcpy(ptr,data,size);
}
unsigned short readshort(unsigned char *input)
{
    unsigned short output = 0;
    int i;
    for(i=0;i<2;i++)
    {
	output += *(input-i+1);
	if(i!=1)
	    output = output << 8;
    }
    return output;
}

unsigned int readint(unsigned char *input)
{
    unsigned int output = 0;
    int i;
    for(i=0;i<4;i++)
    {
	output += *(input-i+3);
	if(i!=3)
	    output = output << 8;
    }
    return output;
}

Pixel** readBMP(FILE* fp,BMP_header* bmpheader)
{
    unsigned char header[54];
    fread(header,sizeof(unsigned char),54,fp);
    bmpheader->magicno = readshort(header);
    bmpheader->filesize = readint(header+2);
    bmpheader->reserved = readint(header+6);
    bmpheader->offset = readint(header+10);
    bmpheader->dib_header_size = readint(header+14);
    bmpheader->width = readint(header+18);
    bmpheader->height = readint(header+22);
    bmpheader->nColorPlanes = readshort(header+26);
    bmpheader->nBits = readshort(header+28);
    bmpheader->BI_RGB = readint(header+30);
    bmpheader->raw_size = readint(header+34);
    bmpheader->resolutionH = readint(header+38);
    bmpheader->resolutionV = readint(header+42);
    bmpheader->nColors = readint(header+46);
    bmpheader->imporColor = readint(header+50);

    const int& height = (int)(bmpheader->height);
    const int& width  = (int)(bmpheader->width);

    // two dimensional image 
    Pixel** pixel = new Pixel*[height];
    for(int h = 0; h < height; ++h)
    {
        pixel[h] = new Pixel[width];
    }

    int nPadding = (bmpheader->width*3/4+1)*4 - bmpheader->width*3;
    if(nPadding == 4)
    {
        nPadding = 0;
    }
    unsigned char *padding = (unsigned char*)malloc(sizeof(unsigned char)*nPadding);
   	memset(padding,0,sizeof(unsigned char)*nPadding); 
    for(int h = 0; h < height; ++h)
    {
        fread(pixel[h] , sizeof(Pixel), bmpheader->width, fp);
        if(nPadding != 0)
            fread(padding,sizeof(unsigned char),nPadding,fp);
    }

    return pixel;
}
void writeBMP2(FILE* fp,BMP_header* bmpheader, Pixel** pixel)
{

    unsigned char header[54];

    writeOn(header,&(bmpheader->magicno),2);//    bmpheader->magicno = readshort(header);
    writeOn(header+2,&(bmpheader->filesize),4);//bmpheader->filesize = readint(header+2);
    writeOn(header+6,&(bmpheader->reserved),4);//bmpheader->reserved = readint(header+6);
   writeOn(header+10,&(bmpheader->offset),4);// bmpheader->offset = readint(header+10);
  writeOn(header+14,&(bmpheader->dib_header_size),4);//    bmpheader->dib_header_size = readint(header+14);
   writeOn(header+18,&(bmpheader->width),4);// bmpheader->width = readint(header+18);
   writeOn(header+22,&(bmpheader->height),4);// bmpheader->height = readint(header+22);
   writeOn(header+26,&(bmpheader->nColorPlanes),2);// bmpheader->nColorPlanes = readshort(header+26);
   writeOn(header+28,&(bmpheader->nBits),2); //bmpheader->nBits = readshort(header+28);
   writeOn(header+30,&(bmpheader->BI_RGB),4); //bmpheader->BI_RGB = readint(header+30);
   writeOn(header+34,&(bmpheader->raw_size),4); //bmpheader->raw_size = readint(header+34);
   writeOn(header+38,&(bmpheader->resolutionH),4); //bmpheader->resolutionH = readint(header+38);
   writeOn(header+42,&(bmpheader->resolutionV),4);// bmpheader->resolutionV = readint(header+42);
   writeOn(header+46,&(bmpheader->nColors),4);// bmpheader->nColors = readint(header+46);
   writeOn(header+50,&(bmpheader->imporColor),4);// bmpheader->imporColor = readint(header+50);
	fwrite(header,sizeof(unsigned char),54,fp);

    const int& height = (int)(bmpheader->height);
    const int& width  = (int)(bmpheader->width);

    // two dimensional image 
    int nPadding = (bmpheader->width*3/4+1)*4 - bmpheader->width*3;
    if(nPadding == 4)
    {
        nPadding = 0;
    }
    unsigned char *padding = (unsigned char*)malloc(sizeof(unsigned char)*nPadding);
   	memset(padding,0,sizeof(unsigned char)*nPadding); 
    
    for(int h = 0; h < height; ++h)
    {
        fwrite(pixel[h] , sizeof(Pixel), bmpheader->width, fp);
        if(nPadding != 0)
            fwrite(padding,sizeof(unsigned char),nPadding,fp);
    }
}

void writeBMP(FILE* frp,FILE* fop,BMP_header bmpheader,Pixel* filtered)
{
    unsigned char originheader[54];
    fseek(frp,0,SEEK_SET);
    fread(originheader,sizeof(unsigned char),54,frp);
    fwrite(originheader,sizeof(unsigned char),54,fop);
    
    int nPadding = (bmpheader.width*3/4+1)*4 - bmpheader.width*3;
    if(nPadding == 4)
	nPadding = 0;
    unsigned char *padding = (unsigned char*)malloc(sizeof(unsigned char)*nPadding);
    int i;
    for(i=0;i<nPadding;i++)
	padding[i] = 0;
    for(i=0;i<(int)bmpheader.height;i++)
    {
	fwrite(filtered+(i*bmpheader.width),sizeof(Pixel),bmpheader.width,fop);
	if(nPadding!=0)
	    fwrite(padding,sizeof(unsigned char),nPadding,fop);
    }
}
//////////////////////////////////////////BMP.CPP end///////////////////

/////////////////////////////////////////VMATH.CPP/////////////////////////////////
//#include "vmath.h"

Patch::Patch(int size){
	_size=size;
	vec = new Vec[_size*_size];
}
void Patch::initializePos(int ii , int jj){
	int patchR = (_size-1)/2;
	for(int i = 0 ; i < _size ; ++i){
		for(int j = 0 ; j < _size; ++j){
			vec[i*_size+j].x = jj-patchR+j;
			vec[i*_size+j].y = ii-patchR+i;
			vec[i*_size+j].z = 1.0;
		}
	}
}
Patch::~Patch(){
	delete vec;
}

inline Vec cvt2DCto3DC(Camera& cmr, double photox, double photoy, double targetDepth){
	Vec tmp;
	double x0 = cmr.K[0*3+2];
	double y0 = cmr.K[1*3+2];
	double fa = cmr.K[0*3+0];
	double fb = cmr.K[1*3+1];

	tmp.z=   targetDepth ;
	tmp.y=   targetDepth* ( photoy - y0 )  /fb;
	tmp.x=   targetDepth* ( photox - x0 )  /fa;
	return tmp;
}
inline Vec cvt3DCto2DC(Camera& cmr,  Vec& posC){
	Vec tmp;
	// 2DC = (K*3DC/Z)
	tmp.x = cmr.K[ 0 * 3 + 0] * posC.x +  cmr.K[ 0 * 3 + 1] * posC.y + cmr.K[ 0 * 3 + 2] * posC.z;
	tmp.y = cmr.K[ 1 * 3 + 0] * posC.x +  cmr.K[ 1 * 3 + 1] * posC.y + cmr.K[ 1 * 3 + 2] * posC.z;
	tmp.z = cmr.K[ 2 * 3 + 0] * posC.x +  cmr.K[ 2 * 3 + 1] * posC.y + cmr.K[ 2 * 3 + 2] * posC.z;
	tmp.x/=tmp.z;
	tmp.y/=tmp.z;
	tmp.z/=tmp.z;
	return tmp;
}

inline Vec cvt3DCto3DG(Camera& cmr, Vec& posC){
	// 3DG = trans(R)* (3DC - T)
	Vec tmp,tmp2;
	tmp2.x=posC.x-cmr.T[0];
	tmp2.y=posC.y-cmr.T[1];
	tmp2.z=posC.z-cmr.T[2];
	tmp.x = cmr.R[ 0 * 3 + 0] * tmp2.x +  cmr.R[ 1 * 3 + 0] * tmp2.y + cmr.R[ 2 * 3 + 0] * tmp2.z ;
	tmp.y = cmr.R[ 0 * 3 + 1] * tmp2.x +  cmr.R[ 1 * 3 + 1] * tmp2.y + cmr.R[ 2 * 3 + 1] * tmp2.z ;
	tmp.z = cmr.R[ 0 * 3 + 2] * tmp2.x +  cmr.R[ 1 * 3 + 2] * tmp2.y + cmr.R[ 2 * 3 + 2] * tmp2.z ;
	return tmp;
}
inline Vec cvt3DGto3DC(Camera& cmr, Vec& posG){
	// 3DC = R* ( 3DG) +T 
	Vec tmp2;
	tmp2.x = cmr.R[ 0 * 3 + 0] * posG.x +  cmr.R[ 0 * 3 + 1] * posG.y + cmr.R[ 0 * 3 + 2] * posG.z + cmr.T[0];
	tmp2.y = cmr.R[ 1 * 3 + 0] * posG.x +  cmr.R[ 1 * 3 + 1] * posG.y + cmr.R[ 1 * 3 + 2] * posG.z + cmr.T[1];
	tmp2.z = cmr.R[ 2 * 3 + 0] * posG.x +  cmr.R[ 2 * 3 + 1] * posG.y + cmr.R[ 2 * 3 + 2] * posG.z + cmr.T[2];
	return tmp2;
}

void setPatch2DCto3DC(Patch& patch, Camera& cmr, Patch& mainPatch, double targetDepth){
	for(int i = 0 ; i < patch._size ; ++i){
		for(int j = 0 ; j < patch._size ; ++j){
			patch(i,j) = cvt2DCto3DC( cmr, mainPatch(i,j).x, mainPatch(i,j).y , targetDepth);//OK
		}
	}
}
void setPatch3DCto3DG(Patch& refPatchG ,Camera&  cmr, Patch& refPatchC){
	for(int i = 0 ; i < refPatchC._size ; ++i){
		for(int j = 0 ; j < refPatchC._size ; ++j){
			refPatchG(i,j) = cvt3DCto3DG( cmr, refPatchC(i,j));
		}
	}
}


void setPatch3DGto2DC(Patch& refPatchC ,Camera&  cmr, Patch& refPatchG){
	Vec tmp;
	for(int i = 0 ; i < refPatchC._size ; ++i){
		for(int j = 0 ; j < refPatchC._size ; ++j){
			tmp = cvt3DGto3DC(cmr, refPatchG(i,j) );
			refPatchC(i,j) = cvt3DCto2DC( cmr, tmp);
		}
	}
}


/////////////////////////////////////////VMATH.CPP end///////////////////////////
#include <iostream>
#include <vector>
#include <cmath>
//#include "bmp.h"
//#include "vmath.h"
#include <iomanip>
using namespace std;

struct SpaceInfo{
	float depthShift;
	float depthStep;
	int MAX_STEP;//256 means:0~255
	int patchR;
};

Pixel getColorITP(SuperBMP& bmp,double i, double j){
	Pixel tmp;
	if(i>0 && i < bmp.bmpHeader.height-1 && j>0 && j<bmp.bmpHeader.width-1 ){
		double lf,rt,up,down,R,G,B;
		lf = j - floor(j);
		rt = 1-lf;
		up = i - floor(i);
		down=1-up;
		Pixel AA = bmp(floor(i),floor(j));
		Pixel BB = bmp(floor(i),ceil(j));
		Pixel CC = bmp(ceil(i),floor(j));
		Pixel DD = bmp(ceil(i), ceil(j));
		R =( (double)AA.R   *rt + (double)BB.R*lf )*down +    ( (double)CC.R*rt + (double)DD.R*lf  )*up;
		G =( (double)AA.G   *rt + (double)BB.G*lf )*down +    ( (double)CC.G*rt + (double)DD.G*lf  )*up;
		B =( (double)AA.B   *rt + (double)BB.B*lf )*down +    ( (double)CC.B*rt + (double)DD.B*lf  )*up;
		tmp.R = (unsigned char)R;
		tmp.G = (unsigned char)G;
		tmp.B = (unsigned char)B;
	}else{
		tmp.R = (unsigned char)255;
		tmp.G = (unsigned char)255;
		tmp.B = (unsigned char)255;
	}
	return tmp;
}
Pixel getColorNEAR(SuperBMP& bmp,double i, double j){
	Pixel tmp;
	if(i>0 && i < bmp.bmpHeader.height-1 && j>0 && j<bmp.bmpHeader.width-1 ){
		int lf = floor(i);
		int up = floor(j);
		if( i-lf>=0.5)i=ceil(i);
		else i = floor(i);
		if( j-up>=0.5)j=ceil(j);
		else j = floor(j);
		tmp = bmp(i,j);
	}else{
		tmp.R = (unsigned char)255;
		tmp.G = (unsigned char)255;
		tmp.B = (unsigned char)255;
	}
	return tmp;
}

void setPatchColor(Patch& color, SuperBMP& bmp, Patch& refPatch){
	int patchR=(color._size-1)/2;
	Pixel tmp;
	for(int i = 0 ; i<color._size ; ++i){
		for(int j = 0 ; j < color._size ; ++j){
//			tmp = getColorITP(bmp,refPatch(i,j).y,refPatch(i,j).x);
			tmp = getColorNEAR(bmp,refPatch(i,j).y,refPatch(i,j).x);
			color(i,j).x =(unsigned char) tmp.R;
			color(i,j).y =(unsigned char) tmp.G;
			color(i,j).z =(unsigned char) tmp.B;
		}
	}
}
double diffSqarColor(Patch& mainPatchColor,Patch& refPatchCimgColor   ){
	double tmp=0,x,y,z;
	for(int i = 0 ; i < mainPatchColor._size;++i){
		for(int j = 0 ; j < mainPatchColor._size;++j){
			x=mainPatchColor(i,j).x-refPatchCimgColor(i,j).x;
			y=mainPatchColor(i,j).y-refPatchCimgColor(i,j).y;
			z=mainPatchColor(i,j).z-refPatchCimgColor(i,j).z;
			tmp+=x*x+y*y+z*z;
		}
	}
	return tmp;
}

int main(int argc, char **argv){
	// 1. load configure file
#ifdef USE_GOOD_IO
	SpaceInfo spaceInfo;
	FILE *fp;
	char str[100];

	char file_name[30];
	int neighbor;
	char bmp_name[30][30];

	char result_name[30];

	if(argc != 2) {
		printf("Usage: ./main [input file]\n");
		return 1;
	}
	fp = fopen(argv[1], "r");
	if(fp == NULL) {
		printf("File \"%s\" does not exist!\n", argv[1]);
		return 1;
	}

	int line = 0;
	int num_bmp = 0;

	while(fgets(str, 100, fp) != NULL) {
		if(line == 0) {
			sscanf(str, "%s", file_name);
		} else if(line == 1) {
			sscanf(str, "%*s %f", &spaceInfo.depthShift);
		} else if(line == 2) {
			sscanf(str, "%*s %f", &spaceInfo.depthStep);
		} else if(line == 3) {
			sscanf(str, "%*s %d", &spaceInfo.MAX_STEP);
		} else if(line == 4) {
			sscanf(str, "%*s %d", &spaceInfo.patchR);
		} else if(line == 5) {
			sscanf(str, "%s", result_name);
		} else if(line == 6) {
			sscanf(str, "%*s %d", &neighbor);
		} else {
			sscanf(str, "%s", bmp_name[num_bmp++]);
		}

		line++;
	}

	fclose(fp);

	const int patchSize =spaceInfo.patchR*2+1;
	// 2. load photos and camera parameter 
	vector<SuperBMP> photo,mask;
	vector<Camera>	 cmr;
	SuperBMP result;
	char prefix_str[] = "./input/";
	char str_photo[30];
	char str_mask[30];
	char str_cmr[30];

	for(int i = 0; i < num_bmp; i++) {
		sprintf(str_photo, "%s%s%s", prefix_str, bmp_name[i], ".bmp");
		sprintf(str_mask, "%s%s%s", prefix_str, bmp_name[i], ".msk.bmp");
		sprintf(str_cmr, "%s%s", prefix_str, file_name);

		photo.push_back(SuperBMP(str_photo));
		mask.push_back(SuperBMP(str_mask));
		cmr.push_back(Camera(str_cmr, bmp_name[i], 0));
	}
//	for(int i = 0 ; i < cmr.size() ; ++i)cmr[i].print();
#else
	SpaceInfo spaceInfo;
	vector<SuperBMP> photo,mask;
	vector<Camera>	 cmr;
	char result_name[30]="resultmmm.bmp";
	SuperBMP result;
	{	
		int numPIC;
		cin>>spaceInfo.depthShift;
		cin>>spaceInfo.depthStep;
		cin>>spaceInfo.MAX_STEP;
		cin>>spaceInfo.patchR;
		cin>>numPIC;
		cmr.resize(numPIC);
		for(int i = 0 ; i < numPIC; ++i){
			for(int j = 0 ; j < 9 ;++j)
				cin>>cmr[i].K[j];
			for(int j = 0 ; j < 9 ;++j)
				cin>>cmr[i].R[j];
			for(int j = 0 ; j < 3 ;++j)
				cin>>cmr[i].T[j];
//			cmr[i].print();
		}
		//photo
		for(int id = 0 ; id < numPIC ; ++id){
			BMP_header header;
			cin>>header.magicno>>header.filesize>>header.reserved>>header.offset>>header.dib_header_size>>header.width
				>>header.height	>>header.nColorPlanes	>>header.nBits	>>header.BI_RGB	>>header.raw_size	>>header.resolutionH
				>>header.resolutionV	>>header.nColors	>>header.imporColor;
			photo.push_back(SuperBMP(header));
			SuperBMP& tmpBmp=photo[id];
			for(int i = 0 ; i < tmpBmp.geth();++i){
				for(int j = 0 ; j < tmpBmp.getw();++j){
					int value;
					cin>>value;tmpBmp(i,j).R=(unsigned char)value;
					cin>>value;tmpBmp(i,j).G=(unsigned char)value;
					cin>>value;tmpBmp(i,j).B=(unsigned char)value;
				}
			}
		}
		//mask
		for(int id = 0 ; id < numPIC ; ++id){
			BMP_header header;
			cin>>header.magicno>>header.filesize>>header.reserved>>header.offset>>header.dib_header_size>>header.width
				>>header.height	>>header.nColorPlanes	>>header.nBits	>>header.BI_RGB	>>header.raw_size	>>header.resolutionH
				>>header.resolutionV	>>header.nColors	>>header.imporColor;
			mask.push_back(SuperBMP(header));
			SuperBMP& tmpBmp=mask[id];
			for(int i = 0 ; i < tmpBmp.geth();++i){
				for(int j = 0 ; j < tmpBmp.getw();++j){
					int value;
					cin>>value;tmpBmp(i,j).R=(unsigned char)value;
					cin>>value;tmpBmp(i,j).G=(unsigned char)value;
					cin>>value;tmpBmp(i,j).B=(unsigned char)value;
				}
			}
		}
	}
	const int patchSize =spaceInfo.patchR*2+1;
#endif

////////////////////////////////////////////////////////DUMPING...
// input: spaceInfo, { num,w,h }{ cmr paras }   photo, mask, cmr 
#ifdef DUMPIN
	ofstream sampleI("sampleIn.txt");
	sampleI<<std::setprecision(20);
	sampleI<<spaceInfo.depthShift<<" "<<spaceInfo.depthStep<<" "<<spaceInfo.MAX_STEP<<" "<<spaceInfo.patchR<<" "<<endl;
	sampleI<<photo.size()<<photo[0].getw()<<" "<<photo[0].geth()<<endl;
	for(int i = 0 ; i < cmr.size() ; ++i){
		for(int j = 0 ; j < 9 ; ++j)
			sampleI<<cmr[i].K[j]<<" ";
		for(int j = 0 ; j < 9 ; ++j)
			sampleI<<cmr[i].R[j]<<" ";
		for(int j = 0 ; j < 3 ; ++j)
			sampleI<<cmr[i].T[j]<<" ";
		sampleI<<endl;
	}
	for(int id = 0 ; id < photo.size() ; ++id){
		for(int i = 0 ; i < photo[id].getw();++i){
			for(int j = 0 ; j < photo[id].geth();++j){
				sampleI<<(int)photo[id](i,j).R<<" ";
				sampleI<<(int)photo[id](i,j).G<<" ";
				sampleI<<(int)photo[id](i,j).B<<" ";
			}
			sampleI<<endl;
		}sampleI<<endl;
	}
	for(int id = 0 ; id < mask.size() ; ++id){
		for(int i = 0 ; i < mask[id].getw();++i){
			for(int j = 0 ; j < mask[id].geth();++j){
				sampleI<<(int)mask[id](i,j).R<<" ";
				sampleI<<(int)mask[id](i,j).G<<" ";
				sampleI<<(int)mask[id](i,j).B<<" ";
			}sampleI<<endl;
		}sampleI<<endl;
	}
	sampleI.close();
#endif

	// 3. create result
	// start body!!!!//////////////////////////////////////////////////////////////////////////
	result.initialize(photo[0]);


	SuperBMP& mainIMG=photo[0];
	for(int i = spaceInfo.patchR  ; i < mainIMG.geth()-spaceInfo.patchR ; ++i){
		Patch mainPatch(patchSize),mainPatchColor(patchSize);
		Patch refPatchC(patchSize);
		Patch refPatchG(patchSize);
		Patch refPatchCimg(patchSize),refPatchCimgColor(patchSize);
		vector<double>  vote ;
		vote.resize(spaceInfo.MAX_STEP,0);
		for(int j = spaceInfo.patchR ; j < mainIMG.getw()-spaceInfo.patchR; ++j){
			if(mask[0](i,j).R <100) {
				Pixel tmpResult;
				tmpResult.R=tmpResult.G=tmpResult.B=(unsigned char)255;
				result(i,j)=tmpResult;
				continue;
			}
			for(int dStep = 0; dStep < spaceInfo.MAX_STEP; ++dStep){
				double targetDepth= spaceInfo.depthStep * dStep + spaceInfo.depthShift;
				//1. set patch position
				mainPatch.initializePos(i,j);
				setPatchColor(mainPatchColor,mainIMG,mainPatch);
				setPatch2DCto3DC(refPatchC,cmr[0], mainPatch,targetDepth);
				setPatch3DCto3DG(refPatchG,cmr[0], refPatchC);
				for(int cmrID=1;cmrID<cmr.size() ; ++cmrID){
					setPatch3DGto2DC(refPatchCimg,cmr[cmrID],refPatchG);
					setPatchColor(refPatchCimgColor ,photo[cmrID],refPatchCimg);
					vote[dStep]+=diffSqarColor( mainPatchColor,refPatchCimgColor   );
				}
			}// for all step
			double minVal = vote[0];
			int minPos    = 0;
			for(int dStep = 0 ; dStep < spaceInfo.MAX_STEP;++dStep){
				if(minVal>vote[dStep]){
					minVal=vote[dStep];
					minPos=dStep;
				}
				vote[dStep]=0;
			}
			result(i,j).R=(unsigned char)minPos;
			result(i,j).G=(unsigned char)minPos;
			result(i,j).B=(unsigned char)minPos;
		}
	}
//#ifdef DUMP_BMP_RESULT

	result.writeImage(result_name);
//#endif
		for(int i = 0 ; i < result.geth();++i){
			for(int j = 0 ; j < result.getw();++j){
				cout<<(int)result(i,j).R<<" ";
				cout<<(int)result(i,j).G<<" ";
				cout<<(int)result(i,j).B<<" ";
			}cout<<endl;
		}cout<<endl;

	return 0;
}
