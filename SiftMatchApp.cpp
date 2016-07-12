//1)too strong block 2)should scale and orientation be used in block match
//
#include "PrjLog.hpp"

#include "../AtMatch/GlobalUtil_AtMatch.h"
#include "../AtMatch/SiftMatch.h"
#include "../AtMatch/StdImage.h"
#include "../AtMatch/MatchFile.h"
#include "../AtMatch/ImagesCoorAssociate.hpp"
#include "algorithms/Pretreatment.h"
#include "algorithms/Transformation.h"

#include "SiftMatchDef.h"

#define PI							3.1415926
#define RELEASE_NAME				"1.2.2"	
#define RELEASE_DATE				"2016/5/11"	

//#define SIFT_PYRM_IMAGE

const char* SiftMatchVersionInfo(){
	return "SiftMatch "
#ifdef SIFT_PYRM_IMAGE
		" (ImagePyramid) "
#endif
		RELEASE_NAME
#ifdef _DEBUG
		" debug "
#else
		" released "
#endif
		RELEASE_DATE;
}
/*  ******************************************************************* */
/*                               Usage()                                */
/* ******************************************************************** */

static void Usage(const char* pszErrorMsg = NULL, int bShort = TRUE)

{
	printf("Usage: SiftMatchEx [--help-general] [--long-usage]\n"
		"       [-enable-gpu(e/m)]\n"
		"       [-ws workspace] [-output match_result_file] [-log log_file]\n"
		"       [-range(G/B) search_range]\n"
		"       [-sift{D/R} distance/ratio]\n"
		"       [-ue] [-refine(h/f)] [-no-mutual]\n"
		"       [-scale ratio_min ratio_max] [-orientation diff_min(r) diff_max(r)]\n"
		"       src_dataset dst_dataset\n"
		"       [-id{s/d} image_id]\n"
		"       [-geo(s/d) (geo_file)]\n"
		"       [-rc(s/d) stCol stRow edCol edRow]\n"// this is only available for extract feature currently (should take exchange image into account).
		"       [-dem dem_file]\n");

	if (!bShort)
	{
		printf("\n%s\n\n", SiftMatchVersionInfo());
		printf("The following format drivers are configured :\n");
		const char* strDriver[] = { "GTiff(*.tif,*tiff)", "HFA(*.img)", "PNG(*.png)", "JPEG(*.jpg)", "PCIDSK(*.pix)","BMP(*.bmp)" };
		for (int iDr = 0; iDr < sizeof(strDriver) / sizeof(const char*); iDr++)
		{
			printf("%s;\t", strDriver[iDr]);
		}
		printf("\nThe following format GeoFile are configured :\n");
		printf("RPC;\n");
		printf("Default value: siftDist(%.2f) siftRatio(%.2f) GeoRange(%d) BlockRange(%d)\n", DEFAULT_SIFT_DIST, DEFAULT_SIFT_RATIO, RANGE_GEO, RANGE_BLOCK);
	}

	if (pszErrorMsg != NULL)
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);

	exit(1);
}

#ifndef _MPT
#define _MPT
typedef struct tagMPt{
	char nameL[POINTID_BUFFER_SIZE];
	char nameR[POINTID_BUFFER_SIZE];
	float	xl;
	float	yl;
	float	xr;
	float	yr;
}MPt;
#endif

MPt*		SiftMatch_Block(
	int* ptSum,
	CStdImage* imgL, const char* lpcstrSiftL, int stBlkColL, int stBlkRowL, int nBlkColsL, int nBlkRowsL,
	CStdImage* imgR, const char* lpcstrSiftR,
	CImagesCoorT*		pCvt,
	float siftDist, float siftRatio, int mutual_best_match
	){
	LogPrint(0, "\n<----SiftMatch_Block starting[range= %d]---->", pCvt->get_range());
	*ptSum = 0;
	CSiftMatch matcher;
	if (!matcher.InitEnvi(GlobalParam::g_sift_match_gpu)) return NULL;

	int nColNum, nRowNum;
	int* rc_split = Pretreatment::split_image(nColNum, nRowNum, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, 8,
		GlobalParam::g_match_buffer_size, Pretreatment::RECT_FIRST);
	int blk_num = nColNum*nRowNum;

	MPt* match_result = NULL;
	int *pSplitRc = rc_split;
	for (int i = 0; i < blk_num; i++, pSplitRc += 4){
		LogPrint(0, ">>>>Block[%d/%d]", i + 1, blk_num);
		int stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR;
		if (!pCvt->overlap_on_right(pSplitRc[0], pSplitRc[1], pSplitRc[0] + pSplitRc[2], pSplitRc[1] + pSplitRc[3],
			&stBlkColR, &stBlkRowR, &nBlkColsR, &nBlkRowsR)){
			LogPrint(0, "No coverage.Continue.");
			continue;
		}
		nBlkColsR = nBlkColsR - stBlkColR + 1;	nBlkRowsR = nBlkRowsR - stBlkRowR + 1;
		LogPrint(0, "(%s)[%d,%d]+[%d,%d]<=>(%s)[%d,%d]+[%d,%d]", imgL->GetImageName(), pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3], imgR->GetImageName(), stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);
		int ptNumL = imgL->ReadSiftFile(lpcstrSiftL, pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3]);
		if (ptNumL < 1)	{ LogPrint(0, "Continue."); continue; }
		int ptNumR = imgR->ReadSiftFile(lpcstrSiftR, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);
		if (ptNumR < 1)	{ LogPrint(0, "Continue."); continue; }

		int(*match_buf)[2] = new int[ptNumL>ptNumR ? ptNumL : ptNumR][2];
		matcher.SetDescriptors(0, ptNumL, imgL->GetSiftDescriptors());
		matcher.SetDescriptors(1, ptNumR, imgR->GetSiftDescriptors());
		
		int match_num = matcher.GetSiftMatch(match_buf, siftDist, siftRatio, mutual_best_match);
		if (match_num>0){
			MPt* mpt = new MPt[*ptSum + match_num];
			if (match_result) { memcpy(mpt, match_result, *ptSum*sizeof(MPt)); delete match_result; }	
			match_result = mpt;
			MPt* pMPT = mpt + *ptSum;	int sz = 0;
			int i;	float lx0 = -99, ly0 = -99, rx0 = -99, ry0 = -99;
			for (i = 0; i < match_num; i++){
				float* xyl = imgL->GetPtLocation() + match_buf[i][0] * 2;
				float* xyr = imgR->GetPtLocation() + match_buf[i][1] * 2;
				pMPT->xl = *xyl;	pMPT->yl = *(xyl + 1);
				pMPT->xr = *xyr;	pMPT->yr = *(xyr + 1);
				if (fabs(pMPT->xl - lx0) < 1e-5 && fabs(pMPT->yl - ly0) < 1e-5
					&& fabs(pMPT->xr - rx0) < 1e-5 && fabs(pMPT->yr - ry0) < 1e-5) continue;
				lx0 = pMPT->xl; ly0 = pMPT->yl;	rx0 = pMPT->xr;	ry0 = pMPT->yr;
				strcpy(pMPT->nameL, imgL->get_pt_name(match_buf[i][0]));
				strcpy(pMPT->nameR, imgR->get_pt_name(match_buf[i][1]));
				pMPT++;	sz++;
			}			
			LogPrint(0, "Remove repeat correspond points.Final Num= %d", sz);
			*ptSum += sz;
		}
		delete[] match_buf;
	}

	delete[] rc_split;
	LogPrint(0, "<----SiftMatch_Block END[total num= %d][state= %s]---->\n", *ptSum, *ptSum<1?"Fail":"OK");
	return match_result;
}
MPt*		SiftMatch_Geo(
	int* ptSum,
	CStdImage* imgL, const char* lpcstrSiftL, int stBlkColL, int stBlkRowL, int nBlkColsL, int nBlkRowsL,
	CStdImage* imgR, const char* lpcstrSiftR,
	CImagesCoorT*		pCvt,
	float H[3][3], float hdistmax,
	float F[3][3], float fdistmax,
	float ratiomin_scale, float ratiomax_scale,
	float ratiomin_ori, float ratiomax_ori,
	float siftDist, float siftRatio,int mutual_best_match
	){
	LogPrint(0, "\n<----SiftMatch_Geo starting[range= %d]---->", pCvt->get_range());

	*ptSum = 0;
	CSiftMatch matcher;
	if (!matcher.InitEnvi(GlobalParam::g_sift_match_gpu)) return NULL;

	bool bScale = CheckLimitedScale(ratiomin_scale, ratiomax_scale);
	bool bOri = CheckLimitedOrientation(ratiomin_ori, ratiomax_ori);
	if (H) LogPrint(0, "Residual error of H = %.2f", hdistmax);
	if (F) LogPrint(0, "Residual error of F = %.2f", fdistmax);
	if (bScale)	LogPrint(0, "Ratio range of sift scale feature = [%.2f,%.2f]", ratiomin_scale, ratiomax_scale);
	if (bOri)	LogPrint(0, "Ratio range of sift orientation feature = [%.2f,%.2f]", ratiomin_ori, ratiomax_ori);

	int nColNum, nRowNum;
	int* rc_split = Pretreatment::split_image(nColNum, nRowNum, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, 8,
		GlobalParam::g_match_buffer_size, Pretreatment::RECT_FIRST);
	int blk_num = nColNum*nRowNum;

	MPt* match_result = NULL;
	int *pSplitRc = rc_split;
	for (int i = 0; i < blk_num; i++, pSplitRc += 4){
		LogPrint(0, ">>>>Block[%d/%d]", i + 1, blk_num);
		int stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR;
		if (!pCvt->overlap_on_right(pSplitRc[0], pSplitRc[1], pSplitRc[0] + pSplitRc[2], pSplitRc[1] + pSplitRc[3],
			&stBlkColR, &stBlkRowR, &nBlkColsR, &nBlkRowsR)){
			LogPrint(0, "No coverage.Continue.");
			continue;
		}
		nBlkColsR = nBlkColsR - stBlkColR + 1;	nBlkRowsR = nBlkRowsR - stBlkRowR + 1;
		LogPrint(0, "(%s)[%d,%d]+[%d,%d]<=>(%s)[%d,%d]+[%d,%d]", imgL->GetImageName(), pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3], imgR->GetImageName(), stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);
		int ptNumL = imgL->ReadSiftFile(lpcstrSiftL, pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3]);
		if (ptNumL < 1)	{ LogPrint(0, "Continue."); continue; }
		int ptNumR = imgR->ReadSiftFile(lpcstrSiftR, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);
		if (ptNumR < 1)	{ LogPrint(0, "Continue."); continue; }

		int(*match_buf)[2] = new int[ptNumL>ptNumR ? ptNumL : ptNumR][2];
		matcher.SetDescriptors(0, ptNumL, imgL->GetSiftDescriptors());
		matcher.SetDescriptors(1, ptNumR, imgR->GetSiftDescriptors());

		matcher.SetLocation(0, imgL->GetPtLocation());
		matcher.SetLocation(1, imgR->GetPtLocation());		

		if (bScale||bOri) {
			matcher.SetFeature(0, imgL->GetSiftFeature());
			matcher.SetFeature(1, imgR->GetSiftFeature());
		}

		int match_num = 0;
		if (!H&&!F) match_num = matcher.GetLimitedXformSiftMatch(match_buf, geo_xform_err, pCvt, pCvt->get_range(),
					ratiomin_scale, ratiomax_scale, ratiomin_ori, ratiomax_ori, siftDist, siftRatio, mutual_best_match);
		else match_num = matcher.GetLimitedGuidedSiftMatch(match_buf, H, F, ratiomin_scale, ratiomax_scale, ratiomin_ori, ratiomax_ori, siftDist, siftRatio, hdistmax, fdistmax, mutual_best_match);

		if (match_num>0){
			MPt* mpt = new MPt[*ptSum + match_num];
			if (match_result) { memcpy(mpt, match_result, *ptSum*sizeof(MPt)); delete match_result; }
			match_result = mpt;
			MPt* pMPT = mpt + *ptSum;	int sz = 0;
			int i;	float lx0 = -99, ly0 = -99, rx0 = -99, ry0 = -99;
			for (i = 0; i < match_num; i++){
				float* xyl = imgL->GetPtLocation() + match_buf[i][0] * 2;
				float* xyr = imgR->GetPtLocation() + match_buf[i][1] * 2;
				pMPT->xl = *xyl;	pMPT->yl = *(xyl + 1);
				pMPT->xr = *xyr;	pMPT->yr = *(xyr + 1);
				if (fabs(pMPT->xl - lx0) < 1e-5 && fabs(pMPT->yl - ly0) < 1e-5
					&& fabs(pMPT->xr - rx0) < 1e-5 && fabs(pMPT->yr - ry0) < 1e-5) continue;
				lx0 = pMPT->xl; ly0 = pMPT->yl;	rx0 = pMPT->xr;	ry0 = pMPT->yr;
				strcpy(pMPT->nameL, imgL->get_pt_name(match_buf[i][0]));
				strcpy(pMPT->nameR, imgR->get_pt_name(match_buf[i][1]));
				pMPT++;	sz++;
			}
			LogPrint(0, "Remove repeat correspond points.Final Num= %d", sz);
			*ptSum += sz;
		}
		delete[] match_buf;
	}

	delete[] rc_split;
	LogPrint(0, "<----SiftMatch_Geo END[total num= %d][state= %s]---->\n", *ptSum, *ptSum<1 ? "Fail" : "OK");
	return match_result;
}

#include <vector>
using namespace std;
inline float CheckSiftScaleMinimum(const float* siftScale,int nBufferSpace, int nFeaNum){
	int num = 1000>nFeaNum ? nFeaNum : 1000;
	float mini = 1000000.0f;
	for (int i = 0; i < num; i++){
		if (*siftScale < mini) mini = *siftScale;
		siftScale += nBufferSpace;
	}
	return mini < 1000000.0f ? mini : 1;
}

inline int	 PickSift_ScaleMini(
	const PtName* name_org, const PtLoc* loc_org, const PtFea* fea_org, const SiftDes* des_org,int ptNum,float scaleMini,
	vector<PtName>&	name, vector<PtLoc>& loc, vector<SiftDes>& des 
	){
	for (int i = 0; i < ptNum; i++ ){
		if (fea_org->scale >= scaleMini){
			name.push_back(*name_org);	loc.push_back(*loc_org);	des.push_back(*des_org);
		}
		name_org++;	loc_org++;	fea_org++;	des_org++;
	}
	return name.size();
}

static int	PickPyrmSift(
	const PtName* name_org, const PtLoc* loc_org, const PtFea* fea_org, const SiftDes* des_org, int ptNum,float zoomRate,
	vector<PtName>&	name, vector<PtLoc>& loc, vector<SiftDes>& des
	){
	float scaleBase = CheckSiftScaleMinimum((const float*)fea_org, 2, ptNum);
	scaleBase = scaleBase / zoomRate;
	name.clear();	name.reserve(ptNum / 5); loc.clear();	loc.reserve(ptNum / 5); des.clear();	des.reserve(ptNum / 5);
	return PickSift_ScaleMini(name_org, loc_org, fea_org, des_org, ptNum, scaleBase, name, loc, des);
}


MPt*		SiftMatch_Pyrm(
	int* ptSum,
	CStdImage* imgL, const char* lpcstrSiftL, int stBlkColL, int stBlkRowL, int nBlkColsL, int nBlkRowsL,
	CStdImage* imgR, const char* lpcstrSiftR, int stBlkColR, int stBlkRowR, int nBlkColsR, int nBlkRowsR,
	CImagesCoorT*		pCvt, float zoomRateL, float zoomRateR,
	float siftDist, float siftRatio, int mutual_best_match,
	const char* lpcstrTempDir
	){
	LogPrint(0, "\n<----SiftMatch_Pyrm starting[range= %d;zoom= %f&%f]---->", pCvt->get_range(), zoomRateL,zoomRateR);
	*ptSum = 0;	MPt* match_result = NULL;
	CSiftMatch matcher;
	if (!matcher.InitEnvi(GlobalParam::g_sift_match_gpu)) return NULL;
	
#ifdef SIFT_PYRM_IMAGE
	CStdImage imgPyrmL, imgPyrmR;
	char	strSiftPyrmL[512], strSiftPyrmR[512];
	int nPyrmColsL, nPyrmRowsL, nPyrmColsR, nPyrmRowsR;

	char strFile[512];	strcpy(strFile, lpcstrTempDir);	
	char* pS = strFile + strlen(strFile);
	{
		sprintf(pS, "/%s_%d_%d_%d_%d_%.f.tif", imgL->GetImageName(),
			stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, 1 / zoomRateL);
		BYTE* pBuf = NULL;	FILE* fp = NULL;
		pBuf = new BYTE[int(nBlkColsL*zoomRateL*nBlkRowsL) + 8];	if (!pBuf) return NULL;
		imgL->ReadGray8(pBuf, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, zoomRateL, &nPyrmColsL, &nPyrmRowsL);
		if (!SaveImageFile(strFile, pBuf, nPyrmColsL, nPyrmRowsL, 1) || !imgPyrmL.Open(strFile))
		{
			delete pBuf;
			return NULL;
		}
		delete pBuf;
		strcpy(strrchr(strFile, '.'), ".info");
		fp = fopen(strFile, "w"); if (!fp) return NULL;	fprintf(fp, "%d\t%d\t%d\t%d\t%f", stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, zoomRateL); fclose(fp);
		strcpy(strrchr(strFile, '.'), ".sift");	strcpy(strSiftPyrmL, strFile);
		if (imgPyrmL.ExtractSift2PtFile(strSiftPyrmL, -1) < 1) return NULL;

		sprintf(pS, "/%s_%d_%d_%d_%d_%.f.tif", imgR->GetImageName(),
			stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, 1 / zoomRateR);
		pBuf = new BYTE[int(nBlkColsR*zoomRateR*nBlkRowsR) + 8];	if (!pBuf) return NULL;
		imgR->ReadGray8(pBuf, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, zoomRateR, &nPyrmColsR, &nPyrmRowsR);
		if (!SaveImageFile(strFile, pBuf, nPyrmColsR, nPyrmRowsR, 1) || !imgPyrmR.Open(strFile))
		{
			delete pBuf;
			return NULL;
		}
		delete pBuf;
		strcpy(strrchr(strFile, '.'), ".info");
		fp = fopen(strFile, "w"); if (!fp) return NULL;	fprintf(fp, "%d\t%d\t%d\t%d\t%f", stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, zoomRateR); fclose(fp);
		strcpy(strrchr(strFile, '.'), ".sift");	strcpy(strSiftPyrmR, strFile);
		if (imgPyrmR.ExtractSift2PtFile(strSiftPyrmR, -1) < 1) return NULL;
	}

	int nColNum, nRowNum;
	int* rc_split = Pretreatment::split_image(nColNum, nRowNum, 0, 0, nPyrmColsL, nPyrmRowsL, 8,
		GlobalParam::g_match_buffer_size, Pretreatment::RECT_FIRST);
	int blk_num = nColNum*nRowNum;

	int *pSplitRc = rc_split;
	for (int i = 0; i < blk_num; i++, pSplitRc += 4){
		LogPrint(0, ">>>>Pyrmblock[%d/%d]", i + 1, blk_num);
		int scr, srr, ncr, nrr;
		if (!pCvt->overlap_on_right(stBlkColL + pSplitRc[0] / zoomRateL, stBlkRowL + pSplitRc[1] / zoomRateL,
			stBlkColL + (pSplitRc[0] + pSplitRc[2]) / zoomRateL, stBlkRowL + (pSplitRc[1] + pSplitRc[3]) / zoomRateL,
			&scr, &srr, &ncr, &nrr)){
			LogPrint(0, "No coverage.Continue.");
			continue;
		}
		ncr = int((ncr - scr + 1)*zoomRateR);	nrr = int((nrr - srr + 1)*zoomRateR);
		scr -= stBlkColR;	srr -= stBlkRowR;
				
// 		LogPrint(0, "left[%d,%d]+[%d,%d]<=>right[%d,%d]+[%d,%d]", 
// 			pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3], stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);

		int ptNumL = imgPyrmL.ReadSiftFile(strSiftPyrmL, pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3]);
		if (ptNumL < 1)	{ LogPrint(0, "Continue."); continue; }
		int ptNumR = imgPyrmR.ReadSiftFile(strSiftPyrmR, scr, srr, ncr, nrr);
		if (ptNumR < 1)	{ LogPrint(0, "Continue."); continue; }

		int(*match_buf)[2] = new int[ptNumL>ptNumR ? ptNumL : ptNumR][2];
		matcher.SetDescriptors(0, ptNumL, imgPyrmL.GetSiftDescriptors());
		matcher.SetDescriptors(1, ptNumR, imgPyrmR.GetSiftDescriptors());

		int match_num = matcher.GetSiftMatch(match_buf, siftDist, siftRatio, mutual_best_match);

		if (match_num > 0){
			MPt* mpt = new MPt[*ptSum + match_num];
			if (match_result) { memcpy(mpt, match_result, *ptSum*sizeof(MPt)); delete match_result; }
			match_result = mpt;
			MPt* pMPT = mpt + *ptSum;	int sz = 0;
			int i;	float lx0 = -99, ly0 = -99, rx0 = -99, ry0 = -99;
			for (i = 0; i < match_num; i++){
				float* xyl = imgPyrmL.GetPtLocation() + match_buf[i][0] * 2;
				float* xyr = imgPyrmR.GetPtLocation() + match_buf[i][1] * 2;
				pMPT->xl = *xyl;	pMPT->yl = *(xyl + 1);
				pMPT->xr = *xyr;	pMPT->yr = *(xyr + 1);
				if (fabs(pMPT->xl - lx0) < 1e-5 && fabs(pMPT->yl - ly0) < 1e-5
					&& fabs(pMPT->xr - rx0) < 1e-5 && fabs(pMPT->yr - ry0) < 1e-5) continue;
				lx0 = pMPT->xl; ly0 = pMPT->yl;	rx0 = pMPT->xr;	ry0 = pMPT->yr;
				pMPT->xl = pMPT->xl / zoomRateL + stBlkColL;	pMPT->yl = pMPT->yl / zoomRateL + stBlkRowL;
				pMPT->xr = pMPT->xr / zoomRateR + stBlkColR;	pMPT->yr = pMPT->yr / zoomRateR + stBlkRowR;
				strcpy(pMPT->nameL, imgPyrmL.get_pt_name(match_buf[i][0]));
				strcpy(pMPT->nameR, imgPyrmR.get_pt_name(match_buf[i][1]));
				pMPT++;	sz++;
			}
			LogPrint(0, "Remove repeat correspond points.Final Num= %d", sz);
			*ptSum += sz;
		}
		delete[] match_buf;
	}
	delete[] rc_split;
#else	
	if (zoomRateL > 1 / 1.5) zoomRateL = 1 / 1.5;	if (zoomRateL < 1 / 10.0) zoomRateL = 1 / 10.0f;
	if (zoomRateR > 1 / 1.5) zoomRateR = 1 / 1.5;	if (zoomRateR < 1 / 10.0) zoomRateR = 1 / 10.0f;
	LogPrint(0, "[%s]zoomRate= %f;[%s]zoomRate= %f;", imgL->GetImageName(), zoomRateL, imgR->GetImageName(), zoomRateR);
	int nColNum, nRowNum;
	int* rc_split = Pretreatment::split_image(nColNum, nRowNum, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, 8,
		int(GlobalParam::g_match_buffer_size / zoomRateL), Pretreatment::RECT_FIRST);
	int blk_num = nColNum*nRowNum;

	vector<PtName>	nameL, nameR;
	vector<PtLoc>	locL, locR;
	vector<SiftDes>	desL, desR;
	int *pSplitRc = rc_split;
	for (int i = 0; i < blk_num; i++, pSplitRc += 4){
		LogPrint(0, ">>>>Block[%d/%d]", i + 1, blk_num);
		int scr, srr, ncr, nrr;
		if (!pCvt->overlap_on_right(pSplitRc[0], pSplitRc[1], pSplitRc[0] + pSplitRc[2], pSplitRc[1] + pSplitRc[3],
			&scr, &srr, &ncr, &nrr)){
			LogPrint(0, "No coverage.Continue.");
			continue;
		}		
		ncr = ncr - scr + 1;	nrr = nrr - srr + 1;
		LogPrint(0, "(%s)[%d,%d]+[%d,%d]<=>(%s)[%d,%d]+[%d,%d]", imgL->GetImageName(), pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3],
			imgR->GetImageName(), scr, srr, ncr, nrr);
		int ptNumL = imgL->ReadSiftFile(lpcstrSiftL, pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3]);
		if (ptNumL < 1)	{ LogPrint(0, "continue."); continue; }		
		ptNumL = PickPyrmSift((const PtName*)imgL->get_pt_name(0), (const PtLoc*)imgL->GetPtLocation(), (const PtFea*)imgL->GetSiftFeature(), (const SiftDes*)imgL->GetSiftDescriptors(),
			ptNumL, zoomRateL, nameL, locL, desL);
		LogPrint(0, "[%s]PyrmSift Num = %d.", imgL->GetImageName(), ptNumL);	if (ptNumL < 1)	{ LogPrint(0, "continue."); continue; }
		imgL->ClearPtBuffer();

		int ptNumR = imgR->ReadSiftFile(lpcstrSiftR, scr, srr, ncr, nrr);
		if (ptNumR < 1)	{ LogPrint(0, "continue."); continue; }
		ptNumR = PickPyrmSift((const PtName*)imgR->get_pt_name(0), (const PtLoc*)imgR->GetPtLocation(), (const PtFea*)imgR->GetSiftFeature(), (const SiftDes*)imgR->GetSiftDescriptors(),
			ptNumR, zoomRateR, nameR, locR, desR);
		LogPrint(0, "[%s]PyrmSift Num = %d.", imgR->GetImageName(), ptNumR);	if (ptNumR < 1)	{ LogPrint(0, "continue."); continue; }
		imgR->ClearPtBuffer();

		int(*match_buf)[2] = new int[ptNumL>ptNumR ? ptNumL : ptNumR][2];
		matcher.SetDescriptors(0, ptNumL, (const unsigned char*)desL.data());
		matcher.SetDescriptors(1, ptNumR, (const unsigned char*)desR.data());

		int match_num = matcher.GetSiftMatch(match_buf, siftDist, siftRatio, mutual_best_match);

		if (match_num > 0){
			MPt* mpt = new MPt[*ptSum + match_num];
			if (match_result) { memcpy(mpt, match_result, *ptSum*sizeof(MPt)); delete match_result; }
			match_result = mpt;
			MPt* pMPT = mpt + *ptSum;	int sz = 0;
			int i;	float lx0 = -99, ly0 = -99, rx0 = -99, ry0 = -99;
			for (i = 0; i < match_num; i++){
				float* xyl = (float*)(locL.data() + match_buf[i][0]);
				float* xyr = (float*)(locR.data() + match_buf[i][1]);
				pMPT->xl = *xyl;	pMPT->yl = *(xyl + 1);
				pMPT->xr = *xyr;	pMPT->yr = *(xyr + 1);
				if (fabs(pMPT->xl - lx0) < 1e-5 && fabs(pMPT->yl - ly0) < 1e-5
					&& fabs(pMPT->xr - rx0) < 1e-5 && fabs(pMPT->yr - ry0) < 1e-5) continue;
				lx0 = pMPT->xl; ly0 = pMPT->yl;	rx0 = pMPT->xr;	ry0 = pMPT->yr;				
				strcpy(pMPT->nameL, (const char*)(nameL.data() + match_buf[i][0]));
				strcpy(pMPT->nameR, (const char*)(nameR.data() + match_buf[i][1])); 
				pMPT++;	sz++;
			}
			LogPrint(0, "Remove repeat correspond points.Final Num= %d", sz);
			*ptSum += sz;
		}
		delete[] match_buf;
	}
	delete[] rc_split;
#endif
	LogPrint(0, "<----SiftMatch_Pyrm END[total num= %d][state= %s]---->\n", *ptSum, *ptSum<1 ? "Fail" : "OK");
	return match_result;
}
void FreeSiftMatchMem(MPt* pt) { delete pt; }

static bool ExtractSiftFile(CStdImage* img, int stBlkCol,int stBlkRow,int nBlkCols,int nBlkRows,
	const char* lpcstrTempDir,bool bUseExist,char* strSiftFile){
	char strSift[512],strTmpSift[512] = "",strTmpSiftInfo[512] = "";
	CStdImage::get_attach_file_path(img->GetImagePath(), AFT_SIFT, strSift);
	if (lpcstrTempDir) {
		sprintf(strTmpSift, "%s/%s_part.sift", lpcstrTempDir, img->GetImageName());
		strcpy(strTmpSiftInfo, strTmpSift);	strcat(strTmpSiftInfo, ".info");
	}

	if (bUseExist) { 
		if (img->CheckSiftFile(strSift, -1) > 0){
			LogPrint(0, "Use exist sift file %s.", strSift);
			strcpy(strSiftFile, strSift);  return true;
		}
		if (lpcstrTempDir) {
			int stBlkColT, stBlkRowT, nBlkColsT, nBlkRowsT;
			FILE* fp = fopen(strTmpSiftInfo, "r");
			if (fp){
				fscanf(fp, "%d%d%d%d", &stBlkColT, &stBlkRowT, &nBlkColsT, &nBlkRowsT);	fclose(fp);
				if (stBlkColT <= stBlkCol&&stBlkRowT <= stBlkRow
					&& (stBlkColT + nBlkColsT >= stBlkCol + nBlkCols) && (stBlkRowT + nBlkRowsT >= stBlkRow + nBlkRows)){
					if (img->CheckSiftFile(strTmpSift, -1) > 0){
						LogPrint(0, "Use exist sift file %s.", strTmpSift);
						strcpy(strSiftFile, strTmpSift);  return true;
					}
				}
			}
		}
	}
	
	float sz = 1.0f*img->GetRows()*img->GetCols();
	float s = (float)nBlkRows*nBlkCols / sz;
	if (s > 0.8 || sz<GlobalParam::g_extract_buffer_size ) {
		if (img->ExtractSift2PtFile(strSift, -1) < 1){
			LogPrint(0, "Can't extract sift feature from %s.",
				img->GetImageName());
			return false;
		}
		strcpy(strSiftFile, strSift);
	}
	else if (!lpcstrTempDir) {
		if (img->ExtractSift2PtFile(strSift, stBlkCol, stBlkRow, nBlkCols, nBlkRows, -1) < 1){
			LogPrint(0, "Can't extract sift feature from %s.",
				img->GetImageName());
			return false;
		}
		strcpy(strSiftFile, strSift);
	}
	else{
		
		if (img->ExtractSift2PtFile(strTmpSift, stBlkCol, stBlkRow, nBlkCols, nBlkRows, -1) < 1){
			LogPrint(0, "Can't extract sift feature from %s.",
				img->GetImageName());
			return false;
		}
		FILE* fp = fopen(strTmpSiftInfo, "w");
		if (fp){
			fprintf(fp, "%d\t%d\t%d\t%d\n", stBlkCol, stBlkRow, nBlkCols, nBlkRows);	fclose(fp);
		}
		strcpy(strSiftFile, strTmpSift);
	}
	return true;
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
		do {if (i + nExtraArg >= argc) \
{ char str[1024]; sprintf(str, "%s option requires %d argument(s)", argv[i], nExtraArg); Usage(str); }} while (0)

#ifndef EQUAL
#  if defined(WIN32) || defined(WIN32CE)
#    define STRCASECMP(a,b)         (stricmp(a,b))
#    define STRNCASECMP(a,b,n)      (strnicmp(a,b,n))
#  else
#    define STRCASECMP(a,b)         (strcasecmp(a,b))
#    define STRNCASECMP(a,b,n)      (strncasecmp(a,b,n))
#  endif
#  define EQUALN(a,b,n)           (STRNCASECMP(a,b,n)==0)
#  define EQUAL(a,b)              (STRCASECMP(a,b)==0)
#endif

int main(int argc, char* argv[])
{
#ifdef _DEBUG
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
//	_CrtSetBreakAlloc(269);
#endif
	char			strWsDir[512];
	bool			bGeoInfo = false;
	bool			bQuiet = false;
	bool			bUseExistSRTM = false, bUseExistSift = false;
	const char		*pszSource = NULL, *pszMatch = NULL, *pszWorkSpace = NULL, *pszMatchFile = NULL, *pszLogFile = NULL;
	const char		*pszSourceGeo = NULL, *pszMatchGeo = NULL,*pszDem = NULL;
	const char		*pszSourceID = NULL, *pszMatchID = NULL;
	int				rcL[4] = { -1 }, rcR[4] = { -1 };	
	char strErrMsg_rc[256];	
	sprintf(strErrMsg_rc, "The rectangle area of source dataset must be available"
		"( stCol < dataset_cols ; stRow < dataset_rows ; 0 < edCol ; 0 < edRow ; edCol-stCol > %d ;  edRow-stRow > %d ).", MIN_OVERLAP_LENGTH, MIN_OVERLAP_LENGTH);
	
	float			rangeG = RANGE_GEO, rangeB = RANGE_BLOCK, range = -1;
	bool			bRefine = false, bUseH = false, bUseF = false;
	float			H[3][3], F[3][3], hdistmax = -1, fdistmax = -1;
	float			siftDist = DEFAULT_SIFT_DIST, siftRatio = DEFAULT_SIFT_RATIO;
	int				mutual_best_match = 1;
	float			ratiomin_scale = -1, ratiomax_scale = -1, diffmin_ori = -1, diffmax_ori = -1;

	/* -------------------------------------------------------------------- */
	/*      Handle command line arguments.                                  */
	/* -------------------------------------------------------------------- */
	int i;
	for (i = 1; i < argc; i++)
	{
		if (EQUAL(argv[i], "--help"))
			Usage();
		else if (EQUAL(argv[i], "--long-usage"))
		{
			Usage(NULL, FALSE);
		}
		else if (EQUAL(argv[i], "-enable-gpu"))
		{
			GlobalParam::g_sift_extract_gpu = GlobalParam::g_sift_match_gpu = true;
		}
		else if (EQUAL(argv[i], "-enable-gpue"))
		{
			GlobalParam::g_sift_extract_gpu = true;
		}
		else if (EQUAL(argv[i], "-enable-gpum"))
		{
			GlobalParam::g_sift_match_gpu = true;
		}
		else if ((EQUAL(argv[i], "-ws") || EQUAL(argv[i], "-workspace")) && i < argc - 1)
		{
			pszWorkSpace = argv[++i];	strcpy(strWsDir, pszWorkSpace);	VERF_SLASH(strWsDir);
			if (IsExist(strWsDir) && !IsDir(strWsDir)) Usage("Workspace must be a directory.");
		}
		else if (EQUAL(argv[i], "-output") && i < argc - 1)
		{
			pszMatchFile = argv[++i];
			if (IsDir(pszMatchFile)) Usage("Output must not be a directory.");
		}
		else if (EQUAL(argv[i], "-log") && i < argc - 1)
		{
			pszLogFile = argv[++i];
			OpenLog(pszLogFile);
		}

		else if (EQUAL(argv[i], "-q") || EQUAL(argv[i], "-quiet"))
		{
			bQuiet = TRUE;
		}

		else if (EQUAL(argv[i], "-range"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			range = (float)atof(argv[++i]);
			if (range <= 0) Usage("Search range must > 0.");
		}
		else if (EQUAL(argv[i], "-rangeG") && i < argc - 1)
		{
			rangeG = atof(argv[++i]);
			if (rangeG<0) Usage("Geo-Search range must >= 0.");
		}
		else if (EQUAL(argv[i], "-rangeB") && i < argc - 1)
		{
			rangeB = atof(argv[++i]);
			if (rangeB<0) Usage("Block-Search range must >= 0.");
		}

		else if (EQUAL(argv[i], "-siftD") )
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			siftDist = atof(argv[++i]);
			if (siftDist < 0 || siftDist > 2 ) Usage("Maximum sift feature match distance must in [0,2].");
		}
		else if (EQUAL(argv[i], "-siftR"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			siftRatio = atof(argv[++i]);
			if (siftRatio < 0 || siftRatio > 1) Usage("Maximum sift feature match distance ratio must in [0,1].");
		}

		else if (EQUAL(argv[i], "-ue"))
		{
			bUseExistSRTM = bUseExistSift = true;
		}
		else if (EQUALN(argv[i], "-refine",7))
		{
			bRefine = true;
			int l = strlen(argv[i]);
			const char* pS = argv[i] + 7;
			char strE[100];	sprintf(strE,"Unknown option name '%s' for 'refine'.", argv[i]);
			if (l == 8) {
				if (*pS == 'h' || *pS == 'H') bUseH = true;
				else if (*pS == 'f' || *pS == 'F') bUseF = true;
				else Usage(strE);
			}else if (l==9){
				if (EQUAL(pS, "FH") || EQUAL(pS, "HF")){
					bUseF = bUseH = true;
				}
				else Usage(strE);
			}
			else if(l > 9 ) Usage(strE);
		}
		else if (EQUAL(argv[i], "-no-mutual"))
		{
			mutual_best_match = 0;
		}
		
		else if (EQUAL(argv[i], "-scale"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			ratiomin_scale = (float)atof(argv[++i]);
			if (ratiomin_scale < 0) Usage("The scale ratio of Sift feature must >= 0.");
			char* endptr = NULL;
			if (argv[i + 1] != NULL
				&& (strtof(argv[i + 1], &endptr) != 0.0 || argv[i + 1][0] == '0'))
			{
				/* Check that last argument is really a number and not a filename */
				/* looking like a number (see ticket #863) */
				if (endptr && *endptr == 0)
				{
					ratiomax_scale = (float)atof(argv[++i]);
					if (ratiomax_scale < 0) Usage("The scale ratio of Sift feature must >= 0.");
					if (ratiomax_scale < ratiomin_scale) Usage("The Maximum of scale ration must larger than the  minimum.");
				}
			}
			
		}
		else if (EQUAL(argv[i], "-orientation"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			diffmin_ori = atof(argv[++i]);
			if (argv[i][strlen(argv[i]) - 1] != 'r') { diffmin_ori = diffmin_ori*PI / 180; }
			if (diffmin_ori < 0 || diffmin_ori > PI ) Usage("The orientation difference of Sift feature must in [0бу,180бу].");
			char* endptr = NULL;
			if (argv[i + 1] != NULL
				&& (strtof(argv[i + 1], &endptr) != 0.0 || argv[i + 1][0] == '0'))
			{
				/* Check that last argument is really a number and not a filename */
				/* looking like a number (see ticket #863) */
				if (endptr && *endptr == 0)
				{
					diffmax_ori = atof(argv[++i]);
					if (argv[i][strlen(argv[i]) - 1] != 'r') { diffmax_ori = diffmax_ori*PI / 180; }
					if (diffmax_ori < 0 || diffmax_ori > PI) Usage("The orientation difference of Sift feature must in [0бу,180бу].");
					if (diffmax_ori < diffmin_ori) Usage("The Maximum of orientation difference must larger than the  minimum.");
				}
			}
		}
		else if (EQUAL(argv[i], "-ids")){
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			pszSourceID = argv[++i];
		}
		else if (EQUAL(argv[i], "-idd")){
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			pszMatchID = argv[++i];
		}

		else if (EQUAL(argv[i], "-geo"))
		{
			bGeoInfo = true;
		}
		else if (EQUAL(argv[i], "-geos"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			pszSourceGeo = argv[++i];
			bGeoInfo = true;
		}
		else if (EQUAL(argv[i], "-geod"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			pszMatchGeo = argv[++i];
			bGeoInfo = true;
		}
		else if (EQUAL(argv[i], "-rc") || EQUAL(argv[i], "-rcs"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(4);
			int t = atoi(argv[++i]);	rcL[0] = t > 0 ? t : 0;
			t = atoi(argv[++i]);	rcL[1] = t > 0 ? t : 0;
			t = atoi(argv[++i]);	if (t <= 0 || t - rcL[0] <= MIN_OVERLAP_LENGTH) Usage(strErrMsg_rc);	rcL[2] = t;
			t = atoi(argv[++i]);	if (t <= 0 || t - rcL[1] <= MIN_OVERLAP_LENGTH) Usage(strErrMsg_rc);	rcL[3] = t;
		}
		else if (EQUAL(argv[i], "-rcd"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(4);
			int t = atoi(argv[++i]);	rcR[0] = t > 0 ? t : 0;
			t = atoi(argv[++i]);	rcR[1] = t > 0 ? t : 0;
			t = atoi(argv[++i]);	if (t <= 0 || t - rcR[0] <= MIN_OVERLAP_LENGTH) Usage(strErrMsg_rc);	rcR[2] = t;
			t = atoi(argv[++i]);	if (t <= 0 || t - rcR[1] <= MIN_OVERLAP_LENGTH) Usage(strErrMsg_rc);	rcR[3] = t;
		}
		else if (EQUAL(argv[i], "-dem") && i < argc - 1)
		{
			pszDem = argv[++i];
		}

		else if (argv[i][0] == '-')
		{
			char strE[512];	sprintf(strE, "Unknown option name '%s'", argv[i]);
			Usage(strE);
		}
		else if (pszSource == NULL)
		{
			pszSource = argv[i];
		}
		else if (pszMatch == NULL)
		{
			pszMatch = argv[i];
		}

		else
		{
			Usage("Too many command options.");
		}
	}

	if (rangeG > rangeB)	Usage("Block-Search range must not less than Geo-Search range.");

#define RC_IMAGE(rc,img)	\
	{ \
	if (rc[0] < 0) { rc[0] = 0; rc[1] = 0;	rc[2] = img.GetCols(); rc[3] = img.GetRows(); }	\
	else if (rc[0] > img.GetCols() - 1 || rc[1] > img.GetRows() - 1) { LogPrint(0, strErrMsg_rc); return 1; }\
	else{\
	if (rc[2] > img.GetCols() - 1) rc[2] = img.GetCols() - 1; \
	if (rc[3] > img.GetRows() - 1) rc[3] = img.GetRows() - 1; \
		rc[2] = rc[2] - rc[0] + 1; rc[3] = rc[3] - rc[1] + 1; if (rc[2] <= MIN_OVERLAP_LENGTH || rc[3] <= MIN_OVERLAP_LENGTH) { LogPrint(0, strErrMsg_rc); return 1; } \
	} }
	
	CGeoImage imgL, imgR;
	int ret_code = 0;
		
	if (pszWorkSpace)	{
		if (!IsExist(strWsDir) && !CreateDir(strWsDir)) Usage("Fail to create workspace.");
		CStdImage::set_workspace_path(strWsDir);
	}
	if (pszMatch == NULL)
	{
		if (pszSource == NULL)
			Usage("No source dataset specified.");
		else{
			if (!imgL.Open(pszSource)){
				LogPrint(0, "Can't open image %s.",
					pszSource);					
				exit(1);
			}
			char strPath[512];
			RC_IMAGE(rcL, imgL);
			if (!ExtractSiftFile(&imgL,rcL[0],rcL[1],rcL[2],rcL[3],NULL,bUseExistSift,strPath) ){
				LogPrint(0,
					"Can't extract sift feature from %s.",
					pszSource);
				return 1;
			}
			return 0;
		}
	}
	
	if (!pszMatchFile){
		if (!pszWorkSpace) Usage("workspace and output are null.At least one of them must be set.");		
	}
	if (!imgL.Open(pszSource)){
		LogPrint(0,
			"Can't open image %s.",
			pszSource);
		exit(1);
	}
	if (!imgR.Open(pszMatch)){
		LogPrint(0,
			"Can't open image %s.",
			pszMatch);
		exit(1);
	}
	CStdDem	dem;
	if (pszDem && !dem.Load4File(pszDem)){
		LogPrint(0,
			"Can't open DEM %s.",
			pszDem);
		exit(1);
	}
	
#define EXCHANGEIMAGE	\
	{	\
		CGeoImage img;	img = imgL;	imgL = imgR; imgR = img;	bImgExchanged = !bImgExchanged;	\
	}	\

	bool				bImgExchanged = false;
	float				fGsd_ratio = 1;	
	CImagesCoorT*		pCvt = NULL;
	CAltitude*			pAlt = NULL;
	if (bGeoInfo){
		if (pszSourceGeo){
			if (!imgL.InitRpcImage(pszSource, pszSourceGeo)){
				LogPrint(0,
					"%s is unsupported GeoFile.",
					pszSourceGeo);
				exit(1);
			}
		}
		else if (!imgL.InitTFWImage(pszSource)){
			LogPrint(0,
				"%s is not GeoImage.",
				pszSource);
			exit(1);
		}
		if (pszMatchGeo){
			if (!imgR.InitRpcImage(pszMatch, pszMatchGeo)){
				LogPrint(0,
					"%s is unsupported GeoFile.",
					pszMatchGeo);
				exit(1);
			}
		}
		else if (!imgR.InitTFWImage(pszMatch)){
			LogPrint(0,
				"%s is not GeoImage.",
				pszMatch);
			exit(1);
		}
		if (pszSourceID) imgL.SetImageID(pszSourceID);	if (pszMatchID) imgR.SetImageID(pszMatchID);
		fGsd_ratio = imgR.GetGsd() / imgL.GetGsd();
		if (fGsd_ratio < 0.9){
			EXCHANGEIMAGE;
			fGsd_ratio = 1 / fGsd_ratio;
		}
		/*
		if (!pszDem){
			char strSRTMPath[512];	CStdImage::get_attach_file_path(imgL.GetImagePath(), AFT_SRTM, strSRTMPath);
			if (bUseExistSRTM&&CDUImage::CheckSRTMFile(strSRTMPath) && dem.Load4File(strSRTMPath)) {
				LogPrint(0, "Use exist SRTM file %s.", strSRTMPath);		pszDem = strSRTMPath;
			}
			else if (CDUImage::SRTM(&imgL, strSRTMPath) && dem.Load4File(strSRTMPath)) {
				pszDem = strSRTMPath;
			}
			bUseExistSRTM = true;
		}
		if (pszDem){
			CDemAltitude* pDem = new CDemAltitude;
			pDem->init(&dem);
			pAlt = pDem;
		}
		else{
			CAvrAltitude* pAvr = new CAvrAltitude;
			pAvr->init(0, 1000);
			pAlt = pAvr;
		}*/
		if (pszDem){
			CDemAltitude* pDem = new CDemAltitude;
			pDem->init(&dem);
			pAlt = pDem;
		}

		CGeoImagesCoorT* cvt = new CGeoImagesCoorT;		
		if (!cvt->init(&imgL, &imgR, pAlt)){
			LogPrint(0, "Fail to initialize the transformation between two images.");
			delete cvt;
			goto end1;
		}
		pCvt = cvt;
	}
	else{
		if (pszSourceID) imgL.SetImageID(pszSourceID);	if (pszMatchID) imgR.SetImageID(pszMatchID);
		if (ratiomax_scale>0 && (ratiomax_scale<1 || ratiomin_scale>1)){
			if (ratiomax_scale < 1) fGsd_ratio = 1 / ratiomax_scale; 
			else{
				EXCHANGEIMAGE;
				fGsd_ratio = ratiomin_scale;				
			}
		}
		CStdImagesCoorT* cvt = new CStdImagesCoorT;
		cvt->init(&imgL, &imgR);
		pCvt = cvt;
		range = -1;
	}
	pCvt->set_range((int)range);

	char strModel[512];	char strTempDir[512];
	int stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL;
	int stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR;
	char strSiftFileL[512], strSiftFileR[512];
	MPt* match_result = NULL;		int match_num;

	if (pszMatchFile) strcpy(strModel, pszMatchFile);
	if (pszWorkSpace){
		sprintf(strTempDir, "%s/temp/%s_%s", CStdImage::get_workspace_path(),imgL.GetImageID(),imgR.GetImageID());	CreateDir(strTempDir);
	}
	else{		
		strcpy(strTempDir, pszMatchFile);	char* pS = strrchr(strTempDir, '.');	if (!pS) pS = strTempDir + strlen(strTempDir);
		strcpy(pS, "_temp");	CreateDir(strTempDir);
	}
refine:
	
	if (fGsd_ratio < 0.9 ){
		EXCHANGEIMAGE;
		CImagesCoorT* cvt = pCvt->CreateInvertT(false);
		delete[] pCvt;	pCvt = cvt;
		fGsd_ratio = 1 / fGsd_ratio;
	}

	if (!pszMatchFile){
		sprintf(strModel, "%s/%s=%s.model", CStdImage::get_workspace_path(), imgL.GetImageName(), imgR.GetImageName());		
	}
	
	if (!pCvt->overlap_on_left(&stBlkColL, &stBlkRowL, &nBlkColsL, &nBlkRowsL)){
		LogPrint(0,
			"No overlap between the two images.");
		ret_code = 2;
		goto end2;
	}
	nBlkColsL = nBlkColsL - stBlkColL + 1;	nBlkRowsL = nBlkRowsL - stBlkRowL + 1;
	if (!pCvt->overlap_on_right(&stBlkColR, &stBlkRowR, &nBlkColsR, &nBlkRowsR)){
		LogPrint(0,
			"No overlap between the two images.");
		ret_code = 2;
		goto end2;
	}
	nBlkColsR = nBlkColsR - stBlkColR + 1;	nBlkRowsR = nBlkRowsR - stBlkRowR + 1;
	LogPrint(0, "\nOverlap=\t(%s)[%d,%d]+[%d,%d]<=>(%s)[%d,%d]+[%d,%d]\n",
		imgL.GetImageName(), stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, imgR.GetImageName(), stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);

	if (nBlkColsL < MIN_OVERLAP_LENGTH || nBlkRowsL < MIN_OVERLAP_LENGTH || nBlkColsR < MIN_OVERLAP_LENGTH || nBlkRowsR < MIN_OVERLAP_LENGTH){
		LogPrint(0, "The overlap on row/col is too small (less than %d).", MIN_OVERLAP_LENGTH);
		ret_code = 2;
		goto end2;
	}
	int search_range = nBlkColsR>nBlkRowsR ? nBlkColsR : nBlkRowsR;
	if (search_range < rangeB) search_range = rangeB;
	if (range < 0 || range > search_range){		
		range = search_range;
	}

#define SAVE_MATCH_FILE(strModel,imgL,imgR)	\
	{	int nb[7] = { sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt), 0 };	\
	CMatchFile::WriteModelFile(strModel, \
	imgL.GetImageID(), imgL.GetImagePath(), \
	imgR.GetImageID(), imgR.GetImagePath(), \
	match_result->nameL, &match_result->xl, &match_result->yl, match_result->nameR, &match_result->xr, &match_result->yr, NULL, match_num, nb);	\
	LogPrint(0, "Save match result to %s.",strModel); }\

#ifndef SIFT_PYRM_IMAGE
	{
		if (!ExtractSiftFile(&imgL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, strTempDir, bUseExistSift, strSiftFileL)) goto end2;
		if (!ExtractSiftFile(&imgR, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, strTempDir, bUseExistSift, strSiftFileR)) goto end2;
		bUseExistSift = true;
	}
#endif

	if (range > rangeB){		
		float zoomRate = rangeB / range;	float length = (float)sqrt(GlobalParam::g_match_buffer_size)*0.25f;
		if (zoomRate*nBlkColsR < length) zoomRate = length / nBlkColsR;	if (zoomRate*nBlkRowsR < length) zoomRate = length / nBlkRowsR;
		if (zoomRate*nBlkColsL / fGsd_ratio < length) zoomRate = length*fGsd_ratio / nBlkColsL;
		if (zoomRate*nBlkRowsL / fGsd_ratio < length) zoomRate = length*fGsd_ratio / nBlkRowsL;
		if (zoomRate < 1) {
			match_result = SiftMatch_Pyrm(&match_num,
				&imgL, strSiftFileL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL,
				&imgR, strSiftFileR, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR,
				pCvt, zoomRate/fGsd_ratio,zoomRate, siftDist, siftRatio, mutual_best_match, strTempDir);
			SAVE_MATCH_FILE(strModel, imgL, imgR);
			if (match_num > 1) {
				char strPath[512];	
				sprintf(strPath, "%s/%s=%s_range%d_zoom%.f.model", strTempDir, imgL.GetImageName(), imgR.GetImageName(), (int)range, 1 / zoomRate);
				CopyFile(strModel, strPath, FALSE);
				int nBufferSpace[4] = { sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt) };
				match_num = pCvt->update(&match_result->xl, &match_result->yl, &match_result->xr, &match_result->yr, nBufferSpace,true, NULL, match_num, NULL, (rangeG + rangeB) / 2);
				if (match_num > 0){
					FreeSiftMatchMem(match_result);
					range = rangeB;	pCvt->set_range(range);	fGsd_ratio = 1 / pCvt->get_gsd_ratio();
					goto refine;
				}
			}
			FreeSiftMatchMem(match_result);
		}
		else range = rangeB;
	}
	if (range <= rangeB){	

#ifdef SIFT_PYRM_IMAGE		
		if (!ExtractSiftFile(&imgL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, strTempDir, bUseExistSift, strSiftFileL)) goto end2;
		if (!ExtractSiftFile(&imgR, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, strTempDir, bUseExistSift, strSiftFileR)) goto end2;
		bUseExistSift = true;		
#endif
		if (range > rangeG){
			match_result = SiftMatch_Block(&match_num,
				&imgL, strSiftFileL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL,
				&imgR, strSiftFileR, pCvt, siftDist, siftRatio, mutual_best_match);
			SAVE_MATCH_FILE(strModel, imgL, imgR);

			if (bRefine && match_num > 1) {
				char strPath[512];	sprintf(strPath, "%s/%s=%s_range%d.model", strTempDir, imgL.GetImageName(), imgR.GetImageName(), (int)range);
				CopyFile(strModel, strPath, FALSE);
				int nBufferSpace[4] = { sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt) };
				if (bUseH || bUseF){
					int h_left_num = 0;	int f_left_num = 0;
					char srange[100];	sprintf(srange, "%f", rangeG);
					char * argv[] = { "method", "ransacA", "bufferType", "byte", "DistanceThreshold", srange };
					int argc = sizeof(argv) / sizeof(char*);
					h_left_num = GeoTransform::EstimateXformMatrix(&match_result->xl, &match_result->yl, &match_result->xr, &match_result->yr, nBufferSpace, match_num, NULL, H[0], argc, argv);
					if (!h_left_num) {
						LogPrint(0, "Fail to calculate the homography matrix between two image.");
						match_num = 0;
					}
					else {
						delete pCvt;
						CStdImagesCoorT* cvt = new CStdImagesCoorT;	cvt->init(&imgL, &imgR, H);
						pCvt = cvt;						
						if (bUseF){
							f_left_num = GeoTransform::EstimateFundamentalMatrix(&match_result->xl, &match_result->yl, &match_result->xr, &match_result->yr, nBufferSpace,match_num, NULL, F[0], argc, argv);
							if (!f_left_num) {
								LogPrint(0, "Fail to calculate the fundamental matrix between two image.");
								match_num = 0;
							}
						}
					}
				}
				else{
					match_num = pCvt->update(&match_result->xl, &match_result->yl, &match_result->xr, &match_result->yr, nBufferSpace,true, NULL, match_num, NULL, rangeG);
				}
				if (match_num > 0){
					range = (int)rangeG;		pCvt->set_range(range);	fGsd_ratio = 1 / pCvt->get_gsd_ratio();
					FreeSiftMatchMem(match_result);
					goto refine;
				}
			}
			FreeSiftMatchMem(match_result);
		}
		if (range <= rangeG){
			if (bUseH) hdistmax = (float)pCvt->get_range();	if (bUseF) fdistmax = (float)pCvt->get_range();
			
			if ((ratiomin_scale >= 0 && ratiomax_scale < 0) || (diffmin_ori >= 0 && diffmax_ori < 0)){
				float affineMatrix[6];
				if (pCvt->calculate_affine_matrix(stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, stBlkColR, stBlkRowR, affineMatrix)){
					float tmp = 0;
					if (ratiomin_scale >= 0 && ratiomax_scale < 0){
						tmp = sqrt((affineMatrix[0] * affineMatrix[0] + affineMatrix[1] * affineMatrix[1] + affineMatrix[3] * affineMatrix[3] + affineMatrix[4] * affineMatrix[4]) / 2);
						if (bImgExchanged&&tmp>0){
							tmp = 1 / tmp;														
						}
						ratiomax_scale = tmp + ratiomin_scale;	ratiomin_scale = tmp - ratiomin_scale; 
						if (ratiomin_scale < 0) ratiomin_scale = 0;
					}
					if (diffmin_ori >= 0 && diffmax_ori < 0){
						tmp = (float)fabs(atan2(affineMatrix[3], affineMatrix[0]));
						diffmax_ori = tmp + diffmin_ori;	if (diffmax_ori > PI) diffmax_ori = (float)PI;
						diffmin_ori = tmp - diffmin_ori;	if (diffmin_ori < 0) diffmin_ori = 0;
					}					
				}				
			}
			if (ratiomax_scale > 0 && bImgExchanged){
				float t = 1 / ratiomax_scale;	ratiomax_scale = ratiomin_scale < 1e-6f ? 1e+7f : 1 / ratiomin_scale;	ratiomin_scale = t;
			}
			match_result = SiftMatch_Geo(&match_num,
				&imgL, strSiftFileL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL,
				&imgR, strSiftFileR, pCvt, bUseH ? H : NULL, hdistmax, bUseF ? F : NULL, fdistmax,
				ratiomin_scale, ratiomax_scale, diffmin_ori, diffmax_ori,
				siftDist, siftRatio, mutual_best_match);
			SAVE_MATCH_FILE(strModel, imgL, imgR);
			FreeSiftMatchMem(match_result);
		}
	}
	
end2:
	if (!pszWorkSpace) RemoveDir(strTempDir);	
end1:
	if (pCvt) delete[] pCvt;
	if (pAlt) delete[] pAlt;
	return ret_code;
}
