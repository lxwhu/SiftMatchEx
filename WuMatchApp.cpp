
#include "PrjLog.hpp"

#ifdef _NO_OPENGL
#include "DPGRID/WuZMch-v1.0.h"
#else
#include "DPGRID/WuZMch-v2.0.h"
#endif

#include "../AtMatch/GlobalUtil_AtMatch.h"
#include "../AtMatch/StdImage.h"
#include "../AtMatch/MatchFile.h"
#include "../AtMatch/ImagesCoorAssociate.hpp"
#include "algorithms/Pretreatment.h"
#include "algorithms/Transformation.h"

#define PI							3.1415926
#define RELEASE_NAME				"1.0.0"	
#define RELEASE_DATE				"2016/4/24"

#define RANGE_BLOCK					500
#define RANGE_GEO					100
#define MIN_OVERLAP_LENGTH			64	

const char* WuMatchVersionInfo(){
	return "WuMatch "RELEASE_NAME
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
	printf("Usage: WuMatchEx [--help-general] [--long-usage]\n"
		"       [-ws workspace] [-output match_result_file] [-log log_file]\n"
		"       [-range(G/B) search_range]\n"
		"       [-sift{D/R} distance/ratio]\n"
		"       [-ue] [-refine(h/f)] [-no-mutual]\n"
		"       [-scale ratio_min ratio_max] [-orientation diff_min(r) diff_max(r)]\n"
		"       src_dataset dst_dataset\n"
		"       [-id{s/d} image_id]\n"
		"       [-geo(s/d) (geo_file)]\n"
		"       [-dem dem_file]\n");

	if (!bShort)
	{
		printf("\n%s\n\n", WuMatchVersionInfo());
		printf("The following format drivers are configured :\n");
		const char* strDriver[] = { "GTiff(*.tif,*tiff)", "HFA(*.img)", "PNG(*.png)", "JPEG(*.jpg)", "PCIDSK(*.pix)", "BMP(*.bmp)" };
		for (int iDr = 0; iDr < sizeof(strDriver) / sizeof(const char*); iDr++)
		{
			printf("%s;\t", strDriver[iDr]);
		}
		printf("\nThe following format GeoFile are configured :\n");
		printf("RPC;\n");
		printf("Default value: siftDist(%.2f) siftRatio(%.2f) GeoRange(%d) BlockRange(%d)\n", SIFT_DIST, SIFT_RATIO, RANGE_GEO, RANGE_BLOCK);
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

MPt*		WuMatch_Block(
	int* ptSum,
	CStdImage* imgL, const char* lpcstrSiftL, int stBlkColL, int stBlkRowL, int nBlkColsL, int nBlkRowsL,
	CStdImage* imgR, const char* lpcstrSiftR,
	CImagesCoorT*		pCvt, float range,
	float siftDist, float siftRatio, int mutual_best_match
	){
	LogPrint(0, "\n<----WuMatch_Block starting[range= %.2f]---->", range);
	*ptSum = 0;

	int nColNum, nRowNum;
	int* rc_split = Pretreatment::split_image(nColNum, nRowNum, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, 8,
		GlobalParam::g_match_buffer_size, Pretreatment::RECT_FIRST);
	int blk_num = nColNum*nRowNum;

	MPt* match_result = NULL;
	int *pSplitRc = rc_split;
	for (int i = 0; i < blk_num; i++, pSplitRc += 4){
		LogPrint(0, ">>>>block[%d/%d]", i + 1, blk_num);
		int stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR;
		if (!pCvt->overlap_on_right(pSplitRc[0], pSplitRc[1], pSplitRc[0] + pSplitRc[2], pSplitRc[1] + pSplitRc[3],
			(int)range, (int)range, &stBlkColR, &stBlkRowR, &nBlkColsR, &nBlkRowsR)){
			LogPrint(0, "no coverage.continue.");
			continue;
		}
		nBlkColsR = nBlkColsR - stBlkColR + 1;	nBlkRowsR = nBlkRowsR - stBlkRowR + 1;
		LogPrint(0, "(%s)[%d,%d]+[%d,%d]<=>(%s)[%d,%d]+[%d,%d]", imgL->GetImageName(), pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3], imgR->GetImageName(), stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);
		int ptNumL = imgL->ReadSiftFile(lpcstrSiftL, pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3]);
		if (ptNumL < 1)	{ LogPrint(0, "continue."); continue; }
		int ptNumR = imgR->ReadSiftFile(lpcstrSiftR, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);
		if (ptNumR < 1)	{ LogPrint(0, "continue."); continue; }

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
	LogPrint(0, "<----SiftMatch_Block END[total num= %d][state= %s]---->\n", *ptSum, *ptSum < 1 ? "Fail" : "OK");
	return match_result;
}

MPt*		WuMatch_Pyrm(
	int* ptSum,
	CStdImage* imgL, int stBlkColL, int stBlkRowL, int nBlkColsL, int nBlkRowsL,
	CStdImage* imgR, int stBlkColR, int stBlkRowR, int nBlkColsR, int nBlkRowsR,
	CImagesCoorT*		pCvt, float range, float zoomRate,
	float siftDist, float siftRatio, int mutual_best_match,
	const char* lpcstrTempDir
	){
	LogPrint(0, "\n<----SiftMatch_Pyrm starting[range= %.2f;zoom= %f]---->", range, zoomRate);
	*ptSum = 0;
	CSiftMatch matcher;
	if (!matcher.InitEnvi(GlobalParam::g_sift_match_gpu)) return NULL;

	CStdImage imgPyrmL, imgPyrmR;
	char	strSiftPyrmL[512], strSiftPyrmR[512];
	int nPyrmColsL, nPyrmRowsL, nPyrmColsR, nPyrmRowsR;
	float gsd_ratio = pCvt->get_gsd_ratio();
	float zoomL = zoomRate*gsd_ratio;	float zoomR = zoomRate;

	char strFile[512];	strcpy(strFile, lpcstrTempDir);
	char* pS = strFile + strlen(strFile);
	{
		sprintf(pS, "/%s_%d_%d_%d_%d_%.f.tif", imgL->GetImageName(),
			stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, 1 / zoomL);
		BYTE* pBuf = NULL;	FILE* fp = NULL;
		pBuf = new BYTE[int(nBlkColsL*nBlkRowsL*zoomL) + 8];	if (!pBuf) return NULL;
		imgL->ReadGray8(pBuf, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, zoomL, &nPyrmColsL, &nPyrmRowsL);
		if (!SaveImageFile(strFile, pBuf, nPyrmColsL, nPyrmRowsL, 1) || !imgPyrmL.Open(strFile))
		{
			delete pBuf;
			return NULL;
		}
		delete pBuf;
		strcpy(strrchr(strFile, '.'), ".info");
		fp = fopen(strFile, "w"); if (!fp) return NULL;	fprintf(fp, "%d\t%d\t%d\t%d\t%f", stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, zoomL); fclose(fp);
		strcpy(strrchr(strFile, '.'), ".sift");	strcpy(strSiftPyrmL, strFile);
		if (imgPyrmL.ExtractSift2PtFile(strSiftPyrmL, -1) < 1) return NULL;

		sprintf(pS, "/%s_%d_%d_%d_%d_%.f.tif", imgR->GetImageName(),
			stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, 1 / zoomR);
		pBuf = new BYTE[int(nBlkColsR*nBlkRowsR*zoomR) + 8];	if (!pBuf) return NULL;
		imgR->ReadGray8(pBuf, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, zoomR, &nPyrmColsR, &nPyrmRowsR);
		if (!SaveImageFile(strFile, pBuf, nPyrmColsR, nPyrmRowsR, 1) || !imgPyrmR.Open(strFile))
		{
			delete pBuf;
			return NULL;
		}
		delete pBuf;
		strcpy(strrchr(strFile, '.'), ".info");
		fp = fopen(strFile, "w"); if (!fp) return NULL;	fprintf(fp, "%d\t%d\t%d\t%d\t%f", stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, zoomR); fclose(fp);
		strcpy(strrchr(strFile, '.'), ".sift");	strcpy(strSiftPyrmR, strFile);
		if (imgPyrmR.ExtractSift2PtFile(strSiftPyrmR, -1) < 1) return NULL;
	}

	int nColNum, nRowNum;
	int* rc_split = Pretreatment::split_image(nColNum, nRowNum, 0, 0, nPyrmColsL, nPyrmRowsL, 8,
		GlobalParam::g_match_buffer_size, Pretreatment::RECT_FIRST);
	int blk_num = nColNum*nRowNum;

	MPt* match_result = NULL;
	int *pSplitRc = rc_split;
	for (int i = 0; i < blk_num; i++, pSplitRc += 4){
		LogPrint(0, ">>>>Pyrmblock<%.6f>[%d/%d]", zoomRate, i + 1, blk_num);
		int scr, srr, ncr, nrr;
		if (!pCvt->overlap_on_right(stBlkColL + pSplitRc[0] / zoomL, stBlkRowL + pSplitRc[1] / zoomL,
			stBlkColL + (pSplitRc[0] + pSplitRc[2]) / zoomL, stBlkRowL + (pSplitRc[1] + pSplitRc[3]) / zoomL,
			(int)range, (int)range, &scr, &srr, &ncr, &nrr)){
			LogPrint(0, "no coverage.continue.");
			continue;
		}
		ncr = int((ncr - scr + 1)*zoomR);	nrr = int((nrr - srr + 1)*zoomR);
		scr -= stBlkColR;	srr -= stBlkRowR;

		// 		LogPrint(0, "left[%d,%d]+[%d,%d]<=>right[%d,%d]+[%d,%d]", 
		// 			pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3], stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR);

		int ptNumL = imgPyrmL.ReadSiftFile(strSiftPyrmL, pSplitRc[0], pSplitRc[1], pSplitRc[2], pSplitRc[3]);
		if (ptNumL < 1)	{ LogPrint(0, "continue."); continue; }
		int ptNumR = imgPyrmR.ReadSiftFile(strSiftPyrmR, scr, srr, ncr, nrr);
		if (ptNumR < 1)	{ LogPrint(0, "continue."); continue; }

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
				pMPT->xl = pMPT->xl / zoomL + stBlkColL;	pMPT->yl = pMPT->yl / zoomL + stBlkRowL;
				pMPT->xr = pMPT->xr / zoomR + stBlkColR;	pMPT->yr = pMPT->yr / zoomR + stBlkRowR;
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
	LogPrint(0, "<----SiftMatch_Pyrm END[total num= %d][state= %s]---->\n", *ptSum, *ptSum<1 ? "Fail" : "OK");
	return match_result;
}
void FreeSiftMatchMem(MPt* pt) { delete pt; }

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
do {
if (i + nExtraArg >= argc) \
{ char str[1024]; sprintf(str, "%s option requires %d argument(s)", argv[i], nExtraArg); Usage(str); }
} while (0)

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
	bool			bGeoInfo = false;
	bool			bQuiet = false;
	bool			bUseExistSRTM = false, bUseExistSift = false;
	const char		*pszSource = NULL, *pszMatch = NULL, *pszWorkSpace = NULL, *pszMatchFile = NULL, *pszLogFile = NULL;
	const char		*pszSourceGeo = NULL, *pszMatchGeo = NULL, *pszDem = NULL;
	const char		*pszSourceID = NULL, *pszMatchID = NULL;

	float			rangeG = RANGE_GEO, rangeB = RANGE_BLOCK, range = -1;
	bool			bRefine = false, bUseH = false, bUseF = false;
	float			H[3][3], F[3][3], hdistmax = -1, fdistmax = -1;
	float			siftDist = SIFT_DIST, siftRatio = SIFT_RATIO;
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
			pszWorkSpace = argv[++i];
			if (IsExist(pszWorkSpace) && !IsDir(pszWorkSpace)) Usage("Workspace must be a directory.");
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
			if (range < 0) Usage("Search range must >= 0.");
		}
		else if (EQUAL(argv[i], "-rangeG") && i < argc - 1)
		{
			rangeG = atof(argv[++i]);
			if (rangeG < 0) Usage("Geo-Search range must >= 0.");
		}
		else if (EQUAL(argv[i], "-rangeB") && i < argc - 1)
		{
			rangeB = atof(argv[++i]);
			if (rangeB < 0) Usage("Block-Search range must >= 0.");
		}

		else if (EQUAL(argv[i], "-siftD"))
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			siftDist = atof(argv[++i]);
			if (siftDist < 0 || siftDist > 2) Usage("Maximum sift feature match distance must in [0,2].");
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
		else if (EQUALN(argv[i], "-refine", 7))
		{
			bRefine = true;
			int l = strlen(argv[i]);
			const char* pS = argv[i] + 7;
			char strE[100];	sprintf(strE, "Unknown option name '%s' for 'refine'.", argv[i]);
			if (l == 8) {
				if (*pS == 'h' || *pS == 'H') bUseH = true;
				else if (*pS == 'f' || *pS == 'F') bUseF = true;
				else Usage(strE);
			}
			else if (l == 9){
				if (EQUAL(pS, "FH") || EQUAL(pS, "HF")){
					bUseF = bUseH = true;
				}
				else Usage(strE);
			}
			else if (l > 9) Usage(strE);
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
			if (diffmin_ori < 0 || diffmin_ori > PI) Usage("The orientation difference of Sift feature must in [0бу,180бу].");
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

	CGeoImage imgL, imgR;
	int ret_code = 0;
	if (pszWorkSpace)	{
		if (!IsExist(pszWorkSpace) && !CreateDir(pszWorkSpace)) Usage("Fail to create workspace.");
		CStdImage::set_workspace_path(pszWorkSpace);
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
			if (ExtractSiftFile(&imgL, 0, 0, imgL.GetCols(), imgL.GetRows(), NULL, bUseExistSift, strPath) < 1){
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
			CGeoImage img;
			img = imgL;	imgL = imgR; imgR = img;
			fGsd_ratio = 1 / fGsd_ratio;
		}
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
		}

		CGeoImagesCoorT* cvt = new CGeoImagesCoorT;
		if (!cvt->init(&imgL, &imgR, pAlt)){
			LogPrint(0, "Fail to initialize the transformation between two images.");
			delete cvt;
			goto end1;
		}
		pCvt = cvt;
		if (range < 0){ range = imgR.GetCols()>imgR.GetRows() ? imgR.GetCols() : imgR.GetRows(); }
	}
	else{
		if (pszSourceID) imgL.SetImageID(pszSourceID);	if (pszMatchID) imgR.SetImageID(pszMatchID);
		CStdImagesCoorT* cvt = new CStdImagesCoorT;
		cvt->init(&imgL, &imgR);
		pCvt = cvt;
		if (range < 0){ range = imgR.GetCols()>imgR.GetRows() ? imgR.GetCols() : imgR.GetRows(); }
	}

	char strModel[512];	char strTempDir[512];
	int stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL;
	int stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR;
	char strSiftFileL[512], strSiftFileR[512];
	MPt* match_result = NULL;		int match_num;

	if (!pszMatchFile){
		sprintf(strTempDir, "%s/TEMP", CStdImage::get_workspace_path());	CreateDir(strTempDir);
	}
	else{
		strcpy(strModel, pszMatchFile);
		strcpy(strTempDir, pszMatchFile);	char* pS = strrchr(strTempDir, '.');	if (!pS) pS = strTempDir + strlen(strTempDir);
		strcpy(pS, "_temp");	CreateDir(strTempDir);
	}
refine:
	fGsd_ratio = 1 / pCvt->get_gsd_ratio();
	if (fGsd_ratio < 0.9){
		CGeoImage img;
		img = imgL;		imgL = imgR;		imgR = img;
		CImagesCoorT* cvt = pCvt->CreateInvertT(false);
		delete[] pCvt;	pCvt = cvt;
		fGsd_ratio = 1 / fGsd_ratio;
	}

	if (!pszMatchFile){
		sprintf(strModel, "%s/%s=%s.model", CStdImage::get_workspace_path(), imgL.GetImageName(), imgR.GetImageName());
	}

	if (!pCvt->overlap_on_left(range, range, &stBlkColL, &stBlkRowL, &nBlkColsL, &nBlkRowsL)){
		LogPrint(0,
			"No overlap between the two images.");
		ret_code = 2;
		goto end2;
	}
	nBlkColsL = nBlkColsL - stBlkColL + 1;	nBlkRowsL = nBlkRowsL - stBlkRowL + 1;
	if (!pCvt->overlap_on_right(range, range, &stBlkColR, &stBlkRowR, &nBlkColsR, &nBlkRowsR)){
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

#define SAVE_MATCH_FILE(strModel,imgL,imgR)	\
	{	int nb[7] = { sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt), 0 };	\
	CMatchFile::WriteModelFile(strModel, \
	imgL.GetImageID(), imgL.GetImagePath(), \
	imgR.GetImageID(), imgR.GetImagePath(), \
	match_result->nameL, &match_result->xl, &match_result->yl, match_result->nameR, &match_result->xr, &match_result->yr, NULL, match_num, nb);	\
	LogPrint(0, "Save match result to %s.", strModel); }\

	if (range > rangeB){
		float zoomRate = rangeB / range;	float length = (float)sqrt(GlobalParam::g_match_buffer_size)*0.5;
		if (zoomRate*imgR.GetCols() < length) zoomRate = length / imgR.GetCols();	if (zoomRate*imgR.GetRows() < length) zoomRate = length / imgR.GetRows();
		if (zoomRate*imgL.GetCols() / fGsd_ratio < length) zoomRate = length*fGsd_ratio / imgL.GetCols();
		if (zoomRate*imgL.GetRows() / fGsd_ratio < length) zoomRate = length*fGsd_ratio / imgL.GetRows();
		if (zoomRate < 1) {
			match_result = SiftMatch_Pyrm(&match_num,
				&imgL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL,
				&imgR, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR,
				pCvt, range, zoomRate, siftDist, siftRatio, mutual_best_match, strTempDir);
			SAVE_MATCH_FILE(strModel, imgL, imgR);
			if (match_num > 1) {
				char strPath[512];	sprintf(strPath, "%s/%s=%s_range%d_zoom%d.model", strTempDir, imgL.GetImageName(), imgR.GetImageName(), (int)range, int(1 / zoomRate));
				CopyFile(strModel, strPath, FALSE);
				int nBufferSpace[4] = { sizeof(MPt), sizeof(MPt), sizeof(MPt), sizeof(MPt) };
				match_num = pCvt->update(&match_result->xl, &match_result->yl, &match_result->xr, &match_result->yr, nBufferSpace, true, NULL, match_num, NULL, (rangeG + rangeB) / 2);
				if (match_num > 0){
					FreeSiftMatchMem(match_result);
					range = rangeB;
					goto refine;
				}
			}
			FreeSiftMatchMem(match_result);
		}
		else range = rangeB;
	}
	if (range <= rangeB){
		if (!ExtractSiftFile(&imgL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, strTempDir, bUseExistSift, strSiftFileL)) goto end2;
		if (!ExtractSiftFile(&imgR, stBlkColR, stBlkRowR, nBlkColsR, nBlkRowsR, strTempDir, bUseExistSift, strSiftFileR)) goto end2;
		bUseExistSift = true;
		if (range > rangeG){
			match_result = SiftMatch_Block(&match_num,
				&imgL, strSiftFileL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL,
				&imgR, strSiftFileR, pCvt, range, siftDist, siftRatio, mutual_best_match);
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
							f_left_num = GeoTransform::EstimateFundamentalMatrix(&match_result->xl, &match_result->yl, &match_result->xr, &match_result->yr, nBufferSpace, match_num, NULL, F[0], argc, argv);
							if (!f_left_num) {
								LogPrint(0, "Fail to calculate the fundamental matrix between two image.");
								match_num = 0;
							}
						}
					}
				}
				else{
					match_num = pCvt->update(&match_result->xl, &match_result->yl, &match_result->xr, &match_result->yr, nBufferSpace, true, NULL, match_num, NULL, rangeG);
				}
				if (match_num > 0){
					range = rangeG;
					FreeSiftMatchMem(match_result);
					goto refine;
				}
			}
			FreeSiftMatchMem(match_result);
		}
		if (range <= rangeG){
			if (bUseH) hdistmax = range;	if (bUseF) fdistmax = range;
			if ((ratiomin_scale >= 0 && ratiomax_scale < 0) || (diffmin_ori >= 0 && diffmax_ori < 0)){
				float affineMatrix[6];
				if (pCvt->calculate_affine_matrix(stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL, stBlkColR, stBlkRowR, affineMatrix)){
					float tmp = 0;
					if (ratiomin_scale >= 0 && ratiomax_scale < 0){
						tmp = sqrt((affineMatrix[0] * affineMatrix[0] + affineMatrix[1] * affineMatrix[1] + affineMatrix[3] * affineMatrix[3] + affineMatrix[4] * affineMatrix[4]) / 2);
						ratiomax_scale = tmp + ratiomin_scale;	ratiomin_scale = tmp - ratiomin_scale;	if (ratiomin_scale < 0) ratiomin_scale = 0;
					}
					if (diffmin_ori >= 0 && diffmax_ori < 0){
						tmp = (float)fabs(atan2(affineMatrix[3], affineMatrix[0]));
						diffmax_ori = tmp + diffmin_ori;	if (diffmax_ori > PI) diffmax_ori = PI;
						diffmin_ori = tmp - diffmin_ori;	if (diffmin_ori < 0) diffmin_ori = 0;
					}
				}
			}
			match_result = SiftMatch_Geo(&match_num,
				&imgL, strSiftFileL, stBlkColL, stBlkRowL, nBlkColsL, nBlkRowsL,
				&imgR, strSiftFileR, pCvt, range, bUseH ? H : NULL, hdistmax, bUseF ? F : NULL, fdistmax,
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
