#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<time.h>
#define SN 4
#define L 25000000
#define CM 10000000
#define M 10000
#define MINI 0.000000001
#define CN 6
#define NBR 310
#define CNVG 4
#define E 2.71828
#define D 10000000.0
#define PT 0.9
#define TMDEP 60.0
#define PLOTP "/vol6/home/bgi_caofei/program/cnv/bin/ssCNVplot.R"
#define PERLP "/vol6/home/bgi_caofei/program/cnv/bin/ssGetGenes.pl"

float depth[L], rowmean, nmean, NB_p, vi[CM][CN][2], dis[CM], pp[CM], D_deno;
// depth, column mean, row mean
int d[CM], ncol, NB_r, copyN[CM], lim[M][5], cnvN, loc[CM][2];
int bestR;
double bestP;
//fixed dep(int), sample_number, colmumn number
const float ini[CN] = {0.01, 0.03, 0.9, 0.03, 0.02, 0.01};
const float endi[CN] = {0.01, 0.03, 0.9, 0.03, 0.02, 0.01};
//dangerous initialization

int batchFix(float dep[], int d[CM], int cn, char *oheader);
int nbinomFit(char *, int d[CM], int cn);
double nbinop(int k, int r, double p);
int fitp(double *, int, int, double *, double *);

int viterbiTraining(int d[CM], int copyN[CM], int	coln);
int acopy(int rS[], int eS[], int len);
int viterbi(float vi[][CN][2], int dep[], float tr[CN][CN], int cn[],	int len);
int mlEstimate(int dep[], int cn[], float tr[CN][CN], int len);
float emi(int copyn, int dep);
int splitCNV(int cnv[CM], float dep[], int len, char *);
int sortArr(float arr[], int arrLen);

int getVar(float d[][CM], float rowm[SN], double var[], int sn, int cn);
int getCor(float d[][CM], float rowm[SN], double var[], float cor[][SN], int sn, int cn);

double fw[CM][CN], bw[CM][CN], fs[CM], bs[CM];
float forwardP(double fw[][CN], float tran[CN][CN], int sdep[], int len);
float backwardP(double bw[][CN], float tran[CN][CN], int sdep[], int len);
float deviAtion(float dep[], int colc, int len, int sn, int rown);
float ddeviAtion(int d[][CM], int colc, int len, int sn, int rown);

int main(int argc, char *argv[])
{
	if(argc != 5) 
	{
		fprintf(stderr, "Usage:%s dep_file loc_file out_header chr\n", argv[0]);
		exit(1);
	}
	char outf[256], binhfn[256];
	sprintf(outf, "%s.cnv", argv[3]);
	FILE *in, *out, *binhf;
	if((binhf = fopen(argv[2], "r")) == NULL)
	{
		printf("Can't open file %s\n", argv[2]);
		exit(1);
	}
	if((in = fopen(argv[1], "r")) == NULL)
	{
		printf("Can't open %s\n", argv[1]);
		exit(1);
	}
	int i, j, binn = 0, binsn, isin = 0;
	char chr[8];
	float gcv;
	printf("running %s:\n", argv[4]);
	while(fscanf(binhf, "%s", chr) == 1)
	{
		fscanf(binhf, "%d", &binsn);
		fscanf(binhf, "%d", &loc[binn][0]);
		fscanf(binhf, "%d", &loc[binn][1]);
		fscanf(binhf, "%f", &gcv);
		fscanf(in, "%f", &depth[binn]);
		if(strcmp(chr,argv[4]) == 0)
		{
			binn++;
			if(isin == 0)
				isin == 1;
		}
		else if(isin == 1)
		{
			break;
		}
	}
	if(fclose(in) != 0 || fclose(binhf) != 0)
	{
		printf("Error closing file\n");
		exit(1);
	}
	ncol = binn;
  batchFix(depth, d, ncol, argv[3]);
	for(i = 1; i < binn; i++)
		dis[i] = (loc[i-1][1] - loc[i][0] + 0.0) / D;
	dis[0] = dis[1];
	nbinomFit(argv[3], d, ncol);
  viterbiTraining(d, copyN, ncol);
	splitCNV(copyN, depth, ncol, argv[3]);
	if(cnvN > 0)
	{
		if((out = fopen(outf, "w")) == NULL)
		{
			printf("Can't write file %s\n", outf);
			exit(1);
		}
		for(j = 0; j < ncol; j++)
		{
			if(copyN[j] > 0)
				nmean += depth[j] / (copyN[j]+0.0);
			else
				nmean += depth[j];
		}
		nmean /= (ncol+0.0);
		for(j = 0; j < ncol; j++)
			depth[j] /= nmean;
		for(j = 0; j < ncol; j++)
			fprintf(out, "%d ", j+1);
		fprintf(out, "\n");
		for(j = 0; j < ncol; j++)		
			fprintf(out, "%.2f ", depth[j]);
		fprintf(out, "\n");
		for(j = 0; j < ncol; j++)
			fprintf(out, "%d ", copyN[j]);
		fprintf(out, "\n");
		for(j = 0; j < ncol; j++)
			fprintf(out, "%.2f ", pp[j]);
		fprintf(out, "\n");
		for(j = 0; j < ncol; j++)
			fprintf(out, "%d ", loc[j][0]);
		fprintf(out, "\n");
		for(j = 0; j < ncol; j++)
			fprintf(out, "%d ", loc[j][1]);
		fprintf(out, "\n");
		if(fclose(out) != 0)
		{
			printf("Error closing out file\n");
			exit(1);
		}
//		char cmd[256];
//		sprintf(cmd, "R --slave --args %s < %s", outf, PLOTP);
//		system(cmd);
//		sprintf(cmd, "perl %s %s %s", PERLP, outf, argv[4]);
//		system(cmd);
	}
  return 0;
}

int batchFix(float dep[], int d[CM], int cn, char *oheader)
{
//	FILE *dmeanf;
//	char outfn[256];
//	sprintf(outfn, "%s.mean", oheader);
//  if((dmeanf = fopen(outfn, "w")) == NULL)
//	{
//		printf("Can't write file %s\n", outfn);
//		exit(1);
//	}
	if(cn > CM) 
	{
		printf("column length exceeds preset %d!\treset it\n",CM);
		exit(1);
	}
	int i, j;
	for(i = 0; i < cn; i++)
	{
		if(dep[i] < 0)
		{
			if(i > 0)
				dep[i] = dep[i-1];
			else
			{
				j = i + 1;
				while(dep[j] < 0)
					j++;
				dep[i] = dep[j];
			}
		}
		rowmean += dep[i];
	}
	rowmean /= (cn+0.0);
	D_deno = rowmean / TMDEP;
	for(i = 0; i < cn; i++)
		d[i] = (int)(0.5 + dep[i]/(D_deno+0.0));
	printf("================================\n");
	printf("got reads count.\ncolumn number: %d\n", cn);
	printf("mean depth: %.2f\n", rowmean);
	printf("dep denominator: %.1f\n", D_deno);
//	fprintf(dmeanf, "%.2f\n", rowmean);
//	if(fclose(dmeanf) != 0)
//	{
//		printf("Error closing file %s\n", outfn);
//		exit(1);
//	}
	return(0);
}

int nbinomFit(char *depf, int d[CM], int cn)
{
	printf("======================================\n");
	printf("Negative binomial fitting:\n");
	printf("NB_r value checked from 2 to %d\n", NBR);
//	FILE *densf;
//	char densfn[256];
//	sprintf(densfn, "%s.dens", depf);
//	if((densf = fopen(densfn, "w")) == NULL)
//	{
//		printf("Can't write density file\n");
//		exit(1);
//	}
	double dens[M], devi[11], p[11];
	int len = cn, max = 0, i, bp, r;
	for(i = 0; i < M; i++)
		dens[i] = 0;
	for(i = 0; i < len; i++)
	{
		if(d[i] > max)
			max = d[i];
		dens[d[i]]++;
	}
	int Max = max;
	if(max >= M)
	{
		printf("data maxium %d exceeds presetting %d\n", max, M);
		exit(1);
	}
	for(i = 0; i <= max; i++)
		dens[i] /= (len + 0.0);
	double mindevi = max;
	for(r = 2; r <= NBR; r++)
	{
		for(i = 0; i < 11; i++)
			p[i] = i * 0.1;
		for(i = 0; i < 4; i++)
			bp = fitp(dens, max, r, p, devi);
		if(devi[bp] < mindevi) 
		{
			mindevi = devi[bp];
			bestR = r;
			bestP = p[bp];
		}
	}
	NB_r = bestR;
	NB_p = bestP;
	mindevi /= (max + 1.0);
	printf("max depth of all: %d\n", max);
	printf("fixed depth fits NB(%d, %f) var: %e\n", bestR, bestP, mindevi);
//	fprintf(densf, "%d\t", bestR);
//	for(i = 0; i <= Max; i++)
//		fprintf(densf, "%d\t", i);
//	fprintf(densf, "\n%f\t", bestP);
//	for(i = 0; i <= Max; i++)
//		fprintf(densf, "%f\t", dens[i]);
//	fprintf(densf, "\n");
//	if(fclose(densf) != 0)
//	{
//		printf("Error closing density file\n");
//		exit(1);
//	}
	return 0;
}

int fitp(double dens[], int max, int r, double p[], double devi[])
{
	int i, j, bp, ans;
	double mindevi = max, gap = p[1] - p[0];
	for(i = 0; i < 11; i++)
	{
		devi[i] = 0;
		for(j = 0; j <= max; j++)
		{
			double deviation = dens[j] - nbinop(j,r,p[i]);
			devi[i] += deviation * deviation;
		}
		if(devi[i] < mindevi)
		{
			mindevi = devi[i];
			bp = i;
		}
	}
	if(bp == 0)
	{
		p[0] = p[0];
		ans = 0;
	}
	else
	{
		if(bp == 10)
		{
			p[0] = p[9];
			ans = 10;
		}
		else
		{
			if(devi[bp-1] < devi[bp+1])
			{
				ans = 10;
				p[0] = p[bp-1];
			}
			else
			{
				ans = 0;
				p[0] = p[bp];
			}
		}
	}
//	printf("\nbest p value: %f\nleast deviation: %f\n", p[bp], devi[bp]);
//	printf("gap\t%f\n", gap);
	gap /= 10;
	for(i = 1; i < 11; i++)
		p[i] = p[0] + gap * i;
	devi[ans] = mindevi;
	return(ans);
}

double nbinop(int k, int r, double p)
{
	double ans = p;
	int i;
	for(i = 0; i < r - 1; i++)
		ans *= p * (k + i + 1.0) / (i + 1.0);
	for(i = 0; i < k; i++)
		ans *= (1-p);
//	printf("nbinop(%d,%d,%f) = %f\n",k,r,p,ans);
	return(ans);
}


int viterbiTraining(int d[CM], int copyN[CM], int coln)
{
	printf("======================================\n");
	printf("Viterbi Training:\n");
  float tr[CN][CN];
		int *dep = d, cnEst[CM], cn[CM];
		int len = coln, i, j;
		float dmean = 0;
		for(i = 0; i < len; i++)
		{
			dmean += dep[i];
			cn[i] = -1;
		}
		dmean /= (len + 0.0);
		for(i = 0; i < len; i++)
		{
			cnEst[i] = 2 * dep[i] / dmean + 0.5;
			if(cnEst[i] >= CN)
				cnEst[i] = CN - 1;
		}
		int max = 2, min = 2;
		for(i = 0; i < len; i++)
		{
			if(cnEst[i] > max)
				max = cnEst[i];
			if(cnEst[i] < min)
				min = cnEst[i];
		}
		mlEstimate(dep, cnEst, tr, len);
		printf("===============================\n");
		printf("estimated transition matrix:\n");
		for(i = 0; i < CN; i++)
		{
			for(j = 0; j < CN; j++)
				printf("%.4f\t", tr[i][j]);
			printf("\n");
		}
		int indi = acopy(cnEst, cn, len);
		printf("length: %d\ndifferences: %d\t",len, indi);
		while(indi > 0)	
		{
			viterbi(vi, dep, tr, cnEst, len);
			mlEstimate(dep, cnEst, tr, len);
			indi = acopy(cnEst, cn, len);
			printf("%d\t", indi);
		}
		printf("\nestimated transition matrix:\n");
		for(i = 0; i < CN; i++)
		{
			for(j = 0; j < CN; j++)
				printf("%.4f\t", tr[i][j]);
			printf("\n");
		}
		float ratio[CN];
		for(i = 0; i < CN; i++)
			ratio[i] = 0;
		for(i = 0; i < len; i++)
			ratio[cn[i]]++;
		printf("copy number proportion:\n");
		for(i = 0; i < CN; i++)
		{
			ratio[i] /= (len + 0.0);
			printf("%.4f\t", ratio[i]);
		}
		printf("\n");
		float fwP = forwardP(fw, tr, dep, len);
    float bwP = backwardP(bw, tr, dep, len);
		float avP = (fwP + bwP) / 2;
		printf("forward possibility: %.6f\n", fwP);
		printf("backward possibility: %.6f\n", bwP);
		int posD = 0;
		for(i = 0; i < len; i++)
		{
			double lpp = 0, maxP = 0, lpp_j;
			int maxCN;
			for(j = 0; j <= i; j++)
				lpp += log10(fs[j]);
			for(j = i; j < len; j++)
				lpp += log10(bs[j]);
			for(j = 0; j < CN; j++)
			{
				lpp_j = lpp + log10(fw[i][j]) + log10(bw[i][j]) - avP;
				lpp_j = pow(10, lpp_j);
				if(lpp_j > maxP)
				{
					maxP = lpp_j;
					maxCN = j;
				}
			}
			pp[i] = maxP;
			copyN[i] = maxCN;
			if(maxCN != cn[i])
				posD++;
		}
		printf("dif between Viterbi and Posterior: %d\n", posD);
	return 0;
}


int mlEstimate(int dep[], int cnEst[], float tr[CN][CN], int len)
{
	int i, j;
  float sum[CN];
	for(i = 0; i < CN; i++)
	{
		sum[i] = 0;
		for(j = 0; j < CN; j++)
			tr[i][j] = 0;
	}
	for(i = 0; i < len - 1; i++)
		tr[cnEst[i]][cnEst[i+1]]++;
	for(i = 0; i < CN; i++)
		for(j = 0; j < CN; j++)
			sum[i] += tr[i][j];
	for(i = 0; i < CN; i++)
		if(sum[i] == 0)
			tr[i][i] = 1;
		else
			for(j = 0; j < CN; j++)
				tr[i][j] = tr[i][j] / sum[i];
	return 0;
}

int viterbi(float vi[][CN][2], int dep[], float tr[CN][CN], int cn[], int len)
{
	int i, j;
	for(i = 0; i < CN; i++)
		for(j = 0; j < CN; j++)
			tr[i][j] = log10(tr[i][j] + MINI);
	for(j = 0; j < CN; j++)
	{
		vi[0][j][0] =log10(emi(j, dep[0]) + MINI) + log10(ini[j]);
		vi[0][j][1] = -1;
	}
	for(i = 1; i < len; i++)
	{
		for(j = 0; j < CN; j++)
		{
			float vt[CN], maxvt = -100000000.0;
			int k, maxk = 0;
			for(k = 0; k < CN; k++)
			{
				float transp;
				if(k != j)
					transp = tr[k][j] * (1 - pow(E, dis[i]));
				else
				{
					transp = tr[k][j];
					int kk;
					for(kk = 0; kk < CN; kk++)
						if(kk != k)
							transp += tr[kk][j] * pow(E, dis[i]);
				}
				vt[k] = vi[i-1][k][0] + tr[k][j];
				if(vt[k] > maxvt)
				{
					maxvt = vt[k];
					maxk = k;
				}
			}
			vi[i][j][0] = maxvt + log10(emi(j, dep[i]) + MINI);
			vi[i][j][1] = maxk;
		}
	}
	int endC = 0;
	float endVt = -1000000000000.0;
  for(i = 0; i < CN; i++)
	{
		if(vi[len-1][i][0] > endVt)
		{
			endVt = vi[len-1][i][0];
			endC = i;
		}
	}
	cn[len-1] = endC;
	for(i = len - 1; i > 0; i--)
		cn[i-1] = vi[i][cn[i]][1];
	return 0;
}

int acopy(int arrA[], int arrB[], int len)
{
	int i, j = 0;
	for(i = 0; i < len; i++)
	{
		if(arrA[i] != arrB[i])
		{
			arrB[i] = arrA[i];
			j++;
		}
	}
	return(j);
}

float emi(int copyn, int dep)
{
	float nbinom_p = NB_p;
	double nb = nbinom_p;
	int i, dd, nbinom_r = NB_r;
	if(copyn == 0)
	{
		if(dep == 0)
			return(1.0 - MINI);
		else
			return(MINI);
	}
	else
		dd = 2 * dep / copyn + 0.5;
	for(i = 0; i < nbinom_r - 1; i++)
		nb *= nbinom_p * (dd + i + 1.0) / (i + 1.0);
	for(i = 0; i < dd; i++)
		nb *= (1 - nbinom_p);
	return(nb);
}

int splitCNV(int cnv[], float dep[], int len, char *oheader)
{
	FILE *locf;
	int i, j, *dept, *cnEst;
	float tr[CN][CN];
		dept = d;
		cnEst = copyN;
		mlEstimate(dept, cnEst, tr, len);
		char cnvdepn[256], locfn[256], cnvdn[256], cnvdmn[256];
		sprintf(cnvdepn, "%s.cnv.dep", oheader);
		sprintf(locfn, "%s.cnv.loc", oheader);
		sprintf(cnvdn, "%s.cnv.d", oheader);
		int copyn = -1, cnvn = 0, isin = 0, cnr[M];
		for(j = 0; j < len; j++)
		{
			if(cnv[j] == 2 && isin == 0)
				continue;
			if(isin == 1)
			{
				if(cnv[j] != copyn)
				{
					lim[cnvn][2] = j - 1;
					cnr[cnvn] = copyn;
					cnvn++;
					if(cnv[j] != 2)
					{
						copyn = cnv[j];
						lim[cnvn][1] = j;
						isin = 1;
					}
					else
						isin = 0;
				}
			}
			else
			{
				if(cnv[j] != 2)
				{
					lim[cnvn][1] = j;
					isin = 1;
					copyn = cnv[j];
				}
			}
		}
		if(isin = 1 && cnv[len-1] == copyn)
		{
			lim[cnvn][2] = len - 1;
			cnr[cnvn] = copyn;
			cnvn++;
		}
		cnvN = cnvn;
	if(cnvN > 0)
	{
		if((locf = fopen(locfn, "w")) == NULL)
		{
			printf("Can't open file %s\n", locfn);
			exit(1);
		}
		int k, m, n;
		for(k = 0; k < cnvn; k++)
		{
      float ppm = 0, ppl = 0, mdd = 0;
			int lim_begin = lim[k][1], lim_end = lim[k][2];
			int cnv_len = lim_end - lim_begin + 1;
			for(m = lim_begin; m <= lim_end; m++)
			{
				ppl += log10(emi(cnr[k], d[m])+MINI);
				ppl += log10(tr[cnr[k]][cnr[k]]+MINI);
				ppm += pp[m];
				mdd += d[m];
			}
			mdd /= (cnv_len+0.0);
			float ddr = mdd / rowmean;
			ppl += log10(pp[lim_begin]);
			ppl /= cnv_len;
			ppm /= cnv_len;
//			if(ppm < PT)
//				continue;
			if(lim_begin < CNVG || lim_begin < cnv_len)
				lim_begin = 0;
			else
				lim_begin -= (CNVG > cnv_len) ? CNVG : cnv_len;
			if(lim_end + CNVG >= len || lim_end + cnv_len >= len)
				lim_end = len - 1;
			else
				lim_end += (CNVG > cnv_len) ? CNVG : cnv_len;
//			if(ppm < 0.8)
//				continue;
//			if(ppl < -2.7 && ppl > -5)
//				continue;
//			if(dr > 0.6 && dr < 1.4 && ddr > 0.6 && ddr < 1.4) 
//				continue;
			fprintf(locf, "%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n",	lim_begin, lim[k][1], lim[k][2], lim_end, cnr[k], ppm, ppl, ddr);
		}
		if(fclose(locf) != 0) 
		{
			printf("Error closing file\n");
			exit(1);
		}
	}
	return(0);
}

int sortArr(float arr[], int arrLen)
{
	if(arrLen < 2)
		return(0);
	int i, j, indi;
	for(i = 0; i < arrLen - 1; i++)
	{
		indi = 0;
		for(j = arrLen - 2; j >= i; j--)
		{
			if(arr[j] > arr[j+1])
			{
				float temp = arr[j];
				arr[j] = arr[j+1];
				arr[j+1] = temp;;
				indi++;
			}
		}
		if(indi == 0)
			return(0);
	}
	return(0);
}

float forwardP(double fw[][CN], float tran[CN][CN], int sdep[], int len)
{
	int i, j;
	for(i = 0; i < len; i++)
		fs[i] = 0;
	for(i = 0; i < CN; i++)
		fw[0][i] = ini[i] * emi(i, sdep[0]);
	for(j = 0; j < CN; j++)
		fs[0] += fw[0][j];
	for(j = 0; j < CN; j++)
		fw[0][j] /= fs[0];
	for(i = 1; i < len; i++) 
	{
		for(j = 0; j < CN; j++)
		{
			int k;
			fw[i][j] = 0;
			for(k = 0; k < CN; k++)
				fw[i][j] += fw[i-1][k] * tran[k][j];
			fw[i][j] *= emi(j, sdep[i]);
			fs[i] += fw[i][j];
		}
		for(j = 0; j < CN; j++)
			fw[i][j] /= fs[i];
	}
	float ans = 0;
	for(j = 0; j < CN; j++)
		ans += fw[len-1][j] * endi[j];
	ans = log10(ans);
	for(i = 0; i < len; i++)
		ans += log10(fs[i]);
	return(ans);
}

float backwardP(double bw[][CN], float tran[CN][CN], int sdep[], int len)
{
	int i, j;
	for(i = 0; i < len; i++)
		bs[i] = 0;
	for(j = 0; j < CN; j++)
		bw[len-1][j] = endi[j];
	for(j = 0; j < CN; j++)
		bs[len-1] += bw[len-1][j];
	for(j = 0; j < CN; j++)
		bw[len-1][j] /= bs[len-1];
	for(i = len - 2; i >= 0; i--)
	{
		for(j = 0; j < CN; j++)
		{
			int k;
			bw[i][j] = 0;
			for(k = 0; k < CN; k++)
				bw[i][j] += tran[j][k] * emi(k, sdep[i+1]) * bw[i+1][k];
			bs[i] += bw[i][j];
		}
		for(j = 0; j < CN; j++)
			bw[i][j] /= bs[i];
	}
	float ans = 0;
	for(j = 0; j < CN; j++)
		ans += bw[0][j] * ini[j] * emi(j, sdep[0]);
	ans = log10(ans);
	for(i = 0; i < len; i++)
		ans += log10(bs[i]); 
	return(ans);
}

int getVar(float d[][CM], float rowm[SN], double var[], int sn, int cn)
{
	int i, j;
	for(i = 0; i < sn; i++)
	{
		var[i] = 0;
		for(j = 0; j < cn; j++)
			var[i] += (d[i][j] - rowm[i]) * (d[i][j] - rowm[i]);
	}
	return 0;
}

int getCor(float d[][CM], float rowm[SN], double var[], float cor[][SN], int sn, int cn)
{
	int i, j;
	for(i = 0; i < sn; i++)
	{
		cor[i][i] = 0;
		for(j = i+1; j < sn; j++)
		{
			int k;
			cor[i][j] = 0;
			for(k = 0; k < cn; k++)
				cor[i][j] += (d[i][k] - rowm[i]) * (d[j][k] - rowm[j]);
			cor[i][j] /= ( sqrt(var[i]) * sqrt(var[j]) );
			cor[j][i] = cor[i][j];
		}
	}
	return 0;
}

float deviAtion(float dep[], int colc, int len, int sn, int rown)
{
	float vari = 0;
	float odm = 0;
	float mdep[SN];
	int i, j;
	for(i = 0; i < sn; i++)
	{
		mdep[i] = 0;
		for(j = 0; j < len; j++)
			mdep[i] += dep[i*len+colc+j];
		mdep[i] /= len;
		odm += mdep[i];
	}
	odm /= sn;
	for(i = 0; i < sn; i++) 
		vari += pow((mdep[i] - odm),2);
	vari = sqrt(vari/(sn+0.0));
	float ans = (mdep[rown] - odm) / vari;
	ans = (ans > 0)?ans:(0-ans);
	return(ans);
}

float ddeviAtion(int d[][CM], int colc, int len, int sn, int rown)
{
	float vari = 0;
	float odm = 0;
	float mdep[SN];
	int i, j;
	for(i = 0; i < sn; i++)
	{
		mdep[i] = 0;
		for(j = 0; j < len; j++)
			mdep[i] += d[i][colc+j];
		mdep[i] /= len;
		odm += mdep[i];
	}
	odm /= sn;
	for(i = 0; i < sn; i++) 
		vari += pow((mdep[i] - odm),2);
	vari = sqrt(vari/(sn+0.0));
	float ans = (mdep[rown] - odm) / vari;
	ans = (ans > 0)?ans:(0-ans);
	return(ans);
}
