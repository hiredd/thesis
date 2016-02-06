#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "svm.h"
#include "model.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

double *samp;
void free();
int from, to, nsig, window, fv_size, nac, ndc, ac1, ac2, ac3, ac4, dc1, dc2, dc3, dc4;
char *filename1, *filename2;

//struct svm_model *model,*createmodel();

double *dbwavelet(double *conv);

int predict(double *coef, int lab, int fv_s);
void flow(double *sig);

int main(int argc, char *argv[]) {

	int gain, baseline, i, j;
	double *full_sig;

	gain = 200;
	baseline = 1024;

	filename1 = argv[2];
	filename2 = argv[3];

	double rp;
	char line[80];
	FILE *fp;

	fp = fopen(argv[1],"r");
	fgets(line,80,fp);
	if (strcmp(line,"coefficients\n")) printf("Error in svmmodel.txt file.\n");
	fgets(line,80,fp);
	ac1 = atoi(line);
	fgets(line,80,fp);
	ac2 = atoi(line);
	fgets(line,80,fp);
	ac3 = atoi(line);
	fgets(line,80,fp);
	ac4 = atoi(line);
	fgets(line,80,fp);
	dc1 = atoi(line);
	fgets(line,80,fp);
	dc2 = atoi(line);
	fgets(line,80,fp);
	dc3 = atoi(line);
	fgets(line,80,fp);
	dc4 = atoi(line);
	fclose(fp);

	model = svm_load_model(argv[5]);	
	//model = createmodel();

	nac = 0;
	ndc = 0;
	if (ac1==1) nac = nac + 130;
	if (ac2==1) nac = nac + 66;
	if (ac3==1) nac = nac + 34;
	if (ac4==1) nac = nac + 18;
	if (dc1==1) ndc = ndc + 130;
	if (dc2==1) ndc = ndc + 66;
	if (dc3==1) ndc = ndc + 34;
	if (dc4==1) ndc = ndc + 18;
	fv_size = nac + ndc;

	full_sig = (double *)malloc((3000)*sizeof(double));

	fp = fopen(argv[4],"r");//("data.txt","r");
	if (fp == 0) printf("Error: could not open data file.\n");
	i=0;
	for (j=0; j<2*3000; j++) {
		fgets(line,80,fp);
		sscanf(line,"%lf", &rp);
		full_sig[i] = (rp - baseline)/gain;
		i = i + 1;
		if (i==3000) { 
			flow(full_sig); //if(j>99*3000)
			i=0;
		}
	}

	fclose(fp);
	free(full_sig);
	//for(i=0; i<(model->l); i++){                                               
          //      free(model->SV[i]);   
           //     model->SV[i]=NULL;   
        //}
	svm_free_and_destroy_model(&model);
	exit(0);

}

void flow(double *sig) {

	int i,j,k,m,n,nl,freq,gain,baseline,label;
	double *wqrs(),*rpeak,*beat,*coef,*noiseremoval(double * a, int b);
	freq = 360;
	gain = 200;
	baseline = 1024;
	nl = 0;
	from = 0;
	to = 3000;
	nsig= 2;
    
	samp = (double *)malloc((to)*sizeof(double));
	
	struct timeval t0, t1;
	double elapsed;

	FILE *fp;
//	char *filename = malloc(13*sizeof(char));
//	strcpy(filename, "output1_");
//	strcat(filename, &cmb);
//	strcat(filename, ".txt");
	fp = fopen(filename1,"a");
	
	int len = 3000;

	/* NOISE REMOVAL */ 
	gettimeofday(&t0, NULL);
	sig = noiseremoval(sig,len);
	gettimeofday(&t1, NULL);
	elapsed = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1000000);
	fprintf(fp,"%f\t\t",elapsed);
//FILE *fp1,*fp2;
//fp1 = fopen("results/sig_filt_101_20_3000.txt", "w+");
//fp2 = fopen("results/samp_filt_101_20_3000.txt", "w+");
	for (i=0; i<len; i++) {
		samp[i] = sig[i]*gain + baseline;
//fprintf(fp1,"%f\n",sig[i]);
//fprintf(fp2,"%f\n",samp[i]);
	}
//fclose(fp1); 
//fclose(fp2);
	
	/* HEARTBEAT DETECTION */
	//rpeak = (double *)malloc(to*sizeof(double));
	gettimeofday(&t0, NULL);
	rpeak = wqrs();
	gettimeofday(&t1, NULL);
	elapsed = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1000000);
	fprintf(fp,"%f\n",elapsed);
	fclose(fp); 
//fp1 = fopen("results/wqrs_101_20_3000.txt", "w+");
//for (i=1; i<=rpeak[0]; i++) {
//fprintf(fp1,"%f\n",rpeak[i]);
//}
//fclose(fp1); 

   //     strcpy(filename, "output2_");
     //   strcat(filename, &cmb);
       // strcat(filename, ".txt");
	fp = fopen(filename2,"a");
//	free(filename);
//fp1 = fopen("results/DWT_SVM_101_20_3000.txt", "w+");
	int normal = 0;
	int abnormal = 0;
	int PR_window = 86;
	int QT_window = 170;
	window = 257;
  	beat = (double *)malloc((PR_window + QT_window+1)*sizeof(double));

	for (i=2; i<rpeak[0]; i++) {//rpeak[0]

		/* HEARTBEAT SEGMENTATION */
		for (j=PR_window; j>0; j--) {
			beat[PR_window - j] = sig[(int)rpeak[i] - j];
		}	
		for (j=0; j<=QT_window; j++) {
			beat[PR_window + j] = sig[(int)rpeak[i] + j];
		}			
//FILE * fp2; fp2 = fopen("results/hbeat_101_20_3000.txt","a"); for (j=0; j<window; j++) fprintf(fp2,"%f\t",beat[j]); fprintf(fp2,"\n"); fclose(fp2);
		/* FEATURE EXTRACTION */
		//coef = (double *)malloc(fv_size*sizeof(double));
		gettimeofday(&t0, NULL);
		coef = dbwavelet(beat);
		gettimeofday(&t1, NULL);
		elapsed = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1000000);
		fprintf(fp,"%f\t\t",elapsed);
//for (j=0; j<fv_size; j++) {
//fprintf(fp1,"%f\t",coef[j]);
//}
		/* CLASSIFICATION */
		gettimeofday(&t0, NULL);		
		label = predict(coef, -1, fv_size);
		gettimeofday(&t1, NULL);
		elapsed = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1000000);
		fprintf(fp,"%f\n",elapsed);

		if (label==1) normal = normal + 1; else abnormal = abnormal + 1;
		free(coef);
//fprintf(fp1,"\nlabel:\t%d\n",label);
	} 

	fclose(fp); 
//fclose(fp1);
	printf("rpeaks: %d\n",(int)rpeak[0]);
	printf("normal: %d\n\t\t\tabnormal: %d\n",normal,abnormal);

	free(samp);
	free(rpeak);
	free(beat);
	free(sig);

	return;
}



double *dbwavelet(double *conv) {
	int i,j,k,m,n,pad;
	double *F,*G,*cA1,*cD1,*coefs,*temp,*conv_d,*conv_sym;
	double Lo_D[4] = {-0.129409522550921, 0.224143868041857, 0.836516303737469, 0.482962913144690};
	double Hi_D[4] = {-0.482962913144690, 0.836516303737469, -0.224143868041857, -0.129409522550921};
	
	m = window + 20;
	n = 4;

	F = (double *)malloc((n+m-1)*sizeof(double));
	G = (double *)malloc((n+m-1)*sizeof(double));
	cA1 = (double *)malloc(((n+m-1)/2)*sizeof(double));
	cD1 = (double *)malloc(((n+m-1)/2)*sizeof(double));
	coefs = (double *)malloc(fv_size*sizeof(double));
	temp = (double *)malloc(ndc*sizeof(double));
	conv_d = (double *)malloc(m*sizeof(double));
	conv_sym = (double *)malloc((m)*sizeof(double));

	int iac = 0;
	int idc = 0;

	m = window;

   for (k=1; k<=4; k++) {

	
	for (i=0; i<10; i++) conv_sym[9-i] = conv[i];
	for (i=0; i<m; i++) conv_sym[10+i] = conv[i];
	for (i=0; i<10; i++) conv_sym[m+10+i] = conv[m-1-i];
	m = m + 20;	
	
	for (i = 1; i<=(n + m -1); i++) {
    		F[i-1] = 0;
    		for (j = MAX(1,i+1-n); j<=MIN(i,m); j++) {
        		F[i-1] = F[i-1] + conv_sym[j-1]*Lo_D[i-j+1-1]; 
    		}
	}

	for (i = 1; i<=(n + m -1); i++) {
    		G[i-1] = 0;
    		for (j = MAX(1,i+1-n); j<=MIN(i,m); j++) {
        		G[i-1] = G[i-1] + conv_sym[j-1]*Hi_D[i-j+1-1]; 
    		}
	}

	j = 0;
	for (i = 1; i<=(n + m -1); i++) {
    		if ((i%2) == 0) {
			cA1[j] = F[i-1];
			cD1[j] = G[i-1];
			j = j + 1;
    		}
	}

	pad = 5;
	for (i = 0; i<j-2*pad; i++) {
		conv[i] = cA1[i+pad];
		conv_d[i] = cD1[i+pad];
	}
	m = j-2*pad;
	
	//free(conv_sym);
		
	if ((k==1 && ac1==1) || (k==2 && ac2==1) || (k==3 && ac3==1) || (k==4 && ac4==1)) {
		for (i = 0; i<m; i++) {
		coefs[iac+i] = conv[i];
		}
		iac = iac + m;
	}

	if ((k==1 && dc1==1) || (k==2 && dc2==1) || (k==3 && dc3==1) || (k==4 && dc4==1)) {
		for (i = 0; i<m; i++) {
		temp[idc + i] = conv_d[i];
		}
		idc = idc + m;
	}

	if ((iac==nac) && (idc==ndc)) break;
   }

	if (ndc!=0) {
		for (i=0; i<ndc; i++) {
			coefs[iac + i] = temp[i];
		}
	}	

	free(F);
	free(G);
	free(cD1);
	free(cA1);
	free(conv_sym);
	free(conv_d);
	free(temp);

	return(coefs);
}


#define BUFLN  16384	/* must be a power of 2, see ltsamp() */
#define EYE_CLS 0.25    /* eye-closing period is set to 0.25 sec (250 ms) */ 
#define MaxQRSw 0.13    /* maximum QRS width (130ms) */                        
#define NDP	 2.5    /* adjust threshold if no QRS found in NDP seconds */
#define PWFreqDEF 60    /* power line (mains) frequency, in Hz (default) */
#define TmDEF	 100	/* minimum threshold value (default) */

double lfsc;		/* length function scale constant */
int *ebuf;
int LPn, LP2n;          /* filter parameters (dependent on sampling rate) */
int LTwindow;           /* LT window size */
int PWFreq;		/* power line (mains) frequency, in Hz */
int Tm;			/* minimum threshold value */
double *lbuf;

double ltsamp(int t)
{
    int dy,sig=0;
    static int Yn, Yn1, Yn2;
    static int tt;
    static int aet, et;
    double v0, v1, v2;

    if (lbuf == NULL) {

	aet = 0;
	tt = (int)-1L;

	lbuf = (double *)malloc((unsigned)BUFLN*sizeof(double));
	ebuf = (int *)malloc((unsigned)BUFLN * sizeof(int));
	
	    for (ebuf[0] = sqrt(lfsc), tt = 1L; tt < BUFLN; tt++)
		ebuf[tt] = ebuf[0];
	    if (t > BUFLN) tt = (int)(t - BUFLN);
	    else tt = (int)-1L;
	    Yn = Yn1 = Yn2 = 0;
    }
    
    while (t > tt) {
	Yn2 = Yn1;
	Yn1 = Yn;
	 if ( (tt<to) && (tt-LPn<to) && (tt-LP2n<to) && (tt>=from) && (tt-LPn>=from) && (tt-LP2n>=from) ) {
	   v0 = samp[tt];
	   v1 = samp[tt-LPn];
	   v2 = samp[tt-LP2n];
	   Yn = 2*Yn1 - Yn2 + v0 - 2*v1 + v2;
	}
	dy = (Yn - Yn1) / LP2n;		/* lowpass derivative of input */
	et = ebuf[(++tt)&(BUFLN-1)] = sqrt(lfsc +dy*dy); /* length transform */
	lbuf[(tt)&(BUFLN-1)] = aet += et - ebuf[(tt-LTwindow)&(BUFLN-1)];
	/* lbuf contains the average of the length-transformed samples over
	   the interval from tt-LTwindow+1 to tt */
    }
    return (lbuf[t&(BUFLN-1)]);
}

double *wqrs(){
 
    lbuf = NULL;
    Tm = TmDEF;
    PWFreq = PWFreqDEF;
    char *p;
    float sps;			     /* sampling frequency, in Hz (SR) */
    float samplingInterval;          /* sampling interval, in milliseconds */
    int i, max, min, minutes = 0, onset, timer, vflag = 0;
    int dflag = 0;		     /* if non-zero, dump raw and filtered
					samples only;  do not run detector */
    int jflag = 0;		     /* if non-zero, annotate J-points */
    int Rflag = 0;		     /* if non-zero, resample at 120 or 150 Hz
				      */
    int EyeClosing;                  /* eye-closing period, related to SR */
    int ExpectPeriod;                /* if no QRS is detected over this period,
					the threshold is automatically reduced
					to a minimum value;  the threshold is
					restored upon a detection */
    double Ta, T0, *tm;		     /* high and low detection thresholds */

    int gain;
    int next_minute, spm, t, tj, tpq, tt, t1;
    
    tm = (double *)malloc(4000*sizeof(double));

    int sig = 0;

    sps = 360.000000;
    gain = 200;

    Tm = 20;	//Tm = muvadu((unsigned)sig, Tm);
    samplingInterval = 1000.0/sps;
    lfsc = 1.25*gain*gain/sps;	/* length function scale constant */
    spm = 60 * sps;
    next_minute = from + spm;
    LPn = sps/PWFreq; 		/* The LP filter will have a notch at the power line (mains) frequency */
    if (LPn > 8)  LPn = 8;	/* avoid filtering too agressively */
    LP2n = 2 * LPn;
    EyeClosing = sps * EYE_CLS; /* set eye-closing period */
    ExpectPeriod = sps * NDP;	/* maximum expected RR interval */
    LTwindow = sps * MaxQRSw;   /* length transform window size */

    if ((t1 = 2880) > BUFLN*0.9)
	t1 = BUFLN/2;
    t1 += from;
    //for (T0 = 0, t = from; t < t1 && sample_valid(); t++)
    for (T0 = 0, t = from; t < t1; t++)
	T0 += ltsamp(t);
    T0 /= t1 - from;
    Ta = 3 * T0;

    /* Main loop */
    int j = 0;
    //for (t = from; t < to || (to == 0L && sample_valid()); t++) {
    for (t = from; t < to; t++) {
	static int learning = 1, T1;
	
	if (learning) {
	    if (t > t1) {
		learning = 0;
		T1 = T0;
		t = from;	/* start over */
	    }
	    else
		T1 = 2*T0;
	}
	
	
	/* Compare a length-transformed sample against T1. */
	if (ltsamp(t) > T1) {	/* found a possible QRS near t */
	    timer = 0; /* used for counting the time after previous QRS */
	    max = min = ltsamp(t);
	    for (tt = t+1; tt < t + EyeClosing/2; tt++)
		if (ltsamp(tt) > max) max = ltsamp(tt);
	    for (tt = t-1; tt > t - EyeClosing/2; tt--)
		if (ltsamp(tt) < min) min = ltsamp(tt);
	    if (max > min+10) { /* There is a QRS near tt */
		/* Find the QRS onset (PQ junction) */
		onset = max/100 + 2;
		tpq = t - 5;
		for (tt = t; tt > t - EyeClosing/2; tt--) {
		    if (ltsamp(tt)   - ltsamp(tt-1) < onset &&
			ltsamp(tt-1) - ltsamp(tt-2) < onset &&
			ltsamp(tt-2) - ltsamp(tt-3) < onset &&
			ltsamp(tt-3) - ltsamp(tt-4) < onset) {
			tpq = tt - LP2n;	/* account for phase shift */
			break;
		    }
		}

		if (!learning) {
		    /* Check that we haven't reached the end of the record. */
		    //(void)sample(sig, tpq);
		    //if (sample_valid() == 0) break;
		    if (tpq>=to) break;
		    tm[j+1] = (double)tpq; 
		    j = j + 1; 
		    if (jflag) {
			/* Find the end of the QRS */
			for (tt = t, tj = t + 5; tt < t + EyeClosing/2; tt++) {
			    if (ltsamp(tt) > max - (max/10)) {
				tj = tt;
				break;
			    }
			}
			//(void)sample(sig, tj);
			//if (sample_valid() == 0) break;
			if (tj>=to) break;
		    }
		}

		/* Adjust thresholds */
		Ta += (max - Ta)/10;
		T1 = Ta / 3;


		/* Lock out further detections during the eye-closing period */
		t += EyeClosing;
	    }
	}
	else if (!learning) {
	    /* Once past the learning period, decrease threshold if no QRS
	       was detected recently. */
	    if (++timer > ExpectPeriod && Ta > Tm) {
		Ta--;
		T1 = Ta / 3;
	    }      
	}

	/* Keep track of progress by printing a dot for each minute analyzed */
	if (t >= next_minute) {
	    next_minute += spm;
	    //(void)fprintf(stderr, ".");
	    //(void)fflush(stderr);
	    if (++minutes >= 60) {
		//(void)fprintf(stderr, " %s\n", timstr(t));
		minutes = 0;
	    }
	}
    }
    //if (minutes) (void)fprintf(stderr, " %s\n", timstr(t));

    (void)free(lbuf);
    (void)free(ebuf);
    tm[0] = (double)j;
  	
    return(tm);
}





