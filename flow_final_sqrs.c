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
int from, to, nsig, window, input_window, fv_size, nac, ndc, ac1, ac2, ac3, ac4, dc1, dc2, dc3, dc4;
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
 
        input_window = 256;

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

	full_sig = (double *)malloc((input_window)*sizeof(double));

	fp = fopen(argv[4],"r");//("data.txt","r");
	if (fp == 0) printf("Error: could not open data file.\n");
	i=0;
	for (j=0; j<1*input_window; j++) {
		fgets(line,80,fp);
		sscanf(line,"%lf", &rp);
		full_sig[i] = (rp - baseline)/gain;
		i = i + 1;
		if (i==input_window) { 
			flow(full_sig); //if(j>99*input_window)
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
	double *sqrs(),*rpeak,*beat,*coef,*noiseremoval(double * a, int b);
	freq = 360;
	gain = 200;
	baseline = 1024;
	nl = 0;
	from = 0;
	to = input_window;
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
	
	int len = input_window;

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
	rpeak = sqrs();
	gettimeofday(&t1, NULL);
	elapsed = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1000000);
	fprintf(fp,"%f\n",elapsed);
	fclose(fp); 
fp = fopen("res/sqrs_output.txt", "w+");
for (i=1; i<=rpeak[0]; i++) {
fprintf(fp,"%f\n",rpeak[i]);
}
fclose(fp); 

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


#define abs(A)	((A) >= 0 ? (A) : -(A))

int vector;
void my_setifreq();
int my_getvec();

double *sqrs()
{
    int filter, i, time = 0,
        slopecrit, sign, maxslope = 0, nslope = 0,
        qtime, maxtime, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9,
        ms160, ms200, s2, scmax, scmin = 500;
    long now, spm;
    double *tm;

    tm = (double *)malloc(40*sizeof(double));
   
    to = to*250/360;

    //setifreq(250.);
    my_setifreq();

    spm = 15000;
    //spm = strtim("1:0");

    scmin = 100;
    //scmin = muvadu((unsigned)signal, scmin);
    slopecrit = scmax = 10 * scmin;
    now = from;

    //ms160 = strtim("0.16"); ms200 = strtim("0.2"); s2 = strtim("2");
    ms160 = 40; ms200 = 50; s2 = 500;

    i = my_getvec(); //v[signal] = vector; printf("vector: %d\n",vector); 
    //(void)getvec(v);
    t9 = t8 = t7 = t6 = t5 = t4 = t3 = t2 = t1 = vector;//v[signal];
int count =1;
    int j = 0;
    do {//v[signal] = vector;
        filter = (t0 = vector) + 4*t1 + 6*t2 + 4*t3 + t4
                - t5         - 4*t6 - 6*t7 - 4*t8 - t9;
//if (count < 100) printf("t0: %d\n",t0); count++;
        if (time % s2 == 0) {
            if (nslope == 0) {
                slopecrit -= slopecrit >> 4;
                if (slopecrit < scmin) slopecrit = scmin;
            }
            else if (nslope >= 5) {
                slopecrit += slopecrit >> 4;
                if (slopecrit > scmax) slopecrit = scmax;
            }
        }
        if (nslope == 0 && abs(filter) > slopecrit) {
            nslope = 1; maxtime = ms160;
            sign = (filter > 0) ? 1 : -1;
            qtime = time;
        }
        if (nslope != 0) {
            if (filter * sign < -slopecrit) {
                sign = -sign;
                maxtime = (++nslope > 4) ? ms200 : ms160;
            }
            else if (filter * sign > slopecrit &&
                     abs(filter) > maxslope)
                maxslope = abs(filter);
            if (maxtime-- < 0) {
                if (2 <= nslope && nslope <= 4) {
                    slopecrit += ((maxslope>>2) - slopecrit) >> 3;
                    if (slopecrit < scmin) slopecrit = scmin;
                    else if (slopecrit > scmax) slopecrit = scmax;
//                    annot.time = now - (time - qtime) - 4;
		    tm[j+1] = round((double)(now - (time - qtime) - 4)*360/250);
//printf("annot.time: %f\n",tm[j+1]);
		    j++;
                    time = 0;
                }
                else if (nslope >= 5) {
//                    annot.time = now - (time - qtime) - 4;
		    tm[j+1] = now - (time - qtime) - 4;
		    j++;
                }
                nslope = 0;
            }
        }
        t9 = t8; t8 = t7; t7 = t6; t6 = t5; t5 = t4;
        t4 = t3; t3 = t2; t2 = t1; t1 = t0; time++; now++;
//i = my_getvec();
    } while (my_getvec() > 0 && (to == 0 || now <= to));
    //wfdbquit();
    tm[0] = (double)j;
   // for (i=0; i<=j; i++) printf("tm: %f\n",tm[i]);
//free(tm);
//free(samples);
    return(tm);
}

static long mticks, nticks, mnticks;
static long rgvtime, gvtime;
static int gv0, gv1, indx;

void my_setifreq()
{
    int sfreq = 360, ifreq = 250;
    int f = ifreq;
    int g = sfreq;
    double error; 
   
    while ((error = f - g) > 0.005 || error < -0.005)
	    if (f > g) f -= g;
	    else g -= f;
    mticks = (long)(sfreq/f + 0.5);
    nticks = (long)(ifreq/f + 0.5);
    mnticks = mticks * nticks;
    gvtime = 0;
    //rgvstat = rgetvec(gv0);
    //rgvstat = rgetvec(gv1);
    indx = from;
    gv0 = samp[indx]; indx++;
    gv1 = samp[indx]; indx++;
    rgvtime = nticks;
    return;
}

int my_getvec()
{
    if (indx>=input_window) return(-1);    
	 
    int i;

    if (rgvtime > mnticks) {
	rgvtime -= mnticks;
	gvtime  -= mnticks;
    }

    while (gvtime > rgvtime) {
	gv0 = gv1;
	//rgvstat = rgetvec(gv1);
        gv1 = samp[indx]; indx++;
	rgvtime += nticks;
    }
    vector = gv0 + (gvtime%nticks)*(gv1-gv0)/nticks;
    gv0 = gv1;
    gvtime += mticks;
    return(1);
}
