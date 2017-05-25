#include "features.h"
#include <string.h>
 
 
void getFeatures(double* signal1, double* signal2, double* Features, int nFFT, int winsize) {
        double Buff1[256], Buff2[256]; /*
        memcpy(Buff1, signal1, winsize);
        memcpy(Buff2, signal2, winsize);*/
        signalProcessing(signal1, Features, nFFT, winsize);
        signalProcessing(signal2, &Features[4], nFFT, winsize);
        Features[8] = (cov(signal1, signal2, winsize) + 1)*0.5;
}
 
void signalProcessing(double* signal, double* Features, int nFFT, int winsize) {
        double* spectrum = (double*)calloc(nFFT, sizeof(double));
        double* a = (double*)calloc(3, sizeof(double));
         
        // Hamming windowing
        int i;
        if(winsize == 256) {
                for(i = 0; i < winsize; i++); {
                        signal[i] = signal[i] * hamming[i];
                }
        }
         
        double maxFreq, maxVal;
        double varyans, rmsVal;
        int l = 3;
                 
        ar_parameter(signal, winsize, l, a);
        varyans = var(signal, winsize);
        rmsVal = rms(signal, winsize);
         
        dft(a,l, nFFT, spectrum);
         
        for(i = 0; i < nFFT; i++) {
                spectrum[i] = sqrt(pow((varyans/spectrum[i]),2));
        }
         
        maxFr(spectrum, (double*)&maxVal, (double*)&maxFreq, nFFT/2);
         
        Features[0] = (maxFreq + 1)/160;
        Features[1] = maxVal;
        Features[2] = rmsVal;
        Features[3] = varyans;
         
        free(spectrum);
        free(a);
}
// Forier Transform Method
void dft(double* signal, int datasize, int Nfft, double* spectrum){
        complex* X = (complex*)calloc(Nfft, sizeof(complex));
         
        int i, j;
        for(i = 0;i < Nfft; i++) {
                for(j = 0;j < datasize; j++) {
                        X[i].real += cos(2*pi*i*j/Nfft) * signal[j];
                        X[i].img  += -sin(2*pi*i*j/Nfft) * signal[j];
                }
                 
                spectrum[i] = sqrt(pow((X[i].real),2)+ pow((X[i].img),2));
        }
         
        free(X);
}
 // find AR(2) parameter
void ar_parameter(double* signal, int datasize, int l, double* a){
 
        double* rxx = (double*)calloc(l, sizeof(double));
                 
        int i = 0;
        int j = 0;
        for(j=0; j < l; j++) {
                for(i = j; i < datasize; i++) {
                        rxx[j] += signal[i-j] * signal[i];
                }
        }
         
        int k = 0;
        int n = 0;
        double detRxx = (pow(rxx[0],2)-pow(rxx[1],2));
        double Rxx[2][2];
        Rxx[0][0] = rxx[0]/detRxx;
        Rxx[0][1] = -rxx[1]/detRxx;
        Rxx[1][1] = rxx[0]/detRxx;
        Rxx[1][0] = -rxx[1]/detRxx;
         
        for(k = 0; k < l-1; k++){
                for(n = 0; n < l-1; n++){
                        a[0] = 1;
                        a[k+1] += -(Rxx[k][n] * rxx[n+1]);
                }
        }
             
        free(rxx);
}
 
 
 // Variance Method
double var(double* signal, int datasize){
        int i;
        double sum = 0;
        double sumpow = 0;
        double average = 0;
        double variance = 0;
        for (i = 0; i < datasize; i++)
        {
                sum += signal[i];
        }
        average = sum / datasize;
 
        for (i = 0; i < datasize; i++)
        {
                sumpow += pow((signal[i] - average), 2);
        }
        variance = sumpow /datasize;
         
        return variance;
}
 // find Maximum Frequence Method
void maxFr(double* array, double* maxValue, double* index, int size) {
        double maximum = array[0];
        int i, j = 1;
        for(i = 1; i < size; i++) {
                if(array[i] > maximum) {
                        maximum = array[i];
                        j = i + 1;
                }
        }
         
        *maxValue = maximum;
        *index = (double)j;
}
 // Root Mean Square Method
double rms(double* array, int size) {
        int i;
        double tmp = 0;
        for(i = 0; i < size; i++) {
                tmp += array[i] * array[i];
        }
         
        return sqrt(tmp/size);
}
 // Covariance Method
double cov(double* array1, double* array2, int size) {
        double* arrayX = (double*)calloc(size, sizeof(double));
        double* arrayY = (double*)calloc(size, sizeof(double));
        int i;
        double sum1 = 0, sum2 = 0, average1, average2, tmp = 0;
        for (i = 0; i < size; i++) {
                sum1 += array1[i];
                sum2 += array2[i];
        }
         
        average1 = sum1 / size;
        average2 = sum2 / size;
         
        for(i = 0; i < size; i++) {
                arrayX[i] = array1[i] - average1;
                arrayY[i] = array2[i] - average2;
                tmp +=  arrayX[i] * arrayY[i];
        }
        free(arrayX); free(arrayY);
        return tmp/size;
}