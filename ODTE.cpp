#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <gtest/gtest.h>
#include <tfhe/tlwe_functions.h>
#include <tfhe/tgsw_functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/polynomials_arithmetic.h>
#include <tfhe/lagrangehalfc_arithmetic.h>
#include "./libtfhe/tgsw-fft-operations.cpp"
#include "./libtfhe/tgsw-functions.cpp"
#include "./libtfhe/polynomials.cpp"
#include <tfhe/tfhe.h>

#define TFHE_TEST_ENVIRONMENT 1
using namespace std;

int main(){
    const int32_t w = 45;
    const int32_t depth = 20;
    double t_all = 0.;
    const int32_t N = 1024;
    const int32_t k = 1;
    const double min_stdev = pow(2.,-25);;
    const double max_stdev = 0.012467;
    const int32_t l = 4;
    const int32_t Bgbit = 8;

    srand(unsigned(time(NULL)));

    // int32_t minimum_lambda = 128;
    // TFheGateBootstrappingParameterSet *gb_params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    // const TGswParams *params = gb_params->tgsw_params;
    // const TLweParams *tlwe_params = params->tlwe_params;
    // cout<< tlwe_params->N <<endl;

    // const TGswParams *params = new_TGswParams(3, 10, new_TLweParams(1024, 1, 0., 1.));
    const TGswParams *params = new_TGswParams(l, Bgbit, new_TLweParams(N, k, min_stdev, max_stdev));
    const TLweParams *tlwe_params = params->tlwe_params;
    TGswKey *rgsw_key = new_TGswKey(params);
    TLweKey *rlwe_key = &rgsw_key->tlwe_key;
    TLweKey *rlwe_key_client = &(new_TGswKey(params))->tlwe_key;

    for (size_t depthi = 0; depthi < depth; depthi++)
    {
        // cout << "---" << depth << "---" << endl;
        TLweSample *lwe = new_TLweSample(tlwe_params);
        TLweSample *result = new_TLweSample(tlwe_params);
        TorusPolynomial *mu0 = new_TorusPolynomial(N);
        for (int32_t zeroi = 0; zeroi < N; ++zeroi) {
            mu0->coefsT[zeroi] = 0;
        }
        tLweNoiselessTrivial(result, mu0, tlwe_params);

        TGswSampleFFT *gswfft = new_TGswSampleFFT(params);
        TGswSample *gsw = new_TGswSample(params);

        int32_t temp_j  = rand() % w;
        int32_t temp_a = 0;
        for (size_t j = 0; j < w; j++)
        {
            // attribute value
            int32_t temp_ni  = rand() % N;
            TorusPolynomial *mu = new_TorusPolynomial(N);
            for (int32_t ni = 0; ni < temp_ni; ++ni) {
                mu->coefsT[ni] = modSwitchToTorus32(3, 4);
            }
            for (int32_t ni = temp_ni; ni < N; ++ni) {
                mu->coefsT[ni] = modSwitchToTorus32(1, 4);
            }
            tLweSymEncrypt(lwe, mu, min_stdev, rlwe_key);
            if (j == temp_j){
                tGswEncryptB(gsw, 1, min_stdev, rgsw_key);
                temp_a = temp_ni;
                // cout << "a:" << temp_a << endl;
            }
            else{
                tGswEncryptB(gsw, 0, min_stdev, rgsw_key);
            }
            tGswToFFTConvert(gswfft, gsw, params);

            time_t t_start0 = clock();    //------TIME BEGIN------
            tGswFFTExternMulToTLwe(lwe, gswfft, params);
            tLweAddTo(result, lwe, tlwe_params);
            time_t t_end0 = clock();    //------TIME END------
            t_all += (t_end0 - t_start0);

            delete_TorusPolynomial(mu);
        }

        // threshold value
        static const Torus32 mubar = modSwitchToTorus32(1, 4);
        TorusPolynomial *deacc = new_TorusPolynomial(N);
        TGswSampleFFT *gswtfft = new_TGswSampleFFT(params);
        TGswSample *gswt = new_TGswSample(params);
        IntPolynomial *mut = new_IntPolynomial(N);
        TorusPolynomial *mur = new_TorusPolynomial(N);
        TLweSample *lwer = new_TLweSample(tlwe_params);
        for (int32_t ni = 0; ni < N; ++ni) {
            mut->coefs[ni] = 0;
        }
        int32_t temp_t  = rand() % N;
        // cout << "t:" << temp_t << endl;
        mut->coefs[N-temp_t] = 1;
        tGswSymEncrypt(gswt, mut, min_stdev, rgsw_key);
        tGswToFFTConvert(gswtfft, gswt, params);
        time_t t_start1 = clock();    //------TIME BEGIN------
        tGswFFTExternMulToTLwe(result, gswtfft, params);
        (&result->a[k])->coefsT[0] += mubar;

        // add masking
        for (int32_t ri = 0; ri < N; ++ri) {
            mur->coefsT[ri] = modSwitchToTorus32(rand() % 2, 2);
        }
        tLweNoiselessTrivial(lwer, mur, tlwe_params);
        tLweAddTo(result, lwer, tlwe_params);

        tLweSymDecrypt(deacc, result, rlwe_key, 2);
        int32_t actual = modSwitchFromTorus32(deacc->coefsT[0]-mur->coefsT[0], 2);
        time_t t_end1 = clock();   //------TIME END------
        t_all += (t_end1 - t_start1);
        // cout << "comparison:" << actual << endl;
        if ((temp_a >= temp_t) != actual)
        {
            cout << depthi << " comparison failed!" << endl;
        }

        // leaf lable
        TorusPolynomial *muleaf = new_TorusPolynomial(N);
        TLweSample *lweleaf = new_TLweSample(tlwe_params);
        TLweSample *lweleaf_client = new_TLweSample(tlwe_params);
        int32_t vu = rand() % 2;
        muleaf->coefsT[0] = modSwitchToTorus32(vu, 2);
        for (int32_t leafi = 1; leafi < N; ++leafi) {
            muleaf->coefsT[leafi] = 0;
        }
        tLweSymEncrypt(lweleaf, muleaf, min_stdev, rlwe_key);
        // add masking
        TorusPolynomial *murl = new_TorusPolynomial(N);
        TLweSample *lwerl = new_TLweSample(tlwe_params);
        TorusPolynomial *deacc1 = new_TorusPolynomial(N);
        time_t t_start2 = clock();    //------TIME BEGIN------
        murl->coefsT[0] = modSwitchToTorus32(rand() % 2, 2);
        for (int32_t rli = 1; rli < N; ++rli) {
            murl->coefsT[rli] = 0;
        }
        tLweNoiselessTrivial(lwerl, murl, tlwe_params);
        tLweAddTo(lweleaf, lwerl, tlwe_params);
        // re-enc
        tLweSymDecrypt(deacc1, lweleaf, rlwe_key, 2);
        // int32_t actual1 = modSwitchFromTorus32(deacc1->coefsT[0], 2);
        tLweSymEncrypt(lweleaf_client, deacc1, min_stdev, rlwe_key_client);
        tLweAddTo(lweleaf_client, lwerl, tlwe_params);
        time_t t_end2 = clock();   //------TIME END------
        t_all += (t_end2 - t_start2);

        // client dec
        TorusPolynomial *deacc2 = new_TorusPolynomial(N);
        tLweSymDecrypt(deacc2, lweleaf_client, rlwe_key_client, 2);
        int32_t actual2 = modSwitchFromTorus32(deacc2->coefsT[0], 2);
        if (actual2 != vu)
        {
            cout << "client dec failed!" << endl;
        }

        delete_TGswSampleFFT(gswfft);
        delete_TGswSampleFFT(gswtfft);
        delete_TGswSample(gsw);
        delete_TGswSample(gswt);
        delete_TLweSample(lwe);
        delete_TLweSample(lwer);
        delete_TLweSample(lweleaf);
        delete_TLweSample(lweleaf_client);
        delete_TLweSample(result);
        delete_TorusPolynomial(deacc);
        delete_TorusPolynomial(deacc1);
        delete_TorusPolynomial(deacc2);
        delete_TorusPolynomial(muleaf);
        delete_TorusPolynomial(murl);
        delete_TorusPolynomial(mur);
        delete_TorusPolynomial(mu0);
        delete_IntPolynomial(mut);

    }
    cout << "(w,d): " << w << ", " << depth << endl;
    cout << "Time: " << t_all/CLOCKS_PER_SEC <<" secs"<< endl;

}