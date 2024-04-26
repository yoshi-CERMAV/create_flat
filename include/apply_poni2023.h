//
//  apply_poni2003.h
//  
//
//  Created by Yoshiharu Nishiyama on 14/12/2023.
//

#ifndef apply_poni2003_h
#define apply_poni2003_h
#include <detector.h>
#include <map>
#include <rotator.h>
#include <Accelerate/Accelerate.h>
#include <vector>

typedef struct {
    float q;
    float beta;
    float intensity;
} plot_data_3d;

typedef struct {
    float q;
    float avg;
    float stdev;
    int count;
} plot_data;

typedef struct {
    float avg;
    float stdev;
    int count;
} average_data;

const double four_pi = M_PI * 4;
class Comparator
{
public:
    Comparator(float *d){dat = d;}
    bool operator() (size_t i, size_t j){return dat[i*2] < dat[j*2];}
private:
    float *dat;
};
struct
{
    bool operator()(double &a, double &b) const {if(isnan(a)) return true; if(isnan(b)) return false; return a < b; }
} customLess;

class apply_poni
{
public:
    apply_poni(detector *in_detector){
        Detector = in_detector;
        allocate(Detector->image_size);
        mask =new char[Detector->image_size];
        for(int i = 0; i <Detector->image_size; i++ ) mask[i] = 0;
    }
    apply_poni(const char poni_filename[], detector *in_detector){
        Detector = in_detector;
        allocate(Detector->image_size);
        read_poni(poni_filename);
        dim0=Detector->dim0;
        init_map();
        if(Detector->flat){
            response = Detector->response; // just copying the pointer
            for(int i = 0; i < Detector->image_size; i++ ){
                correct_response[i] = response[i] * center_corr[i];
            }
        }else{
            for(int i = 0; i < Detector->image_size; i++ ){
                correct_response[i] =  center_corr[i];
                }
        }
    }
    void print_pos(int pos){
        cout << center_qb[pos*2]<<endl;
    }
    void check_flat(){
        if(Detector->flat){
            response = Detector->response; // just copying the pointer
            for(int i = 0; i < Detector->image_size; i++ ){
                correct_response[i] = response[i] * center_corr[i];
            }
        }
        else {
            for(int i = 0; i < Detector->image_size; i++ ){
                correct_response[i] =  center_corr[i];
            }
            response = NULL;
        }
    }
    
    void read_poni(const char poni_filename[])
    {
        char key[255];
        double value;
        char *buffer = new char [1024];
        
        ifstream ponifile(poni_filename);
        
        while(ponifile.getline(buffer, 1024)){
            sscanf(buffer,"%[^:]:%lf",key, &value);
            poni[key]=value;
            cout << key <<": "<< value <<endl;
        }
        cout <<": "<<  poni["Poni1"]<<endl;
        r.init();
        r.roty(-poni["Rot1"]);
        r.rotx(-poni["Rot2"]);
        r.swap(); // swapping x and z;
        center[0] = poni["Distance"];
        center[1] = -poni["Poni1"];
        center[2] = -poni["Poni2"];
        cout << "center "<<center[0] <<" "<<center[1]<<" "<<center[2]<<endl;
        wavelength = poni["Wavelength"];
        init_map();
    }
    
    void allocate(int image_size){
        qb = new float[8*image_size];
        pixel_index = new size_t[image_size];
        center_qb = new float[image_size * 2];
        center_corr = new float[image_size];
        correct_response = new float[image_size];
        work = new double [image_size];
    }
    
    double get_q(int x, int y){
        size_t pos = y*dim0 + x;
        return center_qb[pos*2];
    }
    
    void init_map()
    {
        float *ptr = Detector->pixel_corners;
        float *ptr_qb= qb;
        scale = four_pi / wavelength * 1e-10;
        int num_point = Detector->image_size * 4;
        for(int i = 0; i < num_point ; i++, ptr+=3, ptr_qb+=2){
            add(ptr, center);
            r.apply(ptr);
            double r2 = sqrt(ptr[0] * ptr[0] + ptr[1] * ptr[1]);
            double theta = 0.5 * atan2(r2, ptr[2]);
            ptr_qb[0] = scale * sin(theta);
            ptr_qb[1] = atan2(ptr[0], ptr[1]);
        }
        cout << "checking if on the same side"<<endl;
        float M_PI2 = 2*M_PI;
        for(int i = 0; i < Detector->image_size; i++){
            float *ptr = qb+i*8+1;
            for(int j = 2; j < 8; j+=2){
                if(fabs(ptr[0]-ptr[j])>M_PI){
                    ptr[j] = (ptr[j] > 0 ) ? ptr[j]-M_PI2 :ptr[j]+M_PI2;
                }
            }
        }
        cout << "done"<<endl;
        // calculate center of pixels
        cout << "calculating correction"<<endl;
        //calculate correction factors
        ptr = Detector->center_pos;
        ptr_qb = center_qb;
        float norm[3];
        norm[0] = 1; norm[1] = 0; norm[2] = 0;
        r.apply(norm);
        float temp[3];
        cout << "center "<<center[0] <<" "<<center[1] <<" "<<center[2] << endl;
        for(int i = 0; i < Detector->image_size; i++, ptr+=3, ptr_qb+=2){
            for(int i = 0; i < 3; i++) temp[i] = ptr[i];
            add(temp, center);
            r.apply(temp);
            double r2 = sqrt(temp[0] * temp[0] + temp[1] * temp[1]);
            double Theta = atan2(r2, temp[2]);
            ptr_qb[0] = scale * sin(0.5*Theta);
            double beta = atan2(temp[0], temp[1]);
            ptr_qb[1] =beta;
            
            double cT = cos(Theta);
            double sT = sin(Theta);
            cT*=cT; sT*=sT;
            float corr = 0.5*(1 + cT + 1*cos(beta * 2) * sT); // phase 90˚ diff from -cos..
            corr = 1./corr;
            double r3 =temp[0] * temp[0] + temp[1] * temp[1] + temp[2]*temp[2];
            corr *= fabs(r3*sqrt(r3)/scalar(temp, norm));
            center_corr[i] = corr;
        }
        cout << center_qb[0]<<" "<<center_qb[1]<<endl;
        cout <<"going to sort "<<endl;
        // sort as function of q.
        if(Detector->masked){
            char *mask = Detector->mask;
            valid_count = 0;
            for(int i = 0; i!= Detector->image_size; i++) {
                if(mask[i]) continue;
                else{
                    pixel_index[valid_count] = i;
                    valid_count++;
                }
            }
        }
        else{
            for(int i = 0; i!= Detector->image_size; i++) pixel_index[i]=i;
            valid_count = Detector->image_size;
        }

        cout << "valid pixels" <<valid_count <<endl;
        Comparator comp(center_qb);
        sort(pixel_index, pixel_index + valid_count, comp);
        cout <<"sorted" <<endl;
        auto [q, t] = std::div((int)pixel_index[0], 578);
        cout << "q t "<< q<<" "<<t<<" "<<center_qb[pixel_index[0]*2]<<endl;
        auto [q1, t1] = std::div((int)pixel_index[valid_count-1], 578);
        cout << "q t "<< q1<<" "<<t1<<" "<<center_qb[pixel_index[valid_count-1]*2]<<endl;
    }
    
    double qmin(){return center_qb[pixel_index[0]*2];}
    double qmax(){return center_qb[pixel_index[valid_count-1]*2];}
    void add(float *a, float *b)
    {
        a[0] += b[0];
        a[1] += b[1];
        a[2] += b[2];
    }
    float sq(float *x)
    {
        return x[0] * x[0] + x[1]*x[1] + x[2]*x[2];
    }
    float scalar(float *a, float *b)
    {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    void write(ofstream &fo, double data){
        fo.write(reinterpret_cast<char *>(&data), sizeof(double));
    }
    void write(ofstream &fo, int data){
        fo.write(reinterpret_cast<char *>(&data), sizeof(int));
    }
    
    void select_qrange(double qmin, double qmax, int *dat, float *out_data, int &count)
    {
        check_flat();
        if(Detector->flat) select_qrange(qmin, qmax, dat, out_data, count, &apply_poni::get_flattened_value);
        else select_qrange(qmin, qmax, dat, out_data, count, &apply_poni::get_value);
    }
    typedef  float (apply_poni::*FreeFn)(int p);
    float get_value(int pos){return data[pos];}
    float get_flattened_value(int pos){return data[pos] * response[pos];}

    void select_qrange(double qmin, double qmax, int *dat, float *out_data, int &count, FreeFn func)
    {
        int i = 0;
        float *ptr = out_data;
        int count1 = 0;
        double q;
        while(center_qb[pixel_index[i]*2] < qmin) {
            i++;
            //            if (i> valid_count)  need to implement exception throw
        }
        q = center_qb[pixel_index[i]*2];
        while(q < qmax && i < valid_count){
            int pos =pixel_index[i];
            int pos2 = pos*2;
            ptr[0] = center_qb[pos2];
            ptr[1] = center_qb[pos2+1];
            ptr[2] = (this->*func)(pos);
            count1++;
            //need to check if count1 < count;
        }
        count = count1;
    }
    
    void add_mask(const char *mask1){
        for(int i = 0; i != Detector->image_size; i++) if(mask1[i]) mask[i] = 1;
    }
    
    void mask_pos(int i){mask[i] = 1;}
    
    int mask_line(int n){
        char *ptr = mask + n*dim0;
        for(int i = 0; i != dim0; i++ ) ptr[i] = 1;
        return 0;
    }
    
    void average(size_t *pix, int len,  int *dat)
    {
        reset();
        for(int i = 0; i < len; i++){
            size_t p = pix[i];
            double value = dat[p];
            value *= correct_response[p];
            if(value > 1e6 || value < 1) continue;
            if(isnan(value)) continue;

            add_in(value);
        }
        finish_average();
    }
    
    template<class T>
    void average_sort(size_t *pix, int len,  T *dat, float cutoff)
    {
        reset();
  //      cout << len<<endl;
        int icount = 0;
        for(int i = 0; i < len; i++){
            size_t p = pix[i];
            double value = dat[p];
            value *= correct_response[p];
            work[icount] = value;
                if(value > 1e19 || value < 1) continue;
                if(isnan(value)) continue;
            icount++;
        }
        sort(work, work+icount, customLess);
        int len4 = icount*cutoff;
        int i1 = icount-len4;
        for(int i = len4; i != i1; i++ ){
            add_in(work[i]);
        }
        finish_average();
    }

    template<class T>
    void robust_average(size_t *pix, int len, T *dat, double vmin, double vmax)
    {
        reset();
        float *response = Detector->response;
        for(int i = 0; i < len; i++){
            size_t p = pix[i];
            double value =  dat[p] ;
            value *= correct_response[p];
            if(isnan(value)){mask_pos(p); continue;}

            if(value < vmin){
                mask_pos(p);
                continue;
            }
            if(value > vmax){
                mask_pos(p);
                continue;
            }
            add_in(value);

        }
        finish_average();
    }

    void average(size_t *pix, int len, int *bg, int *dat, double scale)
    {
        reset();
        for(int i = 0; i < len; i++){
            size_t p = pix[i];
            double value = scale * dat[p] - bg[p];
            value *= correct_response[p];
            add_in(value);
        }
        finish_average();
    }
    
    void robust_average(size_t *pix, int len, int *bg, int *dat, double scale, double vmin, double vmax)
    {
        reset();
        float *response = Detector->response;
        for(int i = 0; i < len; i++){
            size_t p = pix[i];
            double value = scale * dat[p] -  bg[p];
            value *= correct_response[p];
            if(isnan(value)){mask_pos(p); continue;}
            if(value < vmin){
                mask_pos(p);
                continue;
            }
            if(value > vmax){
                mask_pos(p);
                continue;
            }
            add_in(value);
        }
        finish_average();
    }
    
    template<class T>
    void integrate(T *dat, average_data *avg)
    {
        cout <<"integrating "<<kmax <<endl;
 //       data = dat;
        int i0 = avg_begin;
        avg_ptr = avg;
        for(int k = 0; k < kmax; k++){
            size_t *pix = pixel_index + i0;
            int len = avg_end[k] - i0;
            i0 = avg_end[k];
            average_sort(pix, len, dat, 0.25);
            double vmin = sum-2*stdev;
            double vmax = sum+2*stdev;

            robust_average(pix, len, dat, vmin, vmax);
 //           cout <<"k "<<k<<" "<< current_q <<" "<<avg_ptr->avg<<endl;

            copy_avg();
            increment_q();
   //         cout <<"k "<<k<<" "<< current_q <<" "<<sum<<endl;

        }
        //     fo.close();
    }
    template <class T>
    void make_flat(float *flat, float q0, float q1,  float qstep, average_data *avg, T *dat)
    {
        for(int i = 0; i < Detector->image_size; i++) flat[i] = 1.;
        int i0 = avg_begin;
        for(int k = 0; k < kmax; k++){
            for(int i = i0; i < avg_end[k]; i++){
                int p=pixel_index[i];
                float qa = (center_qb[p*2]-q0+qstep*0.5)/qstep;
                int qi = (int)(qa);
                float qf = qa-qi;
                float val(0);
                if(qf < 0 && k) val = avg[k].avg * (1+qf) - avg[k-1].avg * qf ;
                if(qf > 0 && k != (kmax-1)) val = avg[k].avg * (1-qf) + avg[k+1].avg * qf ;
                if(val&& dat[p]){
                    flat[p] = val/dat[p]/correct_response[p];
                }
            }
            i0 = avg_end[k];
        }
        ofstream fo("flat.dat");
        fo.write(reinterpret_cast<char *>(flat), sizeof(float)*Detector->image_size);
    }
    
    template <class T>
    void make_flat(float *flat, float *flat_wgt, float q0, float q1,  float qstep, average_data *avg, T *dat)
    {
        for(int i = 0; i < Detector->image_size; i++){
            flat[i] = 1.;
            flat_wgt[i] = 0;
        }
        int i0 = avg_begin;
        for(int k = 0; k < kmax; k++){
            for(int i = i0; i < avg_end[k]; i++){
                int p=pixel_index[i];
                float qa = (center_qb[p*2]-q0+qstep*0.5)/qstep;
                int qi = (int)(qa);
                float qf = qa-qi;
                float val(0);
                if(qf < 0 && k) val = avg[k].avg * (1+qf) - avg[k-1].avg * qf ;
                if(qf > 0 && k != (kmax-1)) val = avg[k].avg * (1-qf) + avg[k+1].avg * qf ;
                if(val&& dat[p]){
                    flat[p] = val/dat[p]/correct_response[p];
                    flat_wgt[i] = sqrt(dat[p]);
                }
            }
            i0 = avg_end[k];
        }
        ofstream fo("flat.dat");
        fo.write(reinterpret_cast<char *>(flat), sizeof(float)*Detector->image_size);
        fo.write(reinterpret_cast<char *>(flat_wgt), sizeof(float)*Detector->image_size);
    }
    

    void integrate(int *bg_in, int *dat, double scale, average_data *avg)
    {
        bg = bg_in;
        data = dat;

        int i0 = avg_begin;
        avg_ptr = avg;
        for(int k = 0; k < kmax; k++){
            size_t *pix = pixel_index + i0;
            int len = avg_end[k] - i0;
            i0 = avg_end[k];
            average(pix, len, bg, dat, scale);
            //        cout << current_q <<" "<<avg<<endl;
            double vmin = sum-2*stdev;
            double vmax = sum+2*stdev;
            robust_average(pix, len, bg, dat, scale, vmin, vmax);
            increment_q();
            copy_avg();
        }
        //     fo.close();
    }

    void integrate(double qmin, double qmax, double qstep,
                   int *bg_in, int *dat, double scale, plot_data *plt, int &size)
    {
        bg = bg_in;
        data = dat;
        check_flat();
        cout << "integrating"<<endl;
        cout << scale <<endl;
        init_integrate(qmin, qmax, qstep, plt);
        int i0 = avg_begin;
        current_q = qmin;
        for(int k = 0; k < kmax; k++){
            size_t *pix = pixel_index + i0;
            int len = avg_end[k] - i0;
            i0 = avg_end[k];
            average(pix, len, bg, dat, scale);
            //        cout << current_q <<" "<<avg<<endl;
            double vmin = sum-3*stdev;
            double vmax = sum+3*stdev;
            robust_average(pix, len, bg, dat, scale, vmin, vmax);
            increment_q();
            copy_plt();
            //          cout << current_q <<" "<< avg<<endl;
            //         write_in();
        }
        size = plt_ptr-plt;
        //     fo.close();
    }
    void integrate(double qmin, double qmax, double qstep,
                   int *bg_in, int *dat, double *flat, double scale, char filename[])
    {
        bg = bg_in;
        data = dat;
        
        double inv_scale = 1./scale;
        ofstream fo(filename);
        double span = qstep *0.5;
        int kmax = (qmax-qmin)/qstep+1;
        double q0 = qmin - span;
        double q1 = qmin + span;
        int i0 = avg_begin;
        for(int k = 0; k < kmax; k++){
            size_t *pix = pixel_index + i0;
            int len = avg_end[k] - i0;
            i0 = avg_end[k];
            double sum = 0;
            double sum2 = 0;
            int count = 0;
            for(int i = 0; i!=len; i++){
                size_t p = pix[i];
                double f = flat[p];
                if(f< 0.5 || f > 2 )continue;
                double my_bg =bg[p];
                double my_data = data[p];
                double v = my_data * inv_scale - my_bg;
                v *= f;
                v *= center_corr[p];
                sum += v;
                sum2 += v*v;
                count++;
            }
            increment_q();
        }
    }
    
    
    // fit against "int" data
    
    double fit(double q0, double q1, int *dat){
        int i0 = 0;
        while (center_qb[pixel_index[i0]*2] < q0) i0++;
        int i1 = i0;
        while (center_qb[pixel_index[i1]*2] < q1) i1++;
        int len = i1-i0;
        double *A = new double [len * 4];
        int *pivot = new int[5];
    }
    
    double fit1(double q0, double q1, int *dat, int *bg, float *flat){
        cout<<"going to fit"<<endl;
        int i0 = 0;
        cout << center_qb[pixel_index[i0]*2]<<endl;
        while (center_qb[pixel_index[i0]*2] < q0) i0++;
        cout << i0<<endl;
        int i1 = i0;
        while (center_qb[pixel_index[i1]*2] < q1) i1++;
        cout << i1 <<endl;
        int len = i1-i0;
        double *A = new double [len * 4];
        int *pivot = new int[5];
        double *y = A + len *3;
        size_t *pix = pixel_index + i0;
        double *A1 = A+len;
        double *A2 = A1+len;
  //      cout <<" len "<<len<<endl;;
        cout << "open"<<endl;
        ofstream fo("temp");
        //       ofstream fo("/Users/yoshi/temp/dump.txt");
        for(int i = 0; i!=len; i++){
            size_t p = pix[i];
            double flat_Value = flat[p];
            //            double flat_Err = flat_err[p];
            
            if(isnan(flat_Value) ||flat_Value < .90 || flat_Value > 1.1  || mask[p]){
                A[i] = y[i] = A1[i] = A2[i] = 0;
                //         cout << "isnan "<<i<<endl;
                continue;
            }else{
                double my_bg =bg[p];
                double my_data = dat[p];
                //                double my_sum = my_bg + my_data;
                
                A[i] = my_data;
                y[i] = my_bg;
                //              fo << i <<" "<<center_qb[p*2] <<" "<<center_qb[p*2+1]<<" "<<A[i]<<" "<<y[i];
                A1[i] = 1./flat_Value;
                A2[i] = center_qb[p*2]/flat_Value;
                fo <<center_qb[p*2]<<" "<< center_qb[p*2+1]<< " "<<my_data << " "<< my_bg<<" "<<flat_Value<<endl;
                //#define AAA
            }
        }
        cout <<"dgels "<<endl;
        int m = 3;
        int one = 1;
        double rcond = 1.e-300;
        int rank1;
        int info;
        double  iwork=-1;
        double *singular= new double[len];
        int *intwork = new int[len];
        int lwork = -1;
        char N = 'N';
        dgelss_(&len, &m, &one, A, &len,
                y, &len, singular, &rcond,
                &rank1, &iwork, &lwork, &info);
        //        dgels_(&N,&len, &m, &one, A, &len,
        //          y, &len, &iwork, &lwork, &info);        cout << iwork <<endl;
        lwork = iwork;
        double *work = new double [lwork];
        dgelss_(&len, &m, &one, A, &len,
                y, &len, singular, &rcond,
                &rank1, work, &lwork, &info);
        //        dgels_(&N, &len, &m, &one, A, &len,
        //          y, &len, work, &lwork, &info);
        cout << y[0] <<" "<< y[1]<<" "<<y[2]<<" "<< rank1<<endl;
        delete [] work;
        delete [] pivot;
        delete [] intwork;
        delete [] singular;
        delete [] A;
        return y[0];
    }
    
    // fit against cleaned "double" data.
    double fit(double q0, double q1){
        double scale =  fit(q0, q1, Detector->clean_data, Detector->clean_bg, Detector->response);
   //     cout <<"scale = "<<scale <<endl;
        return scale;
    }
    double fit(double q0, double q1, double *dat, double *bg, float *flat){
        int i0 = 0;
        while (center_qb[pixel_index[i0]*2] < q0) i0++;
        int i1 = i0;
        while (center_qb[pixel_index[i1]*2] < q1) i1++;
        int len = i1-i0;
        double *A = new double [len * 4];
        int *pivot = new int[5];
        double *y = A + len *3;
        size_t *pix = pixel_index + i0;
        double *A1 = A+len;
        double *A2 = A1+len;
   //      cout <<" len "<<len<<endl;;
        //     ofstream fo("/Users/yoshi/temp/dump.txt");
        for(int i = 0; i!=len; i++){
            size_t p = pix[i];
            double flat_Value = flat[p];
            //        double flat_Err = flat_err[p];
            
            if(isnan(flat_Value) ||flat_Value < .90 || flat_Value > 1.1){
                A[i] = y[i] = A1[i] = A2[i] = 0;
                //         cout << "isnan "<<i<<endl;
                continue;
            }else{
                double my_bg =bg[p];
                double my_data = dat[p];
                double my_sum = my_bg + my_data;
                if(isnan(my_bg) || isnan(my_data)){
                    A[i] = y[i] = A1[i] = A2[i] = 0;
                    continue;
                }
                A[i] = my_data ;//* flat_Value;
                y[i] = my_bg;// * flat_Value;
   //             fo << i <<" "<<center_qb[p*2] <<" "<<center_qb[p*2+1]<<" "<<A[i]<<" "<<y[i];
                A1[i] = 1;
                A2[i] = center_qb[p*2];
            }
        }
        cout <<"dgelsy "<<endl;
        int m = 3;
        int one = 1;
        double rcond = 1e-32;
        int rank1;
        int info;
        double  iwork=-1;
        int lwork = -1;
        dgelsy_(&len, &m, &one, A, &len,
                y, &len, pivot, &rcond,
                &rank1, &iwork, &lwork, &info);
  //      cout << iwork <<endl;
        lwork = iwork;
        double *work = new double [lwork];
        dgelsy_(&len, &m, &one, A, &len,
                y, &len, pivot, &rcond,
                &rank1, work, &lwork, &info);
//        cout << y[0] <<" "<< y[1]<<" "<<y[2]<<" "<< rank1<<endl;
        double scale = y[0];
        delete [] work;
        delete [] pivot;
        delete [] A;
        return scale;
    }
    
    inline void init_integrate(double qmin, double qmax, double qstep_in)
    {
        qstep = qstep_in;
        span = qstep * 0.5;
        kmax = (qmax-qmin)/qstep + 1;
        q0 = qmin - span;
        q1 = qmin + span;
        current_q = qmin;
        
        if(avg_end) delete [] avg_end;
        avg_end = new int[kmax];
        
        int i0 = 0;
        while (center_qb[pixel_index[i0]*2] < q0 && i0 < Detector->image_size) i0++;
        avg_begin=i0;
        for(int k = 0; k < kmax; k++){
            while (center_qb[pixel_index[i0]*2] < q1 && i0 < Detector->image_size) i0++;
            increment_q();
            avg_end[k] = i0;
  //          cout << i0<<endl;
        }
        check_flat();
    }
    
    

protected:
    void reset(){sum = 0; sum2 = 0; count=0;}
    inline void finish_average(){
        if(count){
            sum /= count;
            sum2 /=count;
            stdev = sqrt(sum2-sum*sum);
            avg =sum;
        }
    }
    inline void add_in(double val){
        sum += val;
        sum2 += val*val;
        count++;
    }
    inline void write_in()
    {
        write_float(&current_q);
        write_float(&avg);
        write_float(&stdev);
        write_int(&count);
    }
    inline void copy_plt()
    {
        plt_ptr->q = current_q;
        plt_ptr->avg = avg;
        plt_ptr->stdev = stdev;
        plt_ptr->count = count;
        plt_ptr++;
    }
    
    inline void copy_avg()
    {
        avg_ptr->avg = avg;
        avg_ptr->stdev = stdev;
        avg_ptr->count = count;
        avg_ptr++;
    }
    inline void write_float(float *dat){fo.write(reinterpret_cast<char *>(dat), sizeof(float));}
    inline void write_int(int *dat){fo.write(reinterpret_cast<char *>(dat), sizeof(int));}
    inline void increment_q()
    {
        q0 += qstep;
        q1 += qstep;
        current_q += qstep;
    }
    
    
    inline void init_integrate(double qmin, double qmax, double qstep_in, plot_data *plt)
    {
        init_integrate(qmin, qmax, qstep_in);
        plt_ptr = plt;
        do_write = false;
    }
    
    inline void init_plot(plot_data *plt)
    {
        plt_ptr = plt;
        do_write = false;
    }
    
    inline void init_file(const char filename[])
    {
        fo.open(filename);
        do_write=true;
    }
    
    inline void init_average(average_data *avg)
    {
        avg_ptr = avg;
    }
    
    double *work;
    int *data;
    int *bg;
    char *mask;
    plot_data *plt_ptr;
    
    average_data *avg_ptr;
    int *avg_end = NULL;//last pixel
    int avg_begin;
    
    float scale;
    double sum;
    double sum2;
    float avg;
    float stdev;
    int count;
    
    float *response = NULL;
    float current_q;
    double qstep;
    double q0, q1;
    double span;
    int kmax;
    ofstream fo;
    bool do_write;
    float *correct_response;
    
    float *qb = NULL;
    float *center_qb = NULL;
    float *center_corr = NULL;
    size_t *pixel_index = NULL;
    size_t valid_count;
    int dim0;
    map<string, double> poni;
    rotator<float> r;
    float center[3];
    double wavelength;
    detector *Detector;
};


#endif /* apply_poni2003_h */
