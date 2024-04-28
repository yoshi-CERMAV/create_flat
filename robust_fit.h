//
//  robust_fit.h
//  
//
//  Created by Yoshiharu Nishiyama on 28/04/2024.
//

#ifndef robust_fit_h
#define robust_fit_h


#endif /* robust_fit_h */

class robust_fit:public class apply_poni
{
public:
    
    void init(int n){
        max_data = n;
        num_para = 3;
        lwork = 4 * num_para +1;
        work = new double [lwork];
        A = new double [max_data];
        A1 = new double [max_data];
        y = new double [max_data];
        y1 = new double [max_data];
    }
    void re_alloc(int len){
        max_data = len;
        delete [] y1;
        delete [] y;
        delete [] A1;
        delete [] A;
        A = new double [max_data];
        Aback = new double [max_data];
        y = new double [max_data];
        yback = new double [max_data];
    }
    void fit(double q0, double q1, int *dat)
    {
        int i0, i1;
        int len = find_i_range(q0, q1, i0, i1);
        if (len >  max_data){
            re_alloc(len);
        }
        size_t index = pixel_index+i0;
        double *A1 = A+len;
        double *A2 = A1+len;
        for(int i = 0; i != len; i++)
        {
            size_t p = index[i];
            double q = center_qb[p * 2];
            double flat_value = flat[p];
            double Y = dat[i]
            double sigma = 1./sqrt(Y);
            y[i] = sigma ;
            A[i] = 1;
            A1[i] = q;
            A2[i] = q*q;
        }
        
    }
protected:
    int max_data;
    int num_para;
    int lwork;
    double *work;
    double *A;
    double *Aback;
    double *y;
    double *yback;
}
