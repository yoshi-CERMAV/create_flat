#include <unistd.h>
#include <detector.h>
#include "apply_poni2023.h"
#include <Read_D5WOS.h>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>

char read_h5::D5fmt_[] = "/%d.1/measurement/D5";
char read_h5::WOSfmt_[] = "/%d.1/measurement/WOS";
char read_h5::PM1fmt_[] = "/%d.1/measurement/pm1";
char read_h5::fmtdat_[] = "/%d.1/instrument/epoch/value";

hsize_t vec_sizeD5[2]={578, 960};
hsize_t vec_sizeWOS[2]={1156, 600};

DataSpace read_h5::memspaceD5 (2, vec_sizeD5) ;
DataSpace read_h5::memspaceWOS (2, vec_sizeWOS) ;
enum D5_WOS{D5, WOS};
using namespace std;
detector *Detector;
int *sample_data = 0;
int *mt_data = 0;
double scale = 1;
plot_data *plotdata_WOS;
average_data *Avg_data;
read_h5 *reader ;

double q0 = 0.25;
double q1 = 4;

int qsize;
double wos_range[3] = {0.25, 4, 0.01};

FILE *gp;
apply_poni *Poni = NULL;

template<class T>
void write(const char filename[], T *flat, int size){
    ofstream fo(filename);
    fo.write(reinterpret_cast<char *>(flat), size*sizeof(T));
}

void write(const char fmt[], const char filename[], average_data *ptr, int size)
{
    char outfile[255];
    snprintf(outfile, 255, fmt, filename);
    ofstream fo(outfile);
    fo.write(reinterpret_cast<char *>(ptr), sizeof(average_data)*size);
}
void write(const char fmt[], const char filename[], double *ptr, int size)
{
    char outfile[255];
    snprintf(outfile, 255, fmt, filename);
    ofstream fo(outfile);
    fo.write(reinterpret_cast<char *>(ptr), sizeof(double)*size);
}

void read(const char fmt[], const char filename[], double *ptr, int size)
{
    char infile[255];
    snprintf(infile, 255, fmt, filename);
    ifstream fi(infile);
    fi.read(reinterpret_cast<char *>(ptr), sizeof(double)*size);
}

static void save(const char fmt[], int num, plot_data *dat, int size)
{
    char filename [256];
    snprintf(filename, 256, fmt, num);
    ofstream fo(filename);
    fo.write(reinterpret_cast<char*>(dat), size*16);
}

void read_paths(const char filename[], vector<string> &name)
{
   ifstream fi(filename);
   string s;
   while(fi >>s) name.push_back(s); 
}

void dump_data(const char filename[], void *ptr, int image_size)
{
    ofstream fo(filename);
    fo.write(reinterpret_cast<char *>(ptr), image_size*4);
}

void finalize_flat(float *flat, char *mask, int n, float min, float  max)
{
    int count0 = 0;
    int count1 = 0;
    for(int i = 0; i != n; i++){
        if(flat[i]< min) {flat[i] = 0; mask[i] = 1;count0++;}
        if(flat[i]> max) {flat[i] = 0; mask[i] = 1;count1++;}
    }
    cout << "pixels under "<<min << ": " << count0<< endl;
    cout << "pixels above "<<max << ": " << count1<< endl;
}

template<class T>
void readline(ifstream &fi, T &val)
{
    char buffer[256];
    fi.getline(buffer,256);
    istringstream istr(buffer);
    istr >> val;
}

int main(int argc, char *argv[])
{
    int detector_type;
    char detector_name[256];
    char file_name[255];
    char poni_name[256];
    char flat_fname[256];
    char mask_fname[256];
    char negative_mask[256];
    char profile_fname[256];
    int snap_h, snap_v;
    char snap_fname[256];
    int sample_scan;
    int mt_scan;
    double qstep;
    double scale2;
    double mask_max;
    double mask_min;
    double image_max;
    double image_min;
    ifstream fi(argv[1]);
    readline(fi, detector_name);
    regex wos("WOS", regex_constants::icase);
    if(regex_search(detector_name, wos)) detector_type = D5_WOS::WOS;
    else detector_type=D5_WOS::D5;
    readline(fi, file_name);
    readline(fi, poni_name);
    readline(fi, qstep);
    readline(fi, flat_fname);
    readline(fi, mask_fname);
    readline(fi, profile_fname);
    readline(fi, snap_h);
    readline(fi, snap_v);
    readline(fi, snap_fname);
    readline(fi, sample_scan);
    readline(fi, mt_scan);
    readline(fi, scale2);
    readline(fi, mask_max);
    readline(fi, mask_min);
    readline(fi, image_max);
    readline(fi, image_min);
    cout << "detector_name : "<<detector_name<<endl;
    cout << "snap h : "<<snap_h<<endl;
    cout << "snap v : "<<snap_v<<endl;
    cout << "mask_min : "<<mask_min<<endl;
    cout << "mask_max : "<<mask_max<<endl;
    Detector = new detector(detector_name);
//    readline(fi, negative_mask);
//    cout << negative_mask;
//    Detector->load_negative_mask(negative_mask);
    reader = new read_h5(file_name);
    char buffer[256];
    
    cout << "opening detector "<<endl;
    cout << detector_name <<endl;
    cout << "detector created "<<endl;
    Poni = new apply_poni(Detector);
    Poni->read_poni(poni_name);

    int image_size;
    cout << sample_scan<<endl;
    if(detector_type==D5_WOS::WOS){
        reader->accumulate_scan_WOS(sample_scan, sample_data);
        image_size = vec_sizeWOS[0] * vec_sizeWOS[1];
    }else{
        reader->accumulate_scan_D5(sample_scan, sample_data);
        image_size = vec_sizeD5[0] * vec_sizeD5[1];
    }
    
    sample_data = new int[image_size];
    mt_data  = new int[image_size];
    dump_data("sample_data", sample_data, image_size);
    float pm = reader->accumulate_pm(sample_scan);
    cout << pm<<endl;
    
    cout << mt_scan<<endl;
    if(detector_type==D5_WOS::WOS){
        reader->accumulate_scan_WOS(mt_scan, mt_data);
    }else{
        reader->accumulate_scan_D5(mt_scan, mt_data);
    }
    dump_data("mt_data", mt_data, image_size);
    
    float pm1 = reader->accumulate_pm(mt_scan);
    float scale = pm / pm1;
    cout << "scale "<<scale<<endl;
    float *data = new float[image_size]();
    for(int i = 0; i < image_size; i++){
        if(!Detector->hide(i)) data[i] = sample_data[i] - scale2*scale * mt_data[i];
    }
    dump_data("subtracted_data", data, image_size);

    ///do a robust integration
    q1 = Poni->qmax();
    q0 = Poni->qmin();
    qsize = (q1- q0)/ qstep;\
    
    cout << "qsize = "<< qsize<<endl;
    Avg_data = new average_data[qsize];
    Poni->init_integrate(q0, q1, qstep);
    average_data *Avg_ptr = Avg_data;
    Poni->integrate(data, Avg_ptr);
    write("%s",profile_fname, Avg_data, qsize);

    ///create a flat
    float *flat = new float[image_size];
    Poni->make_flat(flat, q0, q1, qstep, Avg_data, data);
    char *mask = new char[image_size]();
    /// remove outliers and put into mask
    finalize_flat(flat, mask, image_size, mask_min, mask_max);
    
    write(flat_fname, flat, image_size);
    write(mask_fname, mask, image_size);
    int snap_size =snap_h*snap_v;
    unsigned char *detector_image = new unsigned char[snap_size];
    Detector->init_quick_view(snap_h, snap_v);
    Detector->make_image1(detector_image, flat, image_min, image_max);
    write(snap_fname, detector_image, snap_size);
    delete [] Avg_data;
}
