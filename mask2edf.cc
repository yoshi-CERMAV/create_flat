#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;
int hsize = 578;
int vsize = 960;
char header_char[512];
int main(int argc, char *argv[])
{
    snprintf(header_char, 512, "{\nEDF_DataBlockID = 0.Image.Psd ;\nEDF_BinarySize = %d ;\nEDF_HeaderSize =   %d ;\nByteOrder = LowByteFirst ;DataType = SingleByte ;\nDim_1 = %d ;\nDim_2 = %d\n Image = 0 ;HeaderID = EH:000000:000000:000000 ;\n;Size = %d ;\n", hsize*vsize, 512, hsize, vsize, hsize*vsize);
    int n =  strlen(header_char);  
//    for(int i = n-2; i < 510; i++) snprintf(header_char+i,1, " ") ;
    char header_char1[512];
    cout << n << endl;
    snprintf(header_char1,512, "%s%*s}\n", header_char, 509-n, " ") ;
    n =  strlen(header_char1);  
    ifstream fi(argv[1]);
    ofstream fo(argv[2]);
    char *data = new char[hsize*vsize];
    fi.read(data, hsize*vsize);
    fo.write(header_char1, 512);
    fo.write(data, hsize*vsize);
}
