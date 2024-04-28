create_flat : create_flat.cc 
	h5c++ -o create_flat create_flat.cc -I./include  -std=c++11
flat2edf : flat2edf.cc
	c++ -o flat2edf flat2edf.cc
mask2edf : mask2edf.cc
	c++ -o mask2edf mask2edf.cc
all: create_flat flat2edf mask2edf
clean:
	rm create_flat mask2edf flat2edf mt_data sample_data subtracted_data mask_D5_carbon5.dat *.o *.edf carbon_flat5.dat carbon_flat_image_D5_5 carbon+air_D5_5.dat
