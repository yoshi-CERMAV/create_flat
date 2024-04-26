create_flat : create_flat.cc 
	h5c++ -o create_flat create_flat.cc -I./include  -std=c++11
flat2edf : flat2edf.cc
	c++ -o flat2edf flat2edf.cc
mask2edf : mask2edf.cc
	c++ -o mask2edf mask2edf.cc
all: create_flat flat2edf mask2edf
