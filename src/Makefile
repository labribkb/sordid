sordid: SORDID.cpp randomkit.c mystructs.cpp overlap_utils.cpp utils.cpp Monitorer.cpp Monitorer.hpp
	g++-11 -std=c++17 -O3 -funroll-loops SORDID.cpp utils.cpp randomkit.c mystructs.cpp Monitorer.cpp overlap_utils.cpp -o sordid -DCGAL_USE_GMP -DCGAL_USE_GMPXX -DCGAL_USE_MPFR -lgmp -lboost_thread -fopenmp
all : sordid
.PHONY : all
