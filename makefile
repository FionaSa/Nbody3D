all: nbody.g nbody1icc nbody1gcc nbody1clang nbody2icc nbody2gcc nbody2clang nbody3icc nbody3gcc nbody3clang nbody4gcc nbody4icc nbody4clang nbody5clang nbody5gcc nbody5icc nbody6clang nbody6gcc nbody6icc

execute: nbody.g nbody1icc nbody1gcc nbody1clang nbody2icc nbody2gcc nbody2clang nbody3icc nbody3gcc nbody3clang nbody4gcc nbody4icc nbody4clang nbody5clang nbody5gcc nbody5icc nbody6clang nbody6gcc nbody6icc
	taskset -c 6 ./nbody.g
	taskset -c 6 ./nbody1icc
	taskset -c 6 ./nbody1gcc
	taskset -c 6 ./nbody1clang
	taskset -c 6 ./nbody2gcc
	taskset -c 6 ./nbody2icc
	taskset -c 6 ./nbody2clang
	taskset -c 6 ./nbody3icc
	taskset -c 6 ./nbody3gcc
	taskset -c 6 ./nbody3clang
	taskset -c 6 ./nbody4gcc
	taskset -c 6 ./nbody4icc
	taskset -c 6 ./nbody4clang
	taskset -c 6 ./nbody5gcc
	taskset -c 6 ./nbody5icc
	taskset -c 6 ./nbody5clang	
	taskset -c 6 ./nbody6gcc
	taskset -c 6 ./nbody6icc
	taskset -c 6 ./nbody6clang		

executeicx: nbody1icx nbody2icx nbody3icx nbody4icx nbody5icx nbody6icx
	taskset -c 6 ./nbody1icx
	taskset -c 6 ./nbody2icx
	taskset -c 6 ./nbody3icx
	taskset -c 6 ./nbody4icx
	taskset -c 6 ./nbody5icx
	taskset -c 6 ./nbody6icx

nbody.g: nbody.c
	icc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

nbody1icc: nbody1.c
	icc -march=native -mavx2 -Ofast -fopt-info-all=nbody1.gcc.optrpt $< -o $@ -lm -fopenmp
nbody1icx: nbody1.c
	icx -march=native -mavx2 -Ofast $< -o $@ -lm -fopenmp

nbody1gcc: nbody1.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody1.gcc.optrpt $< -o $@ -lm -fopenmp
nbody1clang: nbody1.c
	clang -march=native -mavx2 -Ofast $< -o $@ -lm -fopenmp

nbody2icx: nbody2.c
	icx -march=native -mtune=native -msse2 -mavx2 -Ofast  $< -o $@ -lm -fopenmp
nbody2icc: nbody2.c
	icc -march=native -mtune=native -msse2 -mavx2 -Ofast -fopt-info-all=nbody2.gcc.optrpt $< -o $@ -lm -fopenmp
nbody2gcc: nbody2.c
	gcc -march=native -mtune=native -msse2 -mavx2 -Ofast -fopt-info-all=nbody2.gcc.optrpt $< -o $@ -lm -fopenmp
nbody2clang: nbody2.c
	clang -march=native -mtune=native -msse2 -mavx2 -Ofast $< -o $@ -lm -fopenmp

nbody3icx: nbody3.c
	icx -march=native -mtune=native -msse2 -mavx2 -Ofast  $< -o $@ -lm -fopenmp
nbody3icc: nbody3.c
	icc -march=native -mtune=native -msse2 -mavx2 -Ofast -fopt-info-all=nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
nbody3gcc: nbody3.c
	gcc -march=native -mtune=native -msse2 -mavx2 -Ofast -fopt-info-all=nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
nbody3clang: nbody3.c
	clang -march=native -mtune=native -msse2 -mavx2 -Ofast $< -o $@ -lm -fopenmp

nbody4icx: nbody4.c
	icx -march=native -mtune=native -msse2 -mavx2 -Ofast $< -o $@ -lm -fopenmp
nbody4icc: nbody4.c
	icc -march=native -mtune=native -msse2 -mavx2 -Ofast -ftree-vectorize -finline-functions -funroll-loops -foptimize-sibling-calls  -foptimize-strlen -fopt-info-all=nbody4.gcc.optrpt -fgcse-sm $< -o $@ -lm -fopenmp
nbody4gcc: nbody4.c
	gcc -march=native -mtune=native -msse2 -mavx2 -Ofast -ftree-vectorize -finline-functions -funroll-loops -foptimize-sibling-calls  -foptimize-strlen -fopt-info-all=nbody4.gcc.optrpt -fgcse-sm $< -o $@ -lm -fopenmp
nbody4clang: nbody4.c
	clang -march=native -mtune=native -msse2 -mavx2 -Ofast  $< -o $@ -lm -fopenmp


nbody5icx: nbody5.c
	icx -march=native -mtune=native -mavx2 -Ofast  $< -o $@ -lm -fopenmp
nbody5icc: nbody5.c
	icc -march=native -mtune=native -mavx2 -Ofast --ffp-model=precise -ftree-vectorize -finline-functions -funroll-loops -fopt-info-all=nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
nbody5gcc: nbody5.c
	gcc -march=native -mtune=native -mavx2 -Ofast -ftree-vectorize -finline-functions -funroll-loops  -fopt-info-all=nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
nbody5clang: nbody5.c
	clang -march=native -mtune=native -mavx2 -Ofasts  $< -o $@ -lm -fopenmp


nbody6icx: nbody6.c
	icx -march=native -mtune=native -msse2 -mavx2 -Ofast  $< -o $@ -lm -fopenmp
nbody6icc: nbody6.c
	icc -march=native -mtune=native -msse2 -mavx2 -Ofast --ffp-model=precise -fopt-info-all=nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
nbody6gcc: nbody6.c
	gcc -march=native -mtune=native -msse2 -mavx512f -Ofast  -fopt-info-all=nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
nbody6clang: nbody6.c
	clang -march=native -mtune=native -msse2 -mavx512f -Ofast  $< -o $@ -lm -fopenmp


clean:
	rm -Rf *~ *.g *.i *.optrpt

