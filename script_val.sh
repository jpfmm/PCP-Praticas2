export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/papi/5.3.2/lib/

export OMPP_APPNAME=kxx

export OMPP_CTR1=PAPI_L2_DCW
export OMPP_CTR2=PAPI_L3_DCW
export OMPP_CTR3=PAPI_L2_DRC
export OMPP_CTR4=PAPI_L1_DCM

for th in 1 2 4
do 
echo export OMP_NUM_THREADS=$th
export OMP_NUM_THREADS=$th

export OMPP_APPNAME=kkx-8000x8
echo ./omp_mat_vect_rand_split $th 8000 8
vallgrind --tool=callgrind --cache-sim=yes --separate-threads=yes ./omp_mat_vect_rand_split $th 8000 8

export OMPP_APPNAME=kkx-800x800
echo ./omp_mat_vect_rand_split $th 800 800
vallgrind --tool=callgrind --cache-sim=yes --separate-threads=yes ./omp_mat_vect_rand_split $th 800 800

export OMPP_APPNAME=kkx-8x8000
echo ./omp_mat_vect_rand_split $th 8 8000
vallgrind --tool=callgrind --cache-sim=yes --separate-threads=yes ./omp_mat_vect_rand_split $th 8 8000
done