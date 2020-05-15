g++ staticStackedPWAAbench.cpp -O3

echo "./a.out ../data/fb_200M_uint64 16"
./a.out ../data/fb_200M_uint64 16
echo "./a.out ../data/fb_200M_uint64 64"
./a.out ../data/fb_200M_uint64 64

echo "./a.out ../data/lognormal_200M_uint64 16"
./a.out ../data/lognormal_200M_uint64 16
echo "./a.out ../data/lognormal_200M_uint64 64"
./a.out ../data/lognormal_200M_uint64 64

echo "./a.out ../data/normal_200M_uint64 16"
./a.out ../data/normal_200M_uint64 16
echo "./a.out ../data/normal_200M_uint64 64"
./a.out ../data/normal_200M_uint64 64

echo "./a.out ../data/uniform_dense_200M_uint64 16"
./a.out ../data/uniform_dense_200M_uint64 16
echo "./a.out ../data/uniform_dense_200M_uint64 64"
./a.out ../data/uniform_dense_200M_uint64 64