# Selkie

Space-efficient spaced (l,k)-mer index for computing overlaps between raw
optical mapping data (Rmaps) that are used to assemble consensus optical maps.

## Getting started

In addition to everything included, Selkie also requires SDSL-lite.

```sh
git clone --recursive https://github.com/rikuu/selkie
cd selkie

cd thirdparty/valouev/ && make && cd ..
cd thirdparty/streamvbyte/ && make && cd ..

make

./selkie-index -o index thirdparty/valouev/test/ecoli-2000.valouev
./selkie-candidates thirdparty/valouev/test/ecoli-2000.valouev index > candidates.txt
valouev/ovlp/ovlp2 thirdparty/valouev/test/ecoli-2000.valouev candidates.txt ecoli-2000.ovlps ecoli-2000-detailed.ovlps
```
