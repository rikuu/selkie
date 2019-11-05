# Selkie

Space-efficient spaced (l,k)-mer index for computing overlaps between raw
optical mapping data (Rmaps) that are used to assemble consensus optical maps.

## Getting started

```sh
git clone --recursive https://github.com/rikuu/selkie
cd selkie && make

git clone valouev
make

./selkie-index -o index ecoli-2000.valouev
./selkie-candidates ecoli-2000.valouev index > candidates.txt
valouev/ovlp/ovlp2 ecoli-2000.valouev candidates.txt
```
