#!/bin/bash
g++ -O3 lj3d.cpp

cnt=0
rho1=(5.0 4.0)
rho2=(6.4 0.008)

for x in liq gas; do
echo ${x^^}
########################
echo 1st run
time ./a.out << _EOT_
ini
scale
0.5
${rho1[cnt]}
_EOT_
echo
cp restart.dat start.dat
########################
echo 2nd run
./a.out << _EOT_
restart
scale
0.5
${rho2[cnt]}
_EOT_
echo
cp restart.dat start.dat
########################
echo Fin run
./a.out << _EOT_
restart
no
_EOT_
echo
cp restart.dat ${x}.rst
cp coord.xyz   ${x}.xyz
cp monitor.dat ${x}.out
cnt=$((cnt + 1))
echo
done

g++ rdf.cpp msd.cpp ana.cpp -o b.out
./b.out
