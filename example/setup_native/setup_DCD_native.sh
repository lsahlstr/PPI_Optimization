#!/usr/bin/env bash

mkdir -p DCD

for i in {0..4}
do

    echo $i
    mkdir -p DCD/RUN${i}
    cd DCD/RUN${i}
    ln -s /export/nVerde/users/lsahlstr/go_interact/$1/gopair/optimize/setup_native/n.init${i}.dcd ./${i}.dcd
    cd ../../ 

done
