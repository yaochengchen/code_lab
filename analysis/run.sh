#!/bin/bash
# get all filename in specified path
rm cycabc
g++ -o cycabc -O3 Reconstruction_tchain_new.cc  -L`root-config --libs` `root-config --cflags` `root-config --glibs`

#DATA_Dir=/media/cyc/For_Linux/TAROGE4_DATA/data
DATA_Dir=/media/cyc/1p9TB/TAROGE4_DATA/data

#echo ${DATA_Dir}
#sleep 5000
count=0
#files=$(ls ${DATA_Dir})
#for filename in ${DATA_Dir}/20211208
for filename in ${DATA_Dir}/2019123*
do
   rm ./run_list/${filename##*/}.txt
   rootfiles=$(ls ${filename})
   for rootname in $rootfiles
   do
       echo ${filename}/$rootname>>./run_list/${filename##*/}.txt
   done
   time ./cycabc ${filename##*/} ./run_list/${filename##*/}.txt &
   count=$(($count+1))
   if [ $count -gt 1 ]
   then
      echo $count
      sleep 1800
   fi
   
done

for filename in ${DATA_Dir}/2019121*
do
   rm ./run_list/${filename##*/}.txt
   rootfiles=$(ls ${filename})
   for rootname in $rootfiles
   do
       echo ${filename}/$rootname>>./run_list/${filename##*/}.txt
   done
   time ./cycabc ${filename##*/} ./run_list/${filename##*/}.txt &
   count=$(($count+1))
   if [ $count -gt 8 ]
   then
      echo $count
      sleep 1800
   fi
   
done

for filename in ${DATA_Dir}/2019122*
do
   rm ./run_list/${filename##*/}.txt
   rootfiles=$(ls ${filename})
   for rootname in $rootfiles
   do
       echo ${filename}/$rootname>>./run_list/${filename##*/}.txt
   done
   time ./cycabc ${filename##*/} ./run_list/${filename##*/}.txt &
   count=$(($count+1))
   if [ $count -gt 1 ]
   then
      echo $count
      sleep 1800
   fi
   
done













