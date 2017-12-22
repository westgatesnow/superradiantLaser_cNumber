#This bash file is to run multiple simulations with respect to nAtom
#for the superradiant laser using cumulant theory.

iFile=input.txt
nMax=1
init=10
interval=10

for ((i=0; i<nMax; i+=1)) 
do

nAtom=$(($init+$i*$interval))

#$(echo "10+$i*10" | bc -l);

printf "dt 0.01
tmax 20
nstore 100
nAtom $nAtom
gammac 0.1
repumping 10
name N${nAtom}_repumping10" > $iFile

./superradiantLaser -f $iFile

number=$((1+$i))

echo "Run ${number} of" $nMax

done
