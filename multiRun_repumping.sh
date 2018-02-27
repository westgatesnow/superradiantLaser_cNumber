#This bash file is to run multiple simulations with respect to repumping
#for the superradiant laser using cNumber theory.

iFile=input.txt
nMax=10
init=0.1
interval=0.1

for ((i=0; i<nMax; i+=1)) 
do

w=$(echo "$init + $interval * $i" | bc -l)

if (( $(bc <<< "w < 1") ))
then
	printf "dt 0.01
	tmax 20
	nstore 100
	nTrajectory 40
	nAtom 100
	gammac 0.1
	repumping $w			
	name N100_repumping0${w}_Traj40" > $iFile
else
	printf "dt 0.01
	tmax 20
	nstore 100
	nTrajectory 40
	nAtom 100
	gammac 0.1
	repumping $w
	name N100_repumping${w}_Traj40" > $iFile
fi


./superradiantLaser -f $iFile

number=$((1+$i))

echo "Run ${number} of" $nMax

done
