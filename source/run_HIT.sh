#Shengjie Sun:
#The HIT program need user to provide target ions name and possible number: For example, SOD for sodium ions. which is the atom name in pdb file.
#Users can choose how many ions to select based on the occupancies of possible bound ions, which is provided in the generated txt file in the first running.


echo -e "\n"
echo "##########################################################################"
echo "##################   please cite the article   ###########################" 
echo "# Hybrid method for representing ions in implicit solvation calculations #"
echo "#######   S Sun, C Karki, Y Xie, Y Xian, W Guo, BZ Gao   #################"
echo "##########################################################################"
echo "########################    HIT IS RUNNING   #############################"

while getopts f:i:n:s: flag
do
	case $flag in
		f) filename=$OPTARG;;
		i) ionname=$OPTARG;;
		n) ionnumber=$OPTARG;;
		s) cubesize=$OPTARG;;
	esac
done
echo "target bound ion: $ionname";
echo "number of bound ions: $ionnumber";
if [ -z $cubesize ]
then
 cubesize=3.3
fi
echo "the cubesize: $cubesize"

grep END $filename | wc > frame.number
grep $ionname $filename > $ionname.template

if [ -s "$ionname.template" ]
then
awk '{print substr($0,31,26)}' $ionname.template | awk '{print $1,$2,$3}'> SSS.SSS
./HIT $cubesize  > $ionname.txt

awk '{printf("%8.3f%8.3f%8.3f\n",$1,$2,$3)}' $ionname.txt > Gr.iii
paste $ionname.template Gr.iii | head -$ionnumber > comb.iii
awk '{printf("%s%s%6.2f%6.2f     %s\n",substr($0,1,31),substr($0,81,105),$10,$11,$12)}' comb.iii > $ionname.pdb

rm Gr.iii
rm comb.iii
rm SSS.SSS
rm $ionname.template
rm frame.number
echo "success" > report.txt
echo "$ionname.pdb is the bound ions file" >> report.txt
echo "$ionname.txt contains the targe ions position and correspoding occupancy" >>report.txt
echo "SUCCESS"

else 
	rm frame.number
	rm $ionname.template
	echo "error: no target ions" > report.txt
echo "error:no target ions"
fi


