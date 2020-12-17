# Load file.
n=0
while read fline
do
        n=$((n + 1))
        filearray[$n]=$fline
done < /home/s4430291/chanlab-genomics/jackknifing/jf_scripts/files2runH3H4.txt