nf=`cat nframes.dat`

for (( i=1; i<=$nf; i++ ))
do
    mv pos_${i} temp.dat
    cp pos_template.pov frame_${i}.pov
    povray -V -GA +A0.0 +H600 -J +R5 +W800 frame_${i}.pov
    FILE=frame_${i}.png
    if [ ! -f "$FILE" ]; then
      break
    fi
    rm -f frame_${i}.pov
    rm -f temp.dat
done
