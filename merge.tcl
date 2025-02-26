puts "Please input your psf file,My lord "
gets stdin PSF
puts "your input is $PSF, Now, please input your dcd file"
gets stdin DCD
puts "your input is $PSF and $DCD, please wait, My lord"

set id [mol new "$PSF"]
mol addfile "$DCD" waitfor all molid $id

set n_frames [molinfo $id get numframes]
for {set i 0} {$i< $n_frames} {incr i} { [atomselect top "ions within 10 of protein" frame $i] writepdb $i.ion}

eval exec cat 0.ion > merge.pdb
for {set i 1} {$i < $n_frames} {incr i} {eval exec cat $i.ion >> merge.pdb;eval exec rm $i.ion
}


