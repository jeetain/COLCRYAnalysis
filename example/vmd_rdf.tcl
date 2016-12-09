# Load the system and the trajectory
mol new text.xyz

#set up the atom selections
set sel1 [atomselect top "name A"]
set sel2 [atomselect top "name B"]

#calculate g(r)
set gr [measure gofr $sel1 $sel2 delta .01 rmax 5 usepbc 1 selupdate 1 first 1 last -1 step 1]

#set up the outfile and write out the data
set outfile [open gofr_A_B.dat w]

set r [lindex $gr 0]
set gr2 [lindex $gr 1]
set igr [lindex $gr 2]

set i 0
foreach j $r k $gr2 l $igr {
   puts $outfile "$j $k $l"
}

close $outfile

exit