# Load the system and the trajectory
mol new AlB2/AlB2.xyz

#set up the atom selections
set sel1 [atomselect top "type A"]
set sel2 [atomselect top "type B"]

pbc set {10.0 8.66025403784 5.0}

#calculate g(r)
set gr [measure gofr $sel1 $sel2 delta 0.001 rmax 2.5 usepbc 1 selupdate 1 first 0]

#set up the outfile and write out the data
set outfile [open AlB2/gofr_AlB2_ab.dat w]

set r_ab [lindex $gr 0]
set gr_ab [lindex $gr 1]

foreach j $r_ab k $gr_ab {
   puts $outfile "$j $k"
}

close $outfile

set outfile [open AlB2/gofr_AlB2_aa.dat w]
set gr [measure gofr $sel1 $sel1 delta 0.001 rmax 2.5 usepbc 1 selupdate 1 first 0]

set r_aa [lindex $gr 0]
set gr_aa [lindex $gr 1]

foreach j $r_aa k $gr_aa {
   puts $outfile "$j $k"
}

close $outfile

set outfile [open AlB2/gofr_AlB2_bb.dat w]
set gr [measure gofr $sel2 $sel2 delta 0.001 rmax 2.5 usepbc 1 selupdate 1 first 0]

set r_bb [lindex $gr 0]
set gr_bb [lindex $gr 1]

foreach j $r_bb k $gr_bb {
   puts $outfile "$j $k"
}

close $outfile
exit
