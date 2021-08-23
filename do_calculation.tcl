source get_Ca_angles.tcl

set fr [open "/localdata/zzhang624/data/ompx/api_PDBs_OPM/list.csv" r] ;
set file_name [read $fr] ;
close $fr ;

set dirt {}

foreach pdb $file_name {
 mol new /localdata/zzhang624/data/ompx/PDBs_OPM/addH-toallAA/proteins/$pdb

##calculate and collect
 set angle_list [angles_list]
 foreach amino_acid $angle_list {
  set key [lindex $amino_acid 0]
  set angle [lindex $amino_acid 1]
  set chain_resid [lindex $amino_acid 2]
  if {[lsearch -exact $dirt $key] >= 0} { 
    puts $files($key) "$angle $chain_resid $pdb"
  } else {
    set files($key) [open "data/${key}.dat" w]
    lappend dirt $key
    puts $files($key) "$angle $chain_resid $pdb"
  }
  
 }
##delete dcd and psf
mol delete top
}

foreach key $dirt {
  close $files($key)
}

foreach key $dirt {
  set out_file [open "data/${key}.dat" r]
  set n($key) 0.0
  set sum($key) 0.0
  while {[gets $out_file f] >= 0} {
    set n($key) [expr $n($key) + 1.0]
    set sum($key) [expr $sum($key) + [lindex $f 0]]
  }
  close $out_file
  set average($key) [expr $sum($key) / $n($key)]
}

foreach key $dirt {
  set out_file [open "data/${key}.dat" r]
  set tem_sum 0.0
  while {[gets $out_file f] >= 0} {
    set tem_sum [expr $tem_sum + (([lindex $f 0] - $average($key)) * ([lindex $f 0] - $average($key))) ]
  }
  set stan_devi($key) [expr sqrt($tem_sum / ($n($key) - 1.0))]
  close $out_file
}

set out_file [open "angles.dat" w]
foreach key $dirt {
  puts $out_file "$key $average($key) $stan_devi($key) $n($key)"
}
close $out_file
quit
