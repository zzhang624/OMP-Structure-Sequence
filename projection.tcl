set fr [open "/localdata/zzhang624/data/ompx/PDBs_OPM/ModPDBsv2/files_name" r] ;
set file_name [read $fr] ;
close $fr ;


foreach pdb $file_name {

        mol new /localdata/zzhang624/data/ompx/PDBs_OPM/ModPDBsv2/$pdb
        set sel1 [atomselect top "name CA"]
      
          set f [open xy_dat/${pdb}_xy.dat w]
          set x_set [$sel1 get x]
          set y_set [$sel1 get y]
          foreach x $x_set y $y_set {
            puts $f "$x $y"
          }
        close $f
        mol delete top
}
quit
