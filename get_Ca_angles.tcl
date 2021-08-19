# protein must be only one chain
# # of strand must be even
# 1. all inside resid in structure E
# 2. all pairs from Hbond between name N O and resids between two pairs
# both resids are inside, pair => inside pair
# calculate angle between pairs
# for furture pdb files using, change barrel define
# be careful of beta sheet which does not belong to berrel.

proc angle_inxypla {indexes} {
 set atom1_xy [lreplace [measure center [atomselect top "index [lindex $indexes 0]"]] 2 2]
 set atom2_xy [lreplace [measure center [atomselect top "index [lindex $indexes 1]"]] 2 2]
 set atom3_xy [lreplace [measure center [atomselect top "index [lindex $indexes 2]"]] 2 2]
# set atom2_xy [[atomselect top "index [lindex $indexes 1]"] get {x y}]
 set vector_1 [vecsub $atom1_xy $atom2_xy ]
 set vector_2 [vecsub $atom3_xy $atom2_xy ]
 set angle [expr acos( [vecdot $vector_1 $vector_2]/ ( [veclength $vector_1] * [veclength $vector_2] ) ) /3.1415 *180 ]
 return $angle 
}

proc angles_list {} {
##core code
set Ns [atomselect top "name N and protein"]
set Os [atomselect top "name O and protein"]
set indices [measure hbonds 4.0 45 $Ns $Os]
set length [llength [lindex $indices 0]]
set atoms [atomselect top "index [lindex $indices 0] [lindex $indices 1]"]
set chains_resids [lsort -unique [$atoms get {chain resid}]]
set pairs {}
for {set i 0} {$i < $length}  {incr i} {
  set N_chain_resid [lindex [[atomselect top "index [lindex [lindex $indices 0] $i]"] get {chain resid}] 0]
  set O_chain_resid [lindex [[atomselect top "index [lindex [lindex $indices 1] $i]"] get {chain resid}] 0]
  set pair {}
  lappend pair $N_chain_resid
  lappend pair $O_chain_resid
  if {[lindex $N_chain_resid 1] == [lindex $O_chain_resid 1]} {
   error "Error" "a pair have same resid" 400}
  if { [lindex $pair {0 0}] == [lindex $pair {1 0}]} {
   set pair_in_order [lsort -increasing -index 1 $pair] 
  } else {set pair_in_order [lsort -index 0 $pair] }
  lappend pairs $pair_in_order
}
#pairs show resid in pairs
set pairs [lsort -unique $pairs]
set length_pairs [llength $pairs]
set whole_pairs {}
#find new pairs
#gap 1
foreach pair $pairs {
  lappend whole_pairs $pair
  foreach i {1 -1} {
    set chain_resid_far_1 {}
    set chain_resid_far_2 {}
    set resid_far_1 [expr [lindex $pair {0 1}] + $i * 2]
    set resid_far_2 [expr [lindex $pair {1 1}] - $i * 2]
    set chain_far_1 [lindex $pair {0 0}]
    set chain_far_2 [lindex $pair {1 0}] 
    lappend chain_resid_far_1 $chain_far_1
    lappend chain_resid_far_1 $resid_far_1
    lappend chain_resid_far_2 $chain_far_2
    lappend chain_resid_far_2 $resid_far_2
    set pair_far {}
    lappend pair_far $chain_resid_far_1
    lappend pair_far $chain_resid_far_2
    if { [lindex $pair_far {0 0}] == [lindex $pair_far {1 0}]} {
     set pair_far [lsort -increasing -index 1 $pair_far] 
    } else {set pair_far [lsort -index 0 $pair_far] }
    if {[lsearch -exact $pairs $pair_far] >= 0  && $pair != $pair_far && [[atomselect top "resid [expr $resid_far_1 - $i * 2] [expr $resid_far_1 - $i] $resid_far_1 and structure E and name CA and chain $chain_far_1"] num] == 3 && [[atomselect top "resid [expr $resid_far_2 + $i * 2] [expr $resid_far_2 + $i] $resid_far_2 and structure E and name CA and chain $chain_far_2"] num] == 3} {
      set new_pair_resid_1 {}
      lappend new_pair_resid_1 $chain_far_1
      lappend new_pair_resid_1 [expr [lindex $pair {0 1}] + $i]
      set new_pair_resid_2 {}
      lappend new_pair_resid_2 $chain_far_2
      lappend new_pair_resid_2 [expr [lindex $pair {1 1}] - $i]
      set new_pair {}
      lappend new_pair $new_pair_resid_1
      lappend new_pair $new_pair_resid_2
      if { [lindex $new_pair {0 0}] == [lindex $new_pair {1 0}]} {
       set new_pair [lsort -increasing -index 1 $new_pair]     
      } else {set new_pair [lsort -index 0 $new_pair] }
      lappend whole_pairs $new_pair
      lappend chains_resids $new_pair_resid_1
      lappend chains_resids $new_pair_resid_2
    }
  }
}
set chains_resids [lsort -unique $chains_resids]
set whole_pairs [lsort -unique $whole_pairs]

set inside_chains_resids {}
#find inside_chains_resids
set barrel [atomselect top "name CA and protein"]
set centerPosXY [lreplace [measure center $barrel] 2 2]
foreach chain_resid $chains_resids {
  if {[[atomselect top "chain [lindex $chain_resid 0] and resid [lindex $chain_resid 1] and name CA"] get resname] == "GLY"} {
    set resSC [atomselect top "chain [lindex $chain_resid 0] and resid [lindex $chain_resid 1] and name HA1"]
    set resBB [atomselect top "chain [lindex $chain_resid 0] and resid [lindex $chain_resid 1] and name HA2"]} else {
    set resSC [atomselect top "chain [lindex $chain_resid 0] and resid [lindex $chain_resid 1] and not backbone"] 
    set resBB [atomselect top "chain [lindex $chain_resid 0] and resid [lindex $chain_resid 1] and backbone"]
  }
  set resBBV [lreplace [measure center $resBB] 2 2]
  set resSCV [lreplace [measure center $resSC] 2 2]
  set resBBVL [veclength [vecsub $resBBV ${centerPosXY} ] ] 
  set resSCVL [veclength [vecsub $resSCV ${centerPosXY} ] ]
  if {$resSCVL < $resBBVL} {
  lappend inside_chains_resids $chain_resid
  }
}

set array1 {}
foreach i $inside_chains_resids {
 lappend array1 [lindex $i 1]}

#add lost GLY
#2 missing gap
set missed_chains_resids {}
foreach chain_resid $inside_chains_resids {
  set chain_resid_2fromNow {}
  lappend chain_resid_2fromNow [lindex $chain_resid 0]
  lappend chain_resid_2fromNow [expr [lindex $chain_resid 1] + 2]
  set chain_resid_6fromNow {}
  lappend chain_resid_6fromNow [lindex $chain_resid 0]
  lappend chain_resid_6fromNow [expr [lindex $chain_resid 1] + 6]
  if { [lsearch -exact $inside_chains_resids $chain_resid_2fromNow]< 0 && [lsearch -exact $inside_chains_resids $chain_resid_6fromNow] >= 0 && [[atomselect top "chain [lindex $chain_resid 0] and resid [lindex $chain_resid 1] to [expr [lindex $chain_resid 1] + 6] and name CA"] num] == 7 } { lappend missed_chains_resids $chain_resid_2fromNow}
}
foreach chain_resid $missed_chains_resids {
  lappend inside_chains_resids $chain_resid
}
#1 missing gap
set missed_chains_resids {}
foreach chain_resid $inside_chains_resids {
  set chain_resid_2fromNow {}
  lappend chain_resid_2fromNow [lindex $chain_resid 0]
  lappend chain_resid_2fromNow [expr [lindex $chain_resid 1] + 2]
  set chain_resid_4fromNow {}
  lappend chain_resid_4fromNow [lindex $chain_resid 0]
  lappend chain_resid_4fromNow [expr [lindex $chain_resid 1] + 4]
  if { [lsearch -exact $inside_chains_resids $chain_resid_2fromNow]< 0 && [lsearch -exact $inside_chains_resids $chain_resid_4fromNow] >= 0 && [[atomselect top "chain [lindex $chain_resid 0] and resid [lindex $chain_resid 1] to [expr [lindex $chain_resid 1] + 4] and name CA"] num] == 5 } { lappend missed_chains_resids $chain_resid_2fromNow}
}
foreach chain_resid $missed_chains_resids {
  lappend inside_chains_resids $chain_resid
}

set array2 {}
foreach i $inside_chains_resids {
 lappend array2 [lindex $i 1]}

#delete wrong resid
# 123 remove 2
foreach chain_resid $inside_chains_resids {
  set chain_resid_p1fromNow {}
  lappend chain_resid_p1fromNow [lindex $chain_resid 0]
  lappend chain_resid_p1fromNow [expr [lindex $chain_resid 1] + 1]
  set chain_resid_n1fromNow {}
  lappend chain_resid_n1fromNow [lindex $chain_resid 0]
  lappend chain_resid_n1fromNow [expr [lindex $chain_resid 1] - 1]

  if {[lsearch -exact $inside_chains_resids $chain_resid_p1fromNow] >= 0 && [lsearch -exact $inside_chains_resids $chain_resid_n1fromNow] >= 0} { 
  set idx [lsearch $inside_chains_resids $chain_resid]
  set inside_chains_resids [lreplace $inside_chains_resids $idx $idx]
  }
}
#1356 remove 6
foreach chain_resid $inside_chains_resids {
 foreach i {1 -1} {
  set chain_resid_1fromNow {}
  lappend chain_resid_1fromNow [lindex $chain_resid 0]
  lappend chain_resid_1fromNow [expr [lindex $chain_resid 1] + $i]
  set chain_resid_3fromNow {}
  lappend chain_resid_3fromNow [lindex $chain_resid 0]
  lappend chain_resid_3fromNow [expr [lindex $chain_resid 1] + $i *3]
  set chain_resid_5fromNow {}
  lappend chain_resid_5fromNow [lindex $chain_resid 0]
  lappend chain_resid_5fromNow [expr [lindex $chain_resid 1] + $i *5]

  if {[lsearch -exact $inside_chains_resids $chain_resid_1fromNow] >= 0 && [lsearch -exact $inside_chains_resids $chain_resid_3fromNow] >= 0 && [lsearch -exact $inside_chains_resids $chain_resid_5fromNow] >= 0} {
  set idx [lsearch $inside_chains_resids $chain_resid]
  set inside_chains_resids [lreplace $inside_chains_resids $idx $idx]
  }

 }
}

set array3 {}
foreach i $inside_chains_resids {
 lappend array3 [lindex $i 1]}

#inside_pairs
set inside_pairs {}
set errors_chains_resids {}
foreach pair $whole_pairs {
  if {[lsearch -exact $inside_chains_resids [lindex $pair 0]] >= 0} {
    if {[lsearch -exact $inside_chains_resids [lindex $pair 1]] >= 0} {
      lappend inside_pairs $pair
    } else {lappend errors_chains_resids $pair}
  }
}

#calculate angle
set amino_aicd_angles {}
foreach chain_resid $inside_chains_resids {
  set angle_set {}
  foreach pair $inside_pairs {
    if {[lsearch -exact $pair $chain_resid] >= 0} {
      set i [expr 1 - [lsearch -exact $pair $chain_resid]]
      lappend angle_set [lindex $pair $i]
    }
  }
  if {[llength $angle_set] == 2} {
    set terminal_index1 [[atomselect top "resid [lindex $angle_set {0 1}] and chain [lindex $angle_set {0 0}] and name CA"] get index]
    set terminal_index2 [[atomselect top "resid [lindex $angle_set {1 1}] and chain [lindex $angle_set {1 0}] and name CA"] get index]
    set center_index [[atomselect top "resid [lindex $chain_resid 1] and chain [lindex $chain_resid 0] and name CA"] get index]
    set angle_indexes {}
    lappend angle_indexes $terminal_index1
    lappend angle_indexes $center_index
    lappend angle_indexes $terminal_index2
    set resname [[atomselect top "resid [lindex $chain_resid 1] and chain [lindex $chain_resid 0] and name CA"] get resname]
    set angle [angle_inxypla $angle_indexes]
    set unit {}
    lappend unit $resname
    lappend unit $angle
    lappend unit $chain_resid
    lappend amino_aicd_angles $unit
  }
}
return $amino_aicd_angles
}





#set resids [$N get resid]

#foreach resid $resids {

##select N
#set N [atomselect top "name N and resid $resid and protein"]
##set N_position [$N get {x y z}]
##set HN [atomselect top "name HN and resid $resid and protein"]
##set HN_position [$HN get {x y z}]

#set Os [atomselect top "name O and within 4 of (name N and resid $resid and protein) and protein"]

#set O_resids [$Os get resid]
#  foreach O_resid $O_resids {
#    set O [atomselect top "name O and resid $O_resid and protein"]
#    set indices 
#  }
#}
