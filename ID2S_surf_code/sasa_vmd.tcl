set OUTPUT [open "ID2S_new.txt" w]

# Get list of residues 
set allsel [atomselect top protein]
set residlist [lsort -unique [$allsel get resid]]

# Make atom selection, calculate sasa for each residue 

foreach r $residlist {
	set sel [atomselect top "resid $r"]
	set rsasa [measure sasa 1.4 $allsel -restrict $sel]
	puts $OUTPUT "residue: $r sasa: $rsasa"
}
close $OUTPUT
quit
