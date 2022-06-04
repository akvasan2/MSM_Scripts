proc dist {mono kk_h2 vv_o don_resid acpt_resid frame} {
	set kk_h $kk_h2
	if {$don_resid == 256 || $don_resid == 312} {
		set kk_h {HN}
	} 
	if {$don_resid == 296} {
		set kk_h {HN}
	}
        set min_dist 1000000
	foreach h $kk_h {

		set don [atomselect top "segname $mono and resid $don_resid and name $h" frame $frame]

		set don_com [measure center $don]
		
		$don delete
		foreach o $vv_o {
			set acpt [atomselect top "segname $mono and resid $acpt_resid and name $o" frame $frame]
			#puts "ACPT: $acpt_resid $o [$acpt num]"
			set acpt_com [measure center $acpt]
			$acpt delete
			set dist [veclength [vecsub $don_com $acpt_com]]
			if {$dist < $min_dist} {
				set min_dist $dist
			} else {continue}
		}	
	}
	return $min_dist	
} 

set hydrogens [dict create ASP {HN} GLU {HN} ARG {HH11 HH12 HH22 HH21 HE} LYS {HZ1 HZ2 HZ3} TYR {HH} GLN {HE21 HE22} THR {HG1} GLY {HN}]


set oxygens [dict create GLU {OE1 OE2 O} ASP {OD1 OD2 O} THR {OG1 O} GLY {O} ALA {O} PHE {O} PRO {O} TYR {O OH}]

if {1} {
set DA_pair [dict create 132  {123 121} \
			167  {121 } \
			168  {121 122 } \
			42  {117} \
			82  {121 } \
			16  {117} \
			80  {121} \
			22  {117 118 } \
			310  {117} \
			294  {121 118 } \
			40  {117 } \
			302  {117 } \
			32  {121 } \
			121  {294 } \
			312  {116 } \
			296  {117 118 } \
			119  {22 294 } \
			120  {294 32} \
			122  {256}]
}

set residues [dict create 42 ARG 82 ARG 132 ARG 167 ARG 168 ARG \
				  16 LYS 80 LYS 305 LYS \
				  22 TYR 310 TYR 294 TYR 40 TYR 32 TYR 302 TYR \
				  339 GLN \
				  117 GLU 296 GLU \
				  121 ASP 312 ASP 256 ASP \
				  119 ASP 120 GLY \
				  122 THR \
				  123 ALA \
				  118 PHE \
				  116 PRO ]


for {set run 1} {$run<=5} {incr run} {
	mol load psf ../../build/protein.psf
	mol addfile ../../Run$run/equil_ForAnalysis_total.dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
	set frame [molinfo top get numframes]
	
	foreach mon {PROA PROB PROC} {
	set file1 [open "Run_[expr $run-1]-$mon.dat" w]
                for {set i 0} {$i < $frame} {incr i} {
		        dict for {kk vv} $DA_pair {

			        set kk_res [dict get $residues $kk]
			        set kk_h [dict get $hydrogens $kk_res]

			        foreach v $vv {
				        set vv_res [dict get $residues $v]
				        set vv_o [dict get $oxygens $vv_res]

					set distance [dist $mon $kk_h $vv_o $kk $v $i]
					puts -nonewline $file1 "$distance \t"	

					}
				}
                        puts $file1 ""
		}
        close $file1
	}
mol delete all
}

set file2 [open "DA_pair-list.txt" w] 

dict for {kk vv} $DA_pair {
        set kk_res [dict get $residues $kk]
                                                                 
        foreach v $vv {
	        set vv_res [dict get $residues $v]
		puts  $file2 "$kk $v"		
		}
}
close $file2

exit


