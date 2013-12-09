#vmd -dispdev none -eofexit -e xtc2dcd.tcl test.gro test.xtc test.dcd

# disable unnecessary options
display update off

# get file names from command line
set grofile [lindex $argv 0]
set traj [lindex $argv 1]
set dcdfile [lindex $argv 2]


display update off
mol load gro $grofile xtc $traj
animate write dcd $dcdfile beg 1 waitfor all
exit
