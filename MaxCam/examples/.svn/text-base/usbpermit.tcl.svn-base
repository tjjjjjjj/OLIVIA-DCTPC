#!/opt/apogee/bin/tclsh8.3

set uall [split [exec /sbin/lsusb] \n]
foreach d $uall {
   set id [lindex $d 5]
   if { $id == "125c:0010" } {
      set path [string trim "/proc/bus/usb/[lindex $d 1]/[lindex $d 3]" :]
      puts stdout "Found Apogee ALTA-U : $d, path $path"
      puts stdout "Invoking sudo chmod to alter rw device permissions..."
      puts stdout "sudo chmod a+rw $path"
      exec sudo chmod a+rw $path
   }
   if { $id == "125c:0020" } {
      set path [string trim "/proc/bus/usb/[lindex $d 1]/[lindex $d 3]" :]
      puts stdout "Found Apogee ASCENT : $d, path $path"
      puts stdout "Invoking sudo chmod to alter rw device permissions..."
      puts stdout "sudo chmod a+rw $path"
      exec sudo chmod a+rw $path
   }
}


