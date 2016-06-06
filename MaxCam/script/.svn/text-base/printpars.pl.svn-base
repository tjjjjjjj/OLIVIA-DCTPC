#!/usr/bin/perl

use DM;

use CGI ':standard';
use CGI::Carp qw(fatalsToBrowser);

print "Content-type:text/html\n\n";

#if (!param("table")) {
#  die "Unknown table: $table \n";
#}

$db_handle = DM->connect
    or die("Could not connect! \n");

print '<h3> Status Reports </h3>';
print '<small>';
print '<table border="1">'."\n";
my $datetime = localtime();
my $utctime = gmtime();

print "<tr><td> Item </td><td> Status</td></tr>";
print "<tr><td>Local Time </td><td> $datetime </td></tr> \n\n";
print "<tr><td>UTC Time </td><td> $utctime </td></tr>\n\n";

$slowControlStatus="<span style='background-color:red'>OFF</span>";
@ARGV = ("ps -A | grep DMSlow |");
while (<>) {
  $slowControlStatus="<span style='background-color:green'>ON</span>";
}
print "<tr><td>Slow control</td><td>$slowControlStatus</td></tr>\n";

$sql = "select * from hvstatus ORDER BY timestamp DESC LIMIT 1";
$statement = $db_handle->prepare($sql)
    or die "Couldn't prepare query '$sql': $DBI::errstr\n";
$statement->execute()
    or die "Couldn't execute query '$sql': $DBI::errstr\n";
$row_ref = $statement->fetchrow_hashref();
if ($row_ref->{setval} > 0) {
	print "<tr><td>HV</td><td> <span style='background-color:green'>ON</span></td></tr>\n";
}
else {
	print "<tr><td>HV </td><td><span style='background-color:red'>OFF</span></td></tr>\n";
}
print "</table>\n\n";
print "<h3>Chamber Data<\h3><table border='1'>";
print "<tr><td> Item </td><td> Value</td><td> Set Value </td><td> Status</td></tr>";

#Get Busy status:

$sql = "select * from busy"; 
$statement = $db_handle->prepare($sql)
    or die "Couldn't prepare query '$sql': $DBI::errstr\n";
$statement->execute()
    or die "Couldn't execute query '$sql': $DBI::errstr\n";
$row_ref = $statement->fetchrow_hashref(); 

$meshbusy = $row_ref->{mesh_hv};
$wirebusy = $row_ref->{wire_hv};
$pressurebusy = $row_ref->{pressure}; 

$meshstatus = "RDY";
$wirestatus = "RDY";
$pressurestatus = "RDY";


if ($meshbusy) {$meshstatus = "<b>BSY</b>"; }
if ($wirebusy) {$wirestatus = "<b>BSY</b>"; }
if ($pressurebusy) {$pressurestatus = "<b>BSY</b>"; }

$sql = "select * from pressure ORDER BY timestamp DESC LIMIT 1";
$statement = $db_handle->prepare($sql)
    or die "Couldn't prepare query '$sql': $DBI::errstr\n";
$statement->execute()
    or die "Couldn't execute query '$sql': $DBI::errstr\n";
$row_ref = $statement->fetchrow_hashref();
print  "<tr><td>Pressure </td> <td>$row_ref->{value}</td><td>$row_ref->{setval}</td><td>$pressurestatus</td></tr> \n";    


$sql = "select * from temp0 ORDER BY timestamp DESC LIMIT 1";
$statement = $db_handle->prepare($sql)
    or die "Couldn't prepare query '$sql': $DBI::errstr\n";
$statement->execute()
    or die "Couldn't execute query '$sql': $DBI::errstr\n";
$row_ref = $statement->fetchrow_hashref();
print  "<tr><td>Temperature</td><td> $row_ref->{value}</td><td>$row_ref->{setval}</td><td>TBD</td></tr> \n";    


$sql = "select * from wire_hv ORDER BY timestamp DESC LIMIT 1";
$statement = $db_handle->prepare($sql)
    or die "Couldn't prepare query '$sql': $DBI::errstr\n";
$statement->execute()
    or die "Couldn't execute query '$sql': $DBI::errstr\n";
$row_ref = $statement->fetchrow_hashref();
print  "<tr><td>Anode Voltage</td><td> $row_ref->{value}</td><td>$row_ref->{setval}</td><td>$wirestatus</td></tr> \n";    


$sql = "select * from mesh_hv ORDER BY timestamp DESC LIMIT 1";
$statement = $db_handle->prepare($sql)
    or die "Couldn't prepare query '$sql': $DBI::errstr\n";
$statement->execute()
    or die "Couldn't execute query '$sql': $DBI::errstr\n";
$row_ref = $statement->fetchrow_hashref();
print  "<tr><td>Drift Voltage</td><td>$row_ref->{value}</td><td> $row_ref->{setval}</td><td>$meshstatus</td> \n";    
print "</table>";
print "</small>";
print end_html();
$statement->finish();
$db_handle->disconnect();


