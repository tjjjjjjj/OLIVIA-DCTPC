#!/usr/bin/perl

use DM;

use CGI ':standard';
use CGI::Carp qw(fatalsToBrowser);

print "Content-type:text/html\n\n";

#if (!param("table")) {
#  die "Unknown table: $table \n";
#}

#my $table = param("table");
my $table="pressure";

#print $table\n;


$db_handle = DM->connect
    or die("Could not connect! \n");


$sql = "select * from $table ORDER BY timestamp DESC LIMIT 1";
$statement = $db_handle->prepare($sql)
    or die "Couldn't prepare query '$sql': $DBI::errstr\n";

$statement->execute()
    or die "Couldn't execute query '$sql': $DBI::errstr\n";


while ($row_ref = $statement->fetchrow_hashref())
{
    print  $row_ref->{value};    
}

$db_handle->disconnect();


