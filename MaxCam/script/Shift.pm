package Shift; 


use DM; 

$SALT = "KCN"; 

sub current_shifter_get
{
	$q = $_[1]; 
	$db = DM->connect or die ("could not connect!\n"); 

	$sql = "SELECT shifter FROM shifts ORDER BY idx DESC LIMIT 1"; 
	$stmt = $db->prepare($sql) or die "fail: $sql\n"; 
	$stmt->execute() or die "fail exec $sql\n"; 
	my $row = $stmt->fetchrow_hashref(); 
	my $id = $row->{shifter}; 

	$sql = "SELECT ".$q." FROM shifters WHERE idx = ".$id;
	$stmt = $db->prepare($sql) or die "fail: $sql\n"; 
	$stmt->execute() or die "fail exec $sql\n"; 
	my @row = $stmt->fetchrow_array(); 	
	return $row[0]; 	
}
