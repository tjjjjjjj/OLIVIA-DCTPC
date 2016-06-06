package DM;



use DBI;

#use CGI qw(:standard *table);
#use CGI::Carp qw(fatalsToBrowser);
use Net::SSH::Perl; 

$HOST = "mitdm004.mit.edu";
$DETECTOR = "10L";
$SYNACCESS = "18.146.0.100";
@CAMERA_HOST = ("mitdm00.mit.edu");
@CAMERA_NAME = ("cam0");
$DAQ = "mitdm03.mit.edu"; 
$HTMLDIR = "html/";
$APACHE_ROOT = "/var/www/"; 
$DAQ_IMG_DIR = "/home/dmatter/images";
$DAQ_WORK_DIR = "/home/dmatter/work";

$ADMIN= 'cozzyd@mit.edu';

sub ssh_connect
{
  my $ssh = Net::SSH::Perl->new(@_[1],options=>["protocol 2,1"]);   
  $ssh->login("dmatter","seedark");
  return $ssh; 
}

sub connect {

    my $database = "DM_SLOWCONTROL";
    my $hostname = $HOST; 
    my $username = "dmatter";
    my $password = "seedark";

    my $dsn = "DBI:mysql:database=$database;host=$hostname";

    my $db_handle = DBI->connect($dsn, $username, $password) 
        or die("Could not connect! \n");

    return $db_handle;
}

sub check_hold
{
  $db = DM->connect(); 
  $sql = "SELECT slow_hold FROM busy"; 
  $statement = $db->prepare($sql) or die "can't prepare query '$sql': $DBI::errstr\n"; 
  $statement->execute() or die "can't execute query '$sql': $DBI::errstr\n"; 
  $row_ref= $statement->fetchrow_hashref(); 
  return $row_ref->{slow_hold}; 
}


1;

