#! /bin/awk -f
#Generic script to loop through the lines in a run summary file
#and match all of the desired parameters
#
#Return 1 if count = total number of parameters checked
#0 if count < total checked (some parameter not found or not matched)
#
#
BEGIN {

#Set parameters here
#
  gain = 9.3;
  readnoise = 7.3;
  recoiltype = "neutron";
  camNum = 0;
  nconditions = 3;
  count = 0;
}
{
  if (($3==gain) && ($2 == camNum) && ($1="Gain")) count++;
  if (($3==readnoise) && ($2 == camNum) && ($1=="ReadNoise")) count++;
  if (($3==recoiltype) && ($2 == camNum) && ($1=="EventType")) count++;
}
END{
#1 for pass 0 for not pass
  if (count == nconditions) print 1;
  else print 0;

}
