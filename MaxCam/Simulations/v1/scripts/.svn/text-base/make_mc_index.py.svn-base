#! /usr/bin/env python
#--------------------------------------------------------------------------
# File and Version CVS Information:
# $Id: make_mc_index.py,v 1.4 2011/03/24 02:31:08 jbattat Exp $
#--------------------------------------------------------------------------
import os, stat, grp
import optparse
import sys
import time

# HOW TO RUN?  PUT AN EXAMPLE HERE
def get_field_from_mc_sum_file(sumfilename,keyword,fieldnumber):
    f=open(sumfilename,'r')
    #    triggeroutputfile=open('triggeroutputfile_'+str(runToProcess)+'_'+str(whichRun)+'.dat','w')
    
    a=f.readlines()

    for kk in range(len(a)):
        if a[kk].startswith('#')==False:
            a[kk]=a[kk].replace('\n','')
            asplit=a[kk].split('\t')
            asplit=filter (lambda aitem: aitem!='', asplit)
            if len(asplit)>0:
                if asplit[0]==keyword:
                    return asplit[fieldnumber]
    f.close()
            
# the main function
def main():

    # get a list of .sum files in this directory
    mcsumfileprefix='.sum'
    mcsumfilelist=[]
    for root,dirs,files, in os.walk(os.getcwd()):
        for name in files:
            if mcsumfileprefix in name:
                mcsumfilelist.append(name)
                
        # in order not to recurse
        break

    # sort the sum file list by run-number
    #   get the common suffix of all the .sum files (should just be mcsumfileprefix,
    #   but whatever)
    commonprefix=os.path.commonprefix(mcsumfilelist)
    mcsumfilelist=sorted(mcsumfilelist, key=lambda rc1rdfl : int(rc1rdfl[len(commonprefix):len(rc1rdfl)-len(mcsumfileprefix)]) )
    print mcsumfilelist

    filename = "index.html"
    mcindexfile=open(filename,'w')
    #
    # chmod to give user and group read/write and all write and chgrp to dmtpc
    #
    # S_IXYYY   YYY = GRP (group), USR (user), OTH (other)
    #             X = R (read), W (write), X (execute)
    os.chmod(filename, stat.S_IRGRP | stat.S_IWGRP | stat.S_IWUSR | stat.S_IRUSR | stat.S_IROTH)
    # change the group to dmtpc (but need the numeric group id)
    # chown(filename, uid, gid)  -1 means "do not change"
    os.chown(filename, -1, grp.getgrnam('dmtpc').gr_gid)
    
    mcindexfile.write('<html><head></head><body><table valign="center" align="center" width="90%"><tbody><tr align="center" bgcolor="#cc99cc">')
    mcindexfile.write('<td>File</td> <td>Number of Events</td> <td>Keyword</td> <td>Time Created</td> <td>Description</td></tr>')

    for i in range(len(mcsumfilelist)):
        mcindexfile.write('<tr>')
        mcindexfile.write('<td><a href="'+mcsumfilelist[i]+'"> '+mcsumfilelist[i].replace(mcsumfileprefix,'.root')+'</a></td>')
        mcindexfile.write('<td align="right">'+get_field_from_mc_sum_file(mcsumfilelist[i],'NumberOfEvents',1)+'</td>')
        mcindexfile.write('<td align="right">'+get_field_from_mc_sum_file(mcsumfilelist[i],'EventType',1)+'</td>')
        mcindexfile.write('<td>'+time.ctime(os.path.getctime(mcsumfilelist[i].replace(mcsumfileprefix,'.root')))+'</td>')

        # assumes that the x-y # of pixels is the same, and uses the x values
        ncamerabinsx=float(get_field_from_mc_sum_file(mcsumfilelist[i],'CameraBins',2))
        npixelsperbin=float(get_field_from_mc_sum_file(mcsumfilelist[i],'PixelsPerBin',2))
        mmimagewidthx=float(get_field_from_mc_sum_file(mcsumfilelist[i],'ImageWidths',2))
        mmperpixelx=mmimagewidthx/(npixelsperbin*ncamerabinsx)
        
        mcindexfile.write(('<td>'+
                           'CCD Gain='+
                           get_field_from_mc_sum_file(mcsumfilelist[i],'Gain',2)+
                           '<br>'
                           'CCD Read Noise='+
                           get_field_from_mc_sum_file(mcsumfilelist[i],'ReadNoise',2)+
                           '<br>'+
                           'mm/pixel=%.3f'+
                           '<br>'+
                           'Pressure='+
                           get_field_from_mc_sum_file(mcsumfilelist[i],'Pressure',1)+
                           '<br>'+
                           'Drift Voltage='+
                           get_field_from_mc_sum_file(mcsumfilelist[i],'DriftVoltage',1)+
                           '<br>'+
                           'Anode Voltage='+
                           get_field_from_mc_sum_file(mcsumfilelist[i],'AnodeVoltage',1)+
                           '<br>' )% mmperpixelx )

        if get_field_from_mc_sum_file(mcsumfilelist[i],'EnergyOption',1)=='fixed':
            mcindexfile.write(('KE='+
                               get_field_from_mc_sum_file(mcsumfilelist[i],'FixEnergy',1)+
                               ' keV<br>'))
            
        mcindexfile.write('</td>')
        mcindexfile.write('</tr>')

        mcindexfile.write('<tr valign="top">')
        mcindexfile.write('<td colspan="5" style="background-color: #cc99cc"></td>')
        mcindexfile.write('</tr>')
    
    mcindexfile.write('</tbody></table><table width="80%"></table></body></html>')

if __name__ == '__main__':
    main()
