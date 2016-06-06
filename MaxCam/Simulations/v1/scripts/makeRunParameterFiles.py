#! /usr/bin/python

def ask_ok(prompt,complaint='Please answer yes or no'):
    """Ask if prompt is ok or not.
    Return True for yes, False for no
    """
    while True:
        ok = raw_input(prompt)
        ok = ok.lower()
        if ok in ('yes','y'): return True
        if ok in ('n','no'): return False
        print(complaint)

def enter_integer(prompt,comment = ''):
    """Makes a prompt for entering an integer value
    Loops until an integer is entered
    """
    number = 0
    done = False
    while not done:
        try:
            inp = raw_input(prompt)
            number = int(inp)
            done = True
            print(comment)
        except ValueError:
            print('Please enter an integer value')
    return number

def enter_decimal(prompt,comment = ''):
    """Makes a prompt for entering a decimal number"""
    number = 0
    done = False
    while not done:
        try:
            inp = raw_input(prompt)
            number = float(inp)
            done = True
            print(comment)
        except ValueError:
            print('Please enter a numerical value')
    return number

def string_switch(prompt,options,comments):
    done = False
    str = ' '
    while not done:
        print(prompt)
        print(list(options))
        str = raw_input("")
        if str in options:
            ind = options.index(str)
            print(options[ind]+": "+comments[ind])
            done = True
        else:
            print('Please enter one of the following options:')
            print(list(options))
    return str


    

        

def run_params():
    print("\n")
    print('We will now review all of the parameters for the run.')
    runN = enter_integer('How many events would you like to generate? ')
    print("\n")
    return runN

def recoil_info():
    print("\nWe will now review the type of recoils to carry out")
    rectype = 'wimp'
    projpart = 'wimp'
    projmass = 100E6
    recpart = 'fluorine'
    recmass = 19E6
    srim = 'SRIM_F_in_CF4_100Torr'
    loc = 'mit'
    done = False
    while not done:
        print("\n")
        print('==============================================')
        print('1. Recoil Type: '+rectype)
        print('2. Projectile Particle: '+str(projpart))
        print('3. Projectile Mass [keV]: '+str(projmass))
        print('4. Recoil Particle: '+str(recpart))
        print('5. Recoil Mass [keV]: '+str(recmass))
        print('6. SRIM table: '+srim)
	print('7. Location: ' + loc)
        print('==============================================')
        print("\n")
        ok = ask_ok('Would you like to change any of these? ')
        if ok:
            param = raw_input('Please enter the number of the parameter you would like to change. ')
            if param == '1':
                prompt = 'Please enter the type of recoil you would like to generate'
                opt = ['wimp','alpha','neutron']
                comment = ['','','']
                rectype = string_switch(prompt,opt,comment)
                if rectype == 'wimp':
                    projpart = 'wimp'
                elif rectype == 'neutron':
                    projpart = 'neutron'
                    projmass = 939566
                elif rectype == 'alpha':
                    projmass = 4001506
                    projpart = 'alpha'
                    recmass = projmass
                    recpart = projpart
            elif param == '2':
                prompt = 'Please enter the name of the projectile you would like to use.'
                opt = ['wimp','alpha','neutron']
                comment = ["Don't forget to set the mass",'','']
                projpart = string_switch(prompt,opt,comment)
                if projpart == 'alpha':
                    projmass = 4001506
                    recpart = 'alpha'
                    recmass = projmass
                elif projpart == 'neutron': projmass = 939566
            elif param == '3':
                prompt = 'Please enter the mass [keV] of the projectile particle. '
                projmass = enter_integer(prompt)
            elif param == '4':
                recpart = input('Please enter the name of the recoil particle you would like to use.')
                print("Don't forget to set the mass.")
            elif param == '5':
                prompt = 'Please enter the mass [keV] of the recoil particle. ' 
                recmass = enter_integer(prompt)
            elif param == '6':
                srim = raw_input('Please enter the name (w/out directories) of the SRIM file you would like to use. ')
            elif param == '7':
                prompt = 'Please enter the location from the list below'
                opt = ['mit','wipp','dusel','soudan','gransasso','snolab','kamioka','boulby','fermilab']
                comment = ['','','','','','','','','']
                loc = string_switch(prompt,opt,comment)
        else: done = True
    rec_stuff = [rectype,projpart,projmass,recpart,recmass,srim,loc]
    return rec_stuff

def energy_info():
    print("\nWe will now review the energy options")
    enopt = 'random'
    minEn = 0.0 
    maxEn = 2.0E4
    fixEn = 0.0
    done = False
    while not done:
        print("\n")
        print('==============================================')
        print('1. Energy option: '+enopt)
        print('2. Minimum Energy [keV]: '+str(minEn))
        print('3. Maximum Energy [keV]: '+str(maxEn))
        print('4. Energy for fixed-energy runs [keV]: '+str(fixEn))
        print('==============================================')
        print("\n")
        ok = ask_ok('Would you like to change any of these parameters? ')
        if ok:
            param = raw_input('Please enter the number of the parameter that you would like to change. ')
            if param == '1':
                prompt = 'Enter energy option: '
                opt = ['random','fixed']
                comment = ["Don't forget to set energy limits","Don't forget to set the energy"]
                enopt = string_switch(prompt,opt,comment)
            elif param == '2': minEn = enter_decimal('Enter the minimum recoil energy [keV]: ')
            elif param == '3': maxEn = enter_decimal('Enter the maximum recoil energy [keV]: ')
            elif param == '4': fixEn = enter_decimal('Enter the recoil energy [keV]: ')
        else: done = True
    en_stuff = [enopt,minEn,maxEn,fixEn]
    return en_stuff

def pos_info():
    print("\nWe will now review the recoil position options")
    posopt = 'random'
    minx = -512*0.143
    maxx = 512 * 0.143
    miny = -512*0.143
    maxy = 512 * 0.143
    minz = 0 
    maxz = 200
    x = 0
    y = 0
    z = 0
    done = False
    while not done:
        print ("\n")
        print ('==============================================')
        print ('1. Position Option: '+str(posopt))
        print ('2. x-Limits [mm]: Min: '+str(minx)+' Max: '+str(maxx))
        print ('3. y-Limits [mm]: Min: '+str(miny)+' Max: '+str(maxy))
        print ('4. z-Limits [mm]: Min: '+str(minz)+' Max: '+str(maxz))
        print ('5. Fixed position [mm]: ('+str(x)+', '+str(y)+', '+str(z)+')')
        print ('==============================================')
        print ("\n")
        ok = ask_ok('Would you like to change any of these parameters? ')
        if ok:
            param = raw_input('Please enter the number of the parameter that you would like to change. ')
            if param == '1':
                prompt = 'Enter position option; '
                opt = ['random','fixed']
                comment = ["Don't forget to set the limits","Don't forget to set the position"]
                posopt = string_switch(prompt,opt,comment)
            elif param == '2':
                print('Enter the x-coordinate limits [mm]; ')
                minx = enter_decimal('Min: ')
                maxx = enter_decimal('Max: ')
            elif param == '3':
                print('Enter the y-coordinate limits [mm]: ')
                miny = enter_decimal('Min: ')
                maxy = enter_decimal('Max: ')
            elif param == '4':
                print('Enter the z-coordinate limits [mm]: ')
                minz = enter_decimal('Min: ')
                maxz = enter_decimal('Max: ')
            elif param == '5':
                print('Enter the coordinates [mm] for a fixed position run: ')
                x = enter_decimal('x: ')
                y = enter_decimal('y: ')
                z = enter_decimal('z: ')
        else: done = True
    pos_stuff = [posopt,[minx,miny,minz],[maxx,maxy,maxz],[x,y,z]]
    return pos_stuff

def time_info():
    print("\nWe will now review the recoil time options")
    timopt = 'random'
    mintime = [2009,1,1,0,0,0]
    maxtime = [2009,12,31,23,59,59]
    fixtime = [2009,6,1,12,0,0]
    tstep = 1
    done = False
    while not done:
        print("\n")
        print('==============================================')
        print('1. Time Option: '+timopt)
        print('2. Minimum time: '+str(mintime[3])+':'+str(mintime[4])+':'+str(mintime[5])+' '+str(mintime[2])+'/'+str(mintime[1])+'/'+str(mintime[0]))
        print('3. Maximum time: '+str(maxtime[3])+':'+str(maxtime[4])+':'+str(maxtime[5])+' '+str(maxtime[2])+'/'+str(maxtime[1])+'/'+str(maxtime[0]))
        print('4. Fixed time: '+str(fixtime[3])+':'+str(fixtime[4])+':'+str(fixtime[5])+' '+ str(fixtime[2])+'/'+str(fixtime[1])+'/'+str(fixtime[0]))
        print('5. Time Step [s]: '+str(tstep))
        print('==============================================')
        print("\n")
        ok = ask_ok('Would you like to change any of these parameters? ')
        if ok:
            param = raw_input('Please enter the number of the parameter that you would like to change. ')
            if param == '1':
                prompt = 'Enter time option: '
                opt = ['random','fixed','current','series']
                comment = ["Don't forget to set limits","Don't forget to set the time","Using current time","Don't forget to set time step"]
                timopt = string_switch(prompt,opt,comment)
            elif param == '2':
                print('Enter the minimum (earliest) time at which event will be generated: ')
                mintime[0] = enter_integer('Year: ')
                mintime[1] = enter_integer('Month: ')
                mintime[2] = enter_integer('Day: ')
                mintime[3] = enter_integer('Hour: ')
                mintime[4] = enter_integer('Minute: ')
                mintime[5] = enter_integer('Second: ')
            elif param == '3':
                print('Enter the maximum (latest) time at which event will be generated: ')
                maxtime[0] = enter_integer('Year: ')
                maxtime[1] = enter_integer('Month: ')
                maxtime[2] = enter_integer('Day: ')
                maxtime[3] = enter_integer('Hour: ')
                maxtime[4] = enter_integer('Minute: ')
                maxtime[5] = enter_integer('Second: ')
            elif param == '4':
                print('Enter the time at which event will be generated: ')
                fixtime[0] = enter_integer('Year: ')
                fixtime[1] = enter_integer('Month: ')
                fixtime[2] = enter_integer('Day; ')
                fixtime[3] = enter_integer('Hour: ')
                fixtime[4] = enter_integer('Minute: ')
                fixtime[5] = enter_integer('Second: ')
            elif param == '5':
                tstep = enter_decimal('Enter the time [sec]  between successive events: ')
        else: done = True
    time_stuff = [timopt,mintime,maxtime,fixtime,tstep]
    return time_stuff

def dir_info():
    print("\nWe will now review the projectile direction options")
    diropt = 'isotropic'
    sourcepos = [0.0,0.0,0.0]
    sourcedir = [1.0,0.0,0.0]
    thetamax = 5
    done = False
    while not done:
        print ("\n")
        print ('==============================================')
        print ('1. Direction Option: ' +diropt)
        print ('2. Source position: ('+str(sourcepos[0])+', '+str(sourcepos[1])+', '+str(sourcepos[2])+')')
        print ('3. Source direction: ('+str(sourcedir[0])+', '+str(sourcedir[1])+', '+str(sourcedir[2])+')')
        print ('4. Max. collimated angle: '+str(thetamax))
        print ('==============================================')
        print ("\n")
        ok = ask_ok('Would you like to change any of these parameters? ')
        if ok:
            param = raw_input('Please enter the number of the parameter that you would like to change. ')
            if param == '1':
                prompt = 'Enter energy option: '
                opt = ['isotropic','fixed','source','collimated']
                comment = ['',"Don't forget to set the direction","Don't forget to set the position",
                           "Don't forget to set the position, direction and maximum angular extent of the beam"]
                diropt = string_switch(prompt,opt,comment)
            elif param == '2':
                print('Enter the source position [mm]')
                sourcepos[0] = enter_decimal('x: ')
                sourcepos[1] = enter_decimal('y: ')
                sourcepos[2] = enter_decimal('z: ')
            elif param == '3':
                print('Enter the source direction')
                print('(You may ignore normalization) ')
                sourcedir[0] = enter_decimal('x: ')
                sourcedir[1] = enter_decimal('y: ')
                sourcedir[2] = enter_decimal('z: ')
            elif param == '4':
                thetamax = enter_decimal ('Enter the maximum angle [deg] of a collimated source from its given direction: ')
        else: done = True
    dir_stuff = [diropt,sourcepos,sourcedir,thetamax]
    return dir_stuff

def spec_info():
    print("\nWe will now review some other options")
    print('If you would like to run with a WIMP distribution or using ENDF neutron files, that is entered here')
    specopt = 'none'
    spectrum = 'ENDF_Cf-252_n_spectrum'
    totalxc = 'ENDF_CS_n_on_19F'
    elastScat = 'ENDF_DCS_n_on_19F'
    timeres = 15
    wimpspectrum='infinity'
    done = False
    while not done:
        print ("\n")
        print ('==============================================')
        print ('1. Special option: '+specopt)
        print ('2. Neutron fission spectrum ENDF: '+spectrum)
        print ('3. Total Cross-Section (CS) ENDF: '+totalxc)
        print ('4. Elastic Scattering (DCS) ENDF: '+elastScat)
        print ('5. Electronic timing resolution [ns]: '+str(timeres))
        print ('6. Theoretical Wimp Energy Spectrum: '+wimpspectrum)
        print ('==============================================')
        print ("\n")
        ok = ask_ok('Would you like to change any of these parameters? ')
        if ok:
            param = raw_input('Please enter the number of the parameter that you would like to change. ')
            if param == '1':
                prompt = 'Enter the desired special option: '
                opt = ['wimp','endf','none','theory']
                comment = ['',"Don't forget to set ENDF file names",'','']
                specopt = string_switch(prompt,opt,comment)
            elif param == '2':
                spectrum = raw_input('Please enter the neutron fission spectrum ENDF file: ')
            elif param == '3':
                totalxc = raw_input('Please enter the total cross section (CS) ENDF file: ')
            elif param == '4':
                elastScat = raw_input('Please enter the elastic scattering (DCS) ENDF file: ')
            elif param == '6':
                prompt = 'Enter the desired wimp energy spectrum:'
                opt = ['flat','vesc','infinity']
                comment=['','','']
                wimpspectrum=string_switch(prompt,opt,comment)
        else: done = True
    spec_stuff = [specopt,spectrum,totalxc,elastScat,timeres,wimpspectrum]
    return spec_stuff

def cam_info(camNum,filename):
    print("\nWe will now review the camera properties")
    campos = [0,0]
    cambins = [256,256]
    imagewidths = [0.143*1024,0.143*1024]
    pixperbin = 4
    readnoise = 7.3
    bias = 0
    emgain = 1
    noisefact = 1
    darkcurrent = 0
    camgain = 9.3
    done = False
    while not done:
        print ("\n")
        print ('==============================================')
        print ('1. Camera Position [mm]: ('+str(campos[0])+ ', ' +str(campos[1])+ ')')
        print ('2. Camera Bins: ('+str(cambins[0])+', '+str(cambins[1])+')')
        print ('3. Image widths [mm]: ('+str(imagewidths[0])+', '+str(imagewidths[1])+')')
        print ('4. Pixels per bin: '+str(pixperbin))
        print ('5. Read noise [count/bin]: '+str(readnoise))
        print ('6. Bias [count/bin]: '+str(bias))
        print ('7. Gain: '+str(camgain))
        print ('8. Noise factor from EMCCD: '+str(noisefact))
        print ('9. Dark current [count/bin]: '+str(darkcurrent))
        print ('10.EM Gain: '+str(emgain))
        print ('==============================================')
        print ("\n")
        ok = ask_ok('Would you like to change any of these parameters? ')
        if ok:
            param = raw_input('Please enter the number of the parameter that you would like to change. ')
            if param == '1':
                print('Enter the camera x-y position [mm]: ')
                campos[0] = enter_decimal('x: ')
                campos[1] = enter_decimal('y: ')
            elif param == '2':
                print('Enter the number of bins on each side: ')
                cambins[0] = enter_integer('x: ')
                cambins[1] = enter_integer('y: ')
            elif param == '3':
                print('Enter the image width [mm] on each side: ')
                imagewidths[0] = enter_decimal('x: ')
                imagewidths[1] = enter_decimal('y: ')
            elif param == '4':
                pixperbin = enter_integer('Enter the width of each bin in pixels: ')
            elif param == '5':
                readnoise = enter_decimal('Enter the ccd read noise: ')
            elif param == '6':
                bias = enter_decimal('Enter the ccd bias: ')
            elif param == '7':
                camgain = enter_decimal('Enter the camera gain: ')
            elif param == '8':
                noisefact = enter_decimal('Enter the EMCCD noise factor: ')
            elif param == '9':
                darkcurrent = enter_decimal('Enter the dark current [counts]: ')
            elif param == '10':
                emgain = enter_decimal('Enter the EM gain: ')
        else: done = True
    cam_stuff = [campos,cambins,imagewidths,pixperbin,readnoise,bias,camgain,noisefact,
                 darkcurrent,emgain]
    print("Saving Camera " + str(camNum))
    file = open(filename,"w")
    file.write("#Properties of Camera " + str(camNum))
    file.write("\nCameraPosition\t\t"+str(camNum)+"\t"+str(cam_stuff[0][0])+"\t"+str(cam_stuff[0][1]))
    file.write("\nCameraBins\t\t"+str(camNum)+"\t"+str(cam_stuff[1][0])+"\t"+str(cam_stuff[1][1]))
    file.write("\nImageWidths\t\t"+str(camNum)+"\t"+str(cam_stuff[2][0])+"\t"+str(cam_stuff[2][1]))
    file.write("\nPixelsPerBin\t\t"+str(camNum)+"\t"+str(cam_stuff[3]))
    file.write("\nReadNoise\t\t"+str(camNum)+"\t"+str(cam_stuff[4]))
    file.write("\nBias\t\t\t"+str(camNum)+"\t"+str(cam_stuff[5]))
    file.write("\nGain\t\t\t"+str(camNum)+"\t"+str(cam_stuff[6]))
    file.write("\nNoiseFactor\t\t"+str(camNum)+"\t"+str(cam_stuff[7]))
    file.write("\nDarkCurrent\t\t"+str(camNum)+"\t"+str(cam_stuff[8]))
    file.write("\nEMGain\t\t\t"+str(camNum)+"\t"+str(cam_stuff[9]))
    file.close()
    return ask_ok('Would you like to enter another camera? ')

def chamber_info1():
    print("\nWe will now review the first set of chamber parameters.")
    spacdiam = 2.5
    spacspac = 20
    spacaxis = 'y'
    height = 508
    driftlength = 200
    pressure = 75
    temperature = 300
    driftvoltage = 5000
    anodevoltage = 730
    topbottom = '1'
    orangle = 0
    done = False
    while not done:
        print ("\n")
        print ('==============================================')
        print ('1. Spacer diameter [mm]: '+str(spacdiam))
        print ('2. Spacer spacing [mm]: ' + str(spacspac))
        print ('3. Spacer Axis: ' + spacaxis)
        print ('4. Height [mm]: ' + str(height))
        print ('5. Drift length [mm]: ' + str(driftlength))
        print ('6. Pressure [torr]: ' + str(pressure))
        print ('7. Temperature [K]: ' + str(temperature))
        print ('8. Drift voltage [V]: ' + str(driftvoltage))
        print ('9. Anode voltage [V]: ' + str(anodevoltage))
        print ('10.Orientation angle [deg]: ' +str(orangle))
	print ('11.Top (1) or bottom (0) TPC: ' + topbottom)
        print ('==============================================')
        print ("\n")
        ok = ask_ok('Would you like to change any of these parameters? ')
        if ok:
            param = raw_input('Please enter the number of the parameter that you would like to change. ')
            if param == '1':
                spacdiam = enter_decimal('Enter the spacer diameter [mm]: ')
            elif param == '2':
                spacspac = enter_decimal('Enter the distance between spacers [mm]: ')
            elif param == '3':
                prompt = 'Enter the axis along which spacers are oriented'
                opt = ['x','y']
                comment = ['','']
                spacaxis = string_switch(prompt,opt,comment)
            elif param == '4':
                height = enter_decimal('Enter the height [mm] (distance between mesh and lens): ')
            elif param == '5':
                driftlength = enter_decimal('Enter the drift lenght [mm]: ')
            elif param == '6':
                pressure = enter_decimal('Enter the chamber pressure [torr]: ')
            elif param == '7':
                temperature = enter_decimal('Enter the chamber temperature [K]: ')
            elif param == '8':
                driftvoltage = enter_decimal('Enter the drift voltage [V]: ')
            elif param == '9':
                anodevoltage = enter_decimal('Enter the anode voltage [V]: ')
	    elif param == '10':
		orangle = enter_decimal('Enter the angle of the +y axis in deg. west (counterclockwise) of north: ')
            elif param == '11':
                prompt = 'Enter 1 for top TPC, 0 for bottom TPC'
                opt = ['1','0']
                comment = ['','']
                topbottom = string_switch(prompt,opt,comment)
        else: done = True
    chamb_stuff1 = [spacdiam,spacspac,spacaxis,height,driftlength,pressure,temperature,
                    driftvoltage,anodevoltage,orangle,topbottom]
    return chamb_stuff1

def chamber_info2():
    print("\nNow we will review the second set of chamber parameters")
    diffconst = 0.25
    diffdz = 0.00419
    usediffdz = 0
    elecscint = 1.0
    nuclscint = 0.3
    atten = 0
    elecperkev = 10000
    driftvel = 0.11
    sigmavel = 0.015
    ldiffconst = 400
    ldiffdz = 3
    done = False
    while not done:
        print ("\n")
        print ('==============================================')        
        print ('1. Diffusion constant term [mm^2]: ' + str(diffconst))
        print ('2. Diffusion dz term [mm^2/mm drift]: ' + str(diffdz))
        print ('3. Use diffusion dz term? (Override other related params): '+str(usediffdz) + ' Note: you should leave this off')
        print ('4. Electric scintillation: ' + str(elecscint))
        print ('5. Nuclear scintillation: ' +str(nuclscint))
        print ('6. Attenuation [mm^-1]: '+str(atten))
        print ('7. Electrons per keV (at mesh): ' + str(elecperkev))
        print ('8. Drift Velocity [mm/ns]: ' + str(driftvel))
        print ('9. Longitudinal diffusion constant term [ns^2]: ' + str(ldiffconst))
        print ('10.Longitudinal diffusion dz term [ns^2/mm]: ' + str(ldiffdz))
        print ('==============================================')
        print ("\n")
        ok = ask_ok('Would you like to change any of these parameters? ')
        if ok:
            param = raw_input('Please enter the number of the parameter that you would like to change. ')
            if param == '1':
                diffconst = enter_decimal('Enter the diffusion constant term [mm^2]: ')
            elif param == '2':
                diffdz = enter_decimal('Enter the diffusion dz term [mm^2 / mm drift]')
            elif param == '3':
                prompt = 'Enter whether or not to override the other parameters and just use diffusion dz term as the correct diffusion dz term'
                opt = ['true','false']
                comment = ['','']
                yesno = string_switch(prompt,opt,comment)
                if yesno == 'true': usediffdz = 1
                else: usediffdz = 0
            elif param == '4':
                elecscint = enter_decimal('Enter the fraction of electric scintillation to total electric energy loss: ')
            elif param == '5':
                nuclscint = enter_decimal('Enter the fraction of nuclear scintillation to total nuclear energy loss: ')
            elif param == '6':
                attenuation = enter_decimal('Enter the exponential constant [>0, in mm^-1] governing the exponential fall of drift electrons hitting the mesh: ')
            elif param == '7':
                elecperkev = enter_decimal('Enter the electrons per keV hitting the electronic readout: ')
            elif param == '8':
                driftvel = enter_decimal('Enter the drift velocity in mm/ns: ')
            elif param == '9':
                ldiffconst = enter_decimal('Enter the longitudinal diffusion constant term [ns^2]: ')
            elif param == '10':
                ldiffdz = enter_decimal('Enter the longitudinal diffusion dz term [ns^2/mm]: ')
        else: done = True
    chamb_stuff2 = [diffconst,diffdz,usediffdz,elecscint,nuclscint,atten,elecperkev,driftvel,ldiffconst,ldiffdz]
    return chamb_stuff2

    
if __name__ == "__main__":
    print('Hello!')
    print('If you wish to run using a preexisting text file, please exit this program and execute the appropriate script.')
    print('We will now review all of the parameters necessary to begin this run')
    ok = ask_ok('Would you like to continue? ')
    if ok:
        Nevents = run_params()
        recoil_stuff = recoil_info()
        energy_stuff = energy_info()
        pos_stuff = pos_info()
        time_stuff = time_info()
        dir_stuff = dir_info()
        spec_stuff = spec_info()
        done = False
        camNum = 0
        path = "runParameters/"
        camfilename = "cameraProperties"
        while not done:
            filename = path+camfilename+str(camNum)+".temp"
            done = not cam_info(camNum,filename)
            camNum+=1
        chamb_stuff1 = chamber_info1()
        chamb_stuff2 = chamber_info2()
        runfilename = 'runProperties.temp'
        print ('All parameters set, creating text file ' + runfilename)
        file = open(path+runfilename,'w')
        file.write('Run Parameters')
        file.write("\n===================================")
        file.write("\nNumberOfEvents\t\t" + str(Nevents))
        file.write("\nEventType\t\t" + str(recoil_stuff[0]))
        file.write("\nProjectileParticle\t"+ recoil_stuff[1])
        file.write("\nProjectileMass\t\t" + str(recoil_stuff[2]))
        file.write("\nRecoilParticle\t\t" + recoil_stuff[3])
        file.write("\nRecoilMass\t\t" + str(recoil_stuff[4]))
        file.write("\nSRIMfile\t\t" + recoil_stuff[5])
        file.write("\nLocation\t\t" + recoil_stuff[6])
        file.write("\nEnergyOption\t\t"+energy_stuff[0])
        file.write("\nMinEnergy\t\t"+str(energy_stuff[1]))
        file.write("\nMaxEnergy\t\t"+str(energy_stuff[2]))
        file.write("\nFixEnergy\t\t"+str(energy_stuff[3]))
        file.write("\nPositionOption\t\t"+pos_stuff[0])
        file.write("\nMinPositionLimits\t"+str(pos_stuff[1][0])+"\t"+str(pos_stuff[1][1])+"\t" +str(pos_stuff[1][2]))
        file.write("\nMaxPositionLimits\t"+str(pos_stuff[2][0])+"\t"+str(pos_stuff[2][1])+"\t" +str(pos_stuff[2][2]))
        file.write("\nFixPosition\t\t"+str(pos_stuff[3][0])+"\t"+str(pos_stuff[3][1])+"\t" +str(pos_stuff[3][2]))
        file.write("\nTimeOption\t\t"+time_stuff[0])
        file.write("\nMinTime\t\t\t"+str(time_stuff[1][0])+"\t"+str(time_stuff[1][1])+"\t"+str(time_stuff[1][2])+"\t"+str(time_stuff[1][3])+"\t"+str(time_stuff[1][4])+"\t"+str(time_stuff[1][5]))
        file.write("\nMaxTime\t\t\t"+str(time_stuff[2][0])+"\t"+str(time_stuff[2][1])+"\t"+str(time_stuff[2][2])+"\t"+str(time_stuff[2][3])+"\t"+str(time_stuff[2][4])+"\t"+str(time_stuff[2][5]))
        file.write("\nFixTime\t\t\t"+str(time_stuff[3][0])+"\t"+str(time_stuff[3][1])+"\t"+str(time_stuff[3][2])+"\t"+str(time_stuff[3][3])+"\t"+str(time_stuff[3][4])+"\t"+str(time_stuff[3][5]))
        file.write("\nTimeStep\t\t"+str(time_stuff[4]))
        file.write("\nDirectionOption\t\t"+dir_stuff[0])
        file.write("\nSourcePosition\t\t"+str(dir_stuff[1][0])+"\t"+str(dir_stuff[1][1])+"\t"+str(dir_stuff[1][2]))
        file.write("\nSourceDirection\t\t"+str(dir_stuff[2][0])+"\t"+str(dir_stuff[2][1])+"\t"+str(dir_stuff[2][2]))
        file.write("\nThetaMax\t\t"+str(dir_stuff[3]))
        file.write("\nSpecialOption\t\t"+spec_stuff[0])
        file.write("\nFissionSpectrumENDF\t"+spec_stuff[1])
        file.write("\nTotalCS_ENDF\t\t"+spec_stuff[2])
        file.write("\nElastScatterDCS_ENDF\t"+spec_stuff[3])
        file.write("\nTimeResolution\t\t"+str(spec_stuff[4]))
        file.write("\nWimpEnergyOption\t"+spec_stuff[5])
        file.close()
        chamfilename = "chamberProperties.temp"
        file2 = open(path+chamfilename,"w")
        file2.write("SpacerDiameter\t\t"+str(chamb_stuff1[0]))
        file2.write("\nSpacerSpacing\t\t"+str(chamb_stuff1[1]))
        file2.write("\nSpacerAxis\t\t"+chamb_stuff1[2])
        file2.write("\nHeight\t\t\t" + str(chamb_stuff1[3]))
        file2.write("\nDriftLength\t\t"+str(chamb_stuff1[4]))
        file2.write("\nPressure\t\t"+str(chamb_stuff1[5]))
        file2.write("\nTemperature\t\t"+str(chamb_stuff1[6]))
        file2.write("\nDriftVoltage\t\t"+str(chamb_stuff1[7]))
        file2.write("\nAnodeVoltage\t\t"+str(chamb_stuff1[8]))
        file2.write("\nOrientationAngle\t"+str(chamb_stuff1[9]))
        file2.write("\nTopOrBottomTPC\t\t"+chamb_stuff1[10])
        file2.write("\nDiffusionConstantTerm\t"+str(chamb_stuff2[0]))
        file2.write("\nDiffusionDzTerm\t\t"+str(chamb_stuff2[1])+"\t"+str(chamb_stuff2[2]))
        file2.write("\nElectricScintillation\t"+str(chamb_stuff2[3]))
        file2.write("\nNuclearScintillation\t"+str(chamb_stuff2[4]))
        file2.write("\nAttenuation\t\t"+str(chamb_stuff2[5]))
        file2.write("\nElectronPerkeV\t\t"+str(chamb_stuff2[6]))
        file2.write("\nDriftVelocity\t\t"+str(chamb_stuff2[7]))
        file2.write("\nLongDiffusionConstTerm\t"+str(chamb_stuff2[8]))
        file2.write("\nLongDiffusionDzTerm\t"+str(chamb_stuff2[9]))
        
        file2.close()
        filename = "fileManager.temp"
        file3 = open(path+filename,"w")
        file3.write(runfilename)
        file3.write("\n"+chamfilename)
        for a in range(0,camNum):
            file3.write("\n"+camfilename + str(a) +".temp")
        file3.close()
