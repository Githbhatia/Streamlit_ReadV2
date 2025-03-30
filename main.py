import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import zipfile, io
from itertools import islice
import pandas as pd
from RS_function import RS_function
import math


@st.cache_data
def readFileV2c(_f, f_name):
    for line in islice(f, 1, 2):   
        recTime = line[10:].strip()
        # print(recTime)
    for line in islice(f, 6, 7):
        nameCh=f_name[:f_name.find(".V2c")-7] + " " +line[13:37].strip()
    # print(nameCh)
    for line in islice(f, 15, 16):
        headerPoints = int(line[:5].strip())
        headerLines = int(line[line.find("lines")-4:line.find("lines")].strip())
    # print(headerPoints); print(headerLines)
    header = readchunk15(f,headerLines)
    # print(header)
    latitude = header[0]
    longitude = header[1]
    # print(latitude, longitude)
    dt = header[33]
    # print(dt)
    for line in islice(f, 0, 1):
        commentLines = int(line[:5].strip())
    for line in islice(f, commentLines, commentLines+1):
        numofPoints =int(line[:9].strip())
    # print(numofPoints)
    accel = readchunk15(f,numofPoints)
    # print(accel)
    f.close()
    return recTime,latitude,longitude,nameCh,dt,numofPoints,accel
def chunkstring15(string, length):
    return (float(string[0+i:length+i]) for i in range(0, len(string), length))

def readchunk15(f, numofLines):
    x=[]
    for line in islice(f, 0,  numofLines):
        x = x + list(chunkstring15(line[0:len(line)-1], 15))
    #print(x)
    return x

def chunkstring10(string, length):
    return (float(string[0+i:length+i]) for i in range(0, len(string), length))

def readchunk(f, numofLines):
    x=[]
    for line in islice(f, 0,  numofLines):
        x = x + list(chunkstring10(line[0:len(line)-1], 10))
    #print(x)
    return x
def scaleValue(units):
    if units =="cm/sec2":
        return 1/980.665
    else:
        return 1.0
def lines(points):
    if points % 8 == 0:
        nLines = int(points/8) 
    else:
        nLines = int(points/8)+1
    return nLines

def startlimAccel(z):
    a1 = next(i for i, x in enumerate(accel1) if abs(x) >z)
    startTime  = max(a1*dtAccel1- 2, 0)
    if EOF==0:
        a2 = next(i for i, x in enumerate(accel2) if abs(x) >z)
        a3 = next(i for i, x in enumerate(accel3) if abs(x) >z)
        startTime  = max(min(a1*dtAccel1,a2*dtAccel2,a3*dtAccel3) - 2, 0)
    return round(startTime,2)

def endlimAccel(z):
    a1 = next(i for i, x in reversed(list(enumerate(accel1))) if abs(x) >z)
    endTime  = max(a1*dtAccel1+ 2, 0)
    if EOF==0:
        a2 = next(i for i, x in reversed(list(enumerate(accel2))) if abs(x) >z)
        a3 = next(i for i, x in reversed(list(enumerate(accel3))) if abs(x) >z)
        endTime  = max(min(a1*dtAccel1,a2*dtAccel2,a3*dtAccel3) +2, 0)
    return round(endTime,2)

def maxaccel(x, t):
    ymax = max(x)
    xpos = x.index(ymax); xmax = t[xpos]
    return [xmax, ymax]

def minaccel(x, t):
    ymin = min(x)
    xpos = x.index(ymin); xmin = t[xpos]
    return [xmin, ymin]

def saveFile():
    textstring=""
    j=0
    textstring += "Time(sec)"
    if wch1:
        textstring += ", " + nameCh1.replace(",","_").replace(" ", "_")
    if wch2:
        textstring += ", " + nameCh2.replace(",","_").replace(" ", "_")
    if wch3:
        textstring += ", " + nameCh3.replace(",","_").replace(" ", "_")
    textstring += "\n"
    index = len(T1)
    while j < index:
        textstring += str(round(T1[j],3))
        if wch1:
            textstring += ", " + str(scaledAccel1[j]) 
        if wch2:
            textstring += ", " + str(scaledAccel2[j])
        if wch3:
            textstring += ", " + str(scaledAccel3[j])
        textstring += "\n"
        j+= 1
    return (textstring)

def rsaveFile():
    textstring=""
    j=0
    textstring += "Time_Period(sec)"
    i = 0
    while i < len(xi):
        if rch1:
            textstring += ", " + nameCh1.replace(",","_").replace(" ", "_") + "_" + str(round(xi[i],3))
        if rch2:
            textstring += ", " + nameCh2.replace(",","_").replace(" ", "_") + "_" + str(round(xi[i],3))
        if rch3:
            textstring += ", " + nameCh3.replace(",","_").replace(" ", "_") + "_" + str(round(xi[i],3))
        i+=1
    textstring += "\n"
    index = len(tT)
    while j < index:
        textstring += str(round(tT[j],3))
        i = 0
        while i < len(xi):
            if rch1:
                textstring += ", " + str(sAll[0][i][j]) 
            if rch2:
                textstring += ", " + str(sAll[1][i][j])
            if rch3:
                textstring += ", " + str(sAll[2][i][j])
            i+=1
        textstring += "\n"
        j+= 1
    return (textstring)

def accelim(x,y,z):
    xmax = max([abs(i) for i in x])
    ymax = max([abs(i) for i in y])
    zmax = max([abs(i) for i in z])
    return max([xmax,ymax,zmax])

def tripartitegrids(scale,plt,ax,xl,xr):
    aSeries = np.array([0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100])
    dSeries = np.array([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,100,200,300,400,500,600,700,800,900,1000])
    bset= np.array([0.001,0.01,0.1,1,10,100,1000])
    bset2 = np.array([0.006,0.008,0.009,0.06,0.08,0.09,0.6,0.8,0.9,6,7,8,9,60,80,90,600,800,900])
    periodLimit, velLimit =  ax.transData.inverted().transform(ax.transAxes.transform((0.95,0.95)))
    periodLimit0, velLimit0 =  ax.transData.inverted().transform(ax.transAxes.transform((0.0,0.0)))
   

    for i, items in enumerate(aSeries):
        t0 =2*np.pi*velLimit0/(aSeries[i]*scale)
        t1= 2*np.pi*velLimit/(aSeries[i]*scale)
        
        t=t1;v=velLimit
        m = str(aSeries[i])+ "g"
        if aSeries[i] in bset:
            ax.plot([t0,t],[velLimit0,velLimit], linestyle="--", color= 'k',linewidth=0.3)
            if t < 0.95*xr and t > 1.05*xl:
                ax.annotate(m, xy=(t,v), xytext=(t,v), fontsize=5, color= 'k')
        elif aSeries[i] in bset2:
            ax.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'c',linewidth=0.3)
        else:
            ax.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'c',linewidth=0.3)
            if t < 0.95*xr and t > 1.05*xl:
                ax.annotate(m, xy=(t,v), xytext=(t,v), fontsize=5, color= 'c')

    for i, items in enumerate(dSeries):
        t0 =2*np.pi*dSeries[i]/velLimit0
        t1= 2*np.pi*dSeries[i]/velLimit
        t=t0; v=velLimit0
        m= str(dSeries[i])+"cm"
        if dSeries[i] in bset:        
            ax.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'k',linewidth=0.3)
            if t < 0.95*xr and t > 1.05*xl:
                ax.annotate(m, xy=(t,v), xytext=(t,v), ha='left', va="top",fontsize=5, color= 'k')
        elif dSeries[i] in bset2:
            ax.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'm',linewidth=0.3)
        else:
            ax.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'm',linewidth=0.3)
            if t < 0.95*xr and t > 1.05*xl:
                ax.annotate(m, xy=(t,v), xytext=(t,v), ha='left', va="top",fontsize=5, color= 'm')
    return(1)

def resSpectrafn(accel,ax,rU,rT,channel):
    Sfin=[]
    for i in range(0,len(xi)):
        Sfin= RS_function(accel[int(starttime/dtAccel1):int(endtime/dtAccel1)], df, tT, xi[i], Resp_type = rT)
        if option =="Accel":
            sAll[channel,i,:]=Sfin[0,:]*scaleValue(unitsAccel1)
        else:
            sAll[channel,i,:]=Sfin[0,:]
        amax=[tT[np.argmax(abs(sAll[channel,i,:]))], max(abs(sAll[channel,i,:]))]
        labl = "Damping = "+ str(round(xi[i],3))+ ": Max at "+ str(round(amax[0],3)) +"sec, "+str(round(amax[1],2)) + rU
        ax.plot(tT,sAll[channel,i,:],label=labl,linewidth=1.0)
        
    ax.grid()
    ax.set_ylabel(rL)
    return(1)


def resTripSpectrafn(accel,ax):
    Sfin=[]
    for i in range(0,len(xi)):       
        Sfin= RS_function(accel[int(starttime/dtAccel1):int(endtime/dtAccel1)], df, tT, xi[i], Resp_type = 'SA')
        S=Sfin[0,:]/(2*np.pi/tT)
        labl = 'Damping=' + str(round(xi[i],3))
        ax.plot(tT,S,linewidth=1.0, label = labl)
        

    x_left = np.min(tT); x_right =max(tT)
    ax.set_ylabel('Psuedo Velocity '+ unitsVel1)
    ax.grid()
    ax.set_xscale("log")
    ax.set_yscale("log")
    # x,y =ax.get_ylim()
    # y = y*1.2
    # ax.set_ylim(x,y)
    tripartitegrids(1/scaleValue(unitsAccel1),plt,ax,x_left,x_right)
    ax.set_xlim(x_left,x_right)
    return(1)

def orbitplotfn():
    sLoc = int(starttime/dtAccel1); eLoc = int(endtime/dtAccel1)
    ax.plot(orx[sLoc:eLoc], ory[sLoc:eLoc])
    rotmaxLoc = np.argmax(np.sqrt(np.square(orx[:])+np.square(ory[:])))
    resmax = np.sqrt(np.square(orx[rotmaxLoc])+np.square(ory[rotmaxLoc]))
    resAngle = np.arctan2(orx[rotmaxLoc],ory[rotmaxLoc])
    # print(rotmaxLoc,resAccelmax, xa[rotmaxLoc]*np.cos(resAngle)+ya[rotmaxLoc]*np.sin(resAngle) )
    ax.plot([0,orx[rotmaxLoc]], [0, ory[rotmaxLoc]], color='red',linewidth=2.0 )
    ax.annotate(str(round(resmax,3)) + "@ " +str(round(resAngle*180/math.pi,2))+ r"$^\circ$", xy=(orx[rotmaxLoc], ory[rotmaxLoc]), xytext=(orx[rotmaxLoc], ory[rotmaxLoc]), fontsize=10, color= 'Blue')
    ax.set_xlabel(xRec + ' ' + rT); ax.set_ylabel(yRec + ' ' + rT)
    maxLimit = max(np.max(orx[sLoc:eLoc]), np.max(ory[sLoc:eLoc]),np.abs(np.min(orx[sLoc:eLoc])),np.abs(np.min(ory[sLoc:eLoc])))/0.95
    ax.set_xlim(-maxLimit, maxLimit)
    ax.set_ylim(-maxLimit, maxLimit)
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))
    xlabel=ax.get_xticks()
    zind = np.where(xlabel == 0)[0][0]
    for i in range(zind,len(xlabel)):
        cr = plt.Circle((0, 0), xlabel[i], linestyle="--", color= 'k',linewidth=0.3, fill=False)
        ax.add_patch(cr)
    ax.grid()
    return(1)

def adrs(accel,ax):
    Sfin= RS_function(accel[int(float(starttime/dtAccel1)):int(float(endtime)/dtAccel1)], df, tT, dampoption, Resp_type = 'PSASD')
    S=Sfin[0,:]*scaleValue(unitsAccel1)
    area= round(np.trapezoid(Sfin[0,:],Sfin[1,:])/10000,2)
    ax.set_xlabel('Peak D (cm)')

    ax.plot(Sfin[1,:],S,color= 'Red', linewidth=1.0)
    SfinClosed = np.append(np.insert(Sfin[1,:],0,0.0),0.0)
    SClosed = np.append(np.insert(S,0,0.0),0.0)
    ax.fill(SfinClosed,SClosed, "r", alpha=0.5)
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))
    radialPeriods(1/scaleValue(unitsAccel1), ax)
    ax.text(x_right/3, y_high/3, str(area) + r"$(m/s)^2$", horizontalalignment='center', fontsize=10, color ='Blue')
    ax.text(0.97, 0.97, 'Damping=' + str(round(dampoption,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
    return(1)

def radialPeriods(scale, ax):
    periodSeries = np.concatenate(( np.arange(0.1,1.0,0.1) , np.arange(1.0,2.0,0.5), np.arange(2.0,5.0,1) ))
    #print(periodSeries)

    dispLimit, AccelLimit = ax.transData.inverted().transform(ax.transAxes.transform((0.95,0.95)))

    w = 2*np.pi/periodSeries
    w2 = np.square(w)
    a1 = dispLimit*w2/scale
    d1 = a1*scale/w2
    d2 = AccelLimit*scale/w2
    a2 = d2*w2/scale

    for i, items in enumerate(d1):
        if a1[i] < AccelLimit:
            ax.plot([0, d1[i]],[0, a1[i]], linestyle="--", color= 'Blue',linewidth=0.4)
            ax.annotate(round(periodSeries[i],1), xy=(d1[i], a1[i]), xytext=(d1[i], a1[i]), fontsize=5, color= 'Blue')
        else:
            ax.plot([0, d2[i]],[0, a2[i]], linestyle="--", color= 'Blue',linewidth=0.4)
            ax.annotate(round(periodSeries[i],1), xy=(d2[i], a2[i]), xytext=(d2[i], a2[i]), fontsize=5, color = 'Blue')

    return(1)

# Title

st.title("Vizualize/Plot Recorded Earthquake Ground Motions")
st.write("V2/V2c files are free-field earthquake records that can be downloaded from Center for Earthquake Engineering Strong Motion CESMD webiste.  Download one free-field record at a time and do not unzip.")
st.write("https://www.strongmotioncenter.org/")
st.write("This app helps read the file and show the recording and create spectra from the recordings")
st.write("Orbit plots and Tripartite Spectra options are included.")
filenames=st.file_uploader("Upload V2/V2c zip file",type=[ "zip"])

V2c = V2 = False
f= None
f_all=[];f_name=[]
stationNo = 0
if filenames != None:
    if filenames.name[-4:]==".zip":
        archive = zipfile.ZipFile(filenames, 'r')
        flist = archive.namelist()
        filenames2=io.BytesIO(archive.read(flist[0]))
        if len(flist) > 1:
            for index,vfl in enumerate(flist):
                if vfl[-4:]==".V2c"or vfl[-4:]==".V2C":
                    f_all.append(io.TextIOWrapper(io.BytesIO(archive.read(vfl))))
                    f_name.append(vfl)
                    V2c =True
                if vfl[-3:]==".v2"or vfl[-3:]==".V2":
                    f=io.TextIOWrapper(io.BytesIO(archive.read(vfl)))
                    V2 =True
                    break
            if len(f_all) == 0 and f == None:
                st.write('Error', 'Zip file does not contain freefield .v2 or .V2c file')
                exit()
        else:
            archive2 = zipfile.ZipFile(filenames2, 'r')
            flist2 = archive2.namelist()
            for index,vfl in enumerate(flist2):
                if vfl[-4:]==".V2c"or vfl[-4:]==".V2C":
                    f_all.append(io.TextIOWrapper(io.BytesIO(archive2.read(vfl))))
                    f_name.append(vfl)
                    V2c = True
                if vfl[-3:]==".v2"or vfl[-3:]==".V2":
                    f=io.TextIOWrapper(io.BytesIO(archive2.read(vfl)))
                    V2 = True
                    break
            if len(f_all) == 0 and f == None:
                st.error('Error', 'Zip file does not contain freefield .v2 or .V2c file')
                exit()
    elif filenames[0][-3:]==".v2" or filenames[0][-3:]==".V2":
        f=open(filenames[0])
        V2 = True
    else:
        st.write('Error', 'V2 File not selected, exiting')
        exit()
    EOF =0
    placeholder = st.empty()
    if V2c:
        EOF = 0
        for index,vfl in enumerate(f_name):
            placeholder.write("Reading V2c file " + str(index))
            if ("HNE" in vfl and "acc" in vfl) or ("HN1" in vfl and "acc" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh1,dtAccel1,numofPointsAccel1,accel1 = readFileV2c(f,f_name[index])
            elif ("HNN" in vfl and "acc" in vfl) or ("HN2" in vfl and "acc" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh2,dtAccel2,numofPointsAccel2,accel2 = readFileV2c(f,f_name[index])
            elif ("HNZ" in vfl and "acc" in vfl) or ("HNZ" in vfl and "acc" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh3,dtAccel3,numofPointsAccel3,accel3 = readFileV2c(f,f_name[index])
            elif ("HNE" in vfl and "vel" in vfl) or ("HN1" in vfl and "vel" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh1,dtVel1,numofPointsVel1,vel1 = readFileV2c(f,f_name[index])
            elif ("HNN" in vfl and "vel" in vfl) or ("HN2" in vfl and "vel" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh2,dtVel2,numofPointsVel2,vel2 = readFileV2c(f,f_name[index])
            elif ("HNZ" in vfl and "vel" in vfl) or ("HNZ" in vfl and "vel" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh3,dtVel3,numofPointsVel3,vel3 = readFileV2c(f,f_name[index])
            elif ("HNE" in vfl and "dis" in vfl) or ("HN1" in vfl and "dis" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh1,dtDispl1,numofPointsDispl1,displ1 = readFileV2c(f,f_name[index])
            elif ("HNN" in vfl and "dis" in vfl) or ("HN2" in vfl and "dis" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh2,dtDispl2,numofPointsDispl2,displ2 = readFileV2c(f,f_name[index])
            elif ("HNZ" in vfl and "dis" in vfl) or ("HNZ" in vfl and "dis" in vfl):
                f = f_all[index]
                recTime,latitude,longitude,nameCh3,dtDispl3,numofPointsDispl2,displ3 = readFileV2c(f,f_name[index])
            else:
                st.write("Error", "File not recognized, exiting")
                exit()
        placeholder.write("Completed reading V2c files")
        unitsAccel1 = unitsAccel2 = unitsAccel3 = "cm/sec2"
        unitsVel1 = unitsVel2 = unitsVel3 = "cm/sec"
        unitsDispl1 = unitsDispl2 = unitsDispl3 = "cm"
    elif V2:
        for line in islice(f, 2, 3):   
            recTime = line[:50].strip()
        # print(recTime)

        for line in islice(f, 2, 3):
            stationNo = line[11:17].strip()    
            latlong= line[17:40].strip()
            latitude =float(latlong[:latlong.find(",")-1])
            if latlong[len(latlong)-1:len(latlong)]=="W":
                longitude =float("-"+latlong[latlong.find(",")+1: len(latlong)-1].strip())
            else:
                longitude =float(latlong[latlong.find(",")+1: len(latlong)-1].strip())
            #print(latitude, longitude)
        for line in islice(f, 17, 18):
            nameCh1=line[26:].strip()
        for line in islice(f, 0, 1):
            nameCh1=nameCh1 + line.strip()
            #print(nameCh1)

        for line in islice(f, 20, 21):
            #print(line)
            numofPointsAccel1 = int(line[0: line.find("points")].strip())
            dtAccel1 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsAccel1 = line[line.find(", in") + 4: line.find(". (")].strip()
        numofLines = lines(numofPointsAccel1)
        accel1 = readchunk(f,numofLines)

        for line in islice(f, 0,1):
            #print(line)
            numofPointsVel1 = int(line[0: line.find("points")].strip())
            dtVel1 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsVel1 = line[line.find(", in") + 4: line.find(".  (")].strip()
        numofLines = lines(numofPointsVel1)
        vel1 = readchunk(f,numofLines)

        for line in islice(f, 0,1):
            #print(line)
            numofPointsDispl1 = int(line[0: line.find("points")].strip())
            dtDispl1 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsDispl1 = line[line.find(", in") + 4: line.find(".   ")].strip()
        numofLines = lines(numofPointsDispl1)
        displ1 = readchunk(f,numofLines)

        EOF=1
        for line in islice(f,1,24):
            if line != "":
                EOF = 0

        for line in islice(f, 0, 1):
            nameCh2=line[26:].strip()
        for line in islice(f, 0, 1):
            nameCh2=nameCh2 + line.strip()
            #print(nameCh2)

        for line in islice(f,1,20):
            i=0

        if EOF ==0:
            for line in islice(f, 0, 1):
                #print(line)
                numofPointsAccel2 = int(line[0: line.find("points")].strip())
                dtAccel2 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
                unitsAccel2 = line[line.find(", in") + 4: line.find(". (")].strip()
            numofLines = lines(numofPointsAccel2)
            accel2 = readchunk(f,numofLines)

            for line in islice(f, 0,1):
                #print(line)
                numofPointsVel2 = int(line[0: line.find("points")].strip())
                dtVel2 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
                unitsVel2 = line[line.find(", in") + 4: line.find(".  (")].strip()
            numofLines = lines(numofPointsVel2)
            vel2 = readchunk(f,numofLines)

            for line in islice(f, 0,1):
                #print(line)
                numofPointsDispl2 = int(line[0: line.find("points")].strip())
                dtDispl2 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
                unitsDispl2 = line[line.find(", in") + 4: line.find(".   ")].strip()
            numofLines = lines(numofPointsDispl2)
            displ2 = readchunk(f,numofLines)  

            for line in islice(f,1,24):
                i=0
            for line in islice(f, 0, 1):
                nameCh3=line[26:].strip()
            for line in islice(f, 0, 1):
                nameCh3=nameCh3 + line.strip()
                #print(nameCh3)
            for line in islice(f,1,20):
                i=0    

            for line in islice(f, 0, 1):
                #print(line)
                numofPointsAccel3 = int(line[0: line.find("points")].strip())
                dtAccel3 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
                unitsAccel3 = line[line.find(", in") + 4: line.find(". (")].strip()
            numofLines = lines(numofPointsAccel3)
            accel3 = readchunk(f,numofLines)

            for line in islice(f, 0,1):
                #print(line)
                numofPointsVel3 = int(line[0: line.find("points")].strip())
                dtVel3 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
                unitsVel3 = line[line.find(", in") + 4: line.find(".  (")].strip()
            numofLines = lines(numofPointsVel3)
            vel3 = readchunk(f,numofLines)

            for line in islice(f, 0,1):
                #print(line)
                numofPointsDispl3 = int(line[0: line.find("points")].strip())
                dtDispl3 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
                unitsDispl3 = line[line.find(", in") + 4: line.find(".   ")].strip()
            numofLines = lines(numofPointsDispl3)
            displ3 = readchunk(f,numofLines)  
        f.close()


    

    T1 = np.arange(0.0,numofPointsAccel1*dtAccel1, dtAccel1)
    scale = scaleValue(unitsAccel1) 
    scaledAccel1 = [value*scale for value in accel1]
    if EOF == 0:
        T2 = np.arange(0.0,numofPointsAccel2*dtAccel2, dtAccel2)
        T3 = np.arange(0.0,numofPointsAccel1*dtAccel3, dtAccel3)
        scale = scaleValue(unitsAccel2) 
        scaledAccel2 = [value*scale for value in accel2]
        scale = scaleValue(unitsAccel3)  
        scaledAccel3 = [value*scale for value in accel3]
    

    st.header(recTime)
    df = pd.DataFrame({"lat":[float(latitude)], "lon":[float(longitude)]})
    st.map(df)   
    c1, c2 =st.columns(2)
    with c1:
        if stationNo != 0 and stationNo.isnumeric():
            st.link_button("See Instrument Details", 'https://www.strongmotioncenter.org/cgi-bin/CESMD/stationhtml.pl?stationID=CE'+stationNo+'&network=CGS')  
    with c2:
        st.link_button("See location of instrument in Google Maps", 'http://www.google.com/maps/place/'+ str(latitude) +','+str(longitude)+'/@'+ str(latitude) +','+str(longitude)+',12z')
    st.subheader("Recorded Values")
    trigger = min(abs(max(accel1, key=abs))/10,abs(max(accel2, key=abs))/10,abs(max(accel3, key=abs))/10,5)

    values = st.sidebar.slider("Select range of times to use", 0.0, dtAccel1*numofPointsAccel1, (startlimAccel(trigger), endlimAccel(trigger)), step= 0.1)
    st.sidebar.caption("*Range autoselected using a trigger of " + str(round(trigger*scaleValue(unitsAccel1),3)) + "g")
    starttime, endtime = values
    width = st.sidebar.slider("plot width", 1, 25, 10)
    height = st.sidebar.slider("plot height", 1, 10, 8)
    doption = st.selectbox("Plot",("Accel", "Vel", "Disp"),)

    if EOF == 1:
        if doption =="Disp":
            rT='Disp (cm)'
            yV = displ1
        elif doption =="Vel":
            rT ='Vel (cm/sec)'
            yV = vel1
        else:
            rT ='Accel (g)'
            yV = scaledAccel1
        fig, ax = plt.subplots(1,1,sharex='col',sharey='all',figsize=(width, height))

        ax.set_xlabel('Time (secs)')
        ax.set_ylabel(rT)
        ax.set_xlim(starttime,endtime)

        ax.set_title(nameCh1)
        ax.plot(T1,yV, label="Channel1", color= 'Red', linewidth=1.0)
        amax=maxaccel(yV, T1); ax.annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(yV, T1); ax.annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
        st.pyplot(fig)
    else:
        fig, ax = plt.subplots(3,1,sharex='col',sharey='all',figsize=(width, height))
        if doption =="Disp":
            rT='Disp (cm)'
            yV1 = displ1; yV2 = displ2; yV3 = displ3
        elif doption =="Vel":
            rT ='Vel (cm/sec)'
            yV1 = vel1; yV2 = vel2; yV3 = vel3
        else:
            rT ='Accel (g)'
            yV1 = scaledAccel1; yV2 = scaledAccel2; yV3 = scaledAccel3; 
        ax[2].set_xlabel('Time (secs)')
        ax[0].set_ylabel(rT)
        ax[1].set_ylabel(rT)
        ax[2].set_ylabel(rT)
        ax[0].set_xlim(starttime,endtime)

        ax[0].set_title(nameCh1)
        ax[0].grid()
        ax[0].plot(T1,yV1, label="Channel1", color= 'Red', linewidth=1.0)
        amax=maxaccel(yV1, T1); ax[0].annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(yV1, T1); ax[0].annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')

        ax[1].set_title(nameCh2)
        ax[1].grid()
        ax[1].plot(T1,yV2, label="Channel2", color= 'Red', linewidth=1.0)
        amax=maxaccel(yV2, T1); ax[1].annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(yV2, T1); ax[1].annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')

        ax[2].set_title(nameCh3)
        ax[2].grid()
        ax[2].plot(T1,yV3, label="Channel3", color= 'Red', linewidth=1.0)
        amax=maxaccel(yV3, T1); ax[2].annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(yV3, T1); ax[2].annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
        st.pyplot(fig)

        st.subheader("Orbit Plots")
        orbitplot = st.checkbox("Create Orbit Plots")
        if orbitplot:
            ooption = st.selectbox("Orbit Plot Type",("Accel", "Vel", "Disp"),)


            if "Up" in nameCh1 or "HNZ" in nameCh1:
                if "360" in nameCh2 or "180" in nameCh2:
                    xa = scaledAccel3.copy(); ya = scaledAccel2.copy(); za = scaledAccel1.copy()
                    xv = vel3.copy(); yv = vel2.copy(); zv = vel1.copy()
                    x = displ3.copy(); y = displ2.copy(); z = displ1.copy()
                    xRec=nameCh3;yRec=nameCh2;zRec=nameCh1
                else:
                    xa = scaledAccel2.copy(); ya = scaledAccel3.copy(); za = scaledAccel1.copy()
                    xv = vel2.copy(); yv = vel3.copy(); zv = vel1.copy()
                    x = displ2.copy(); y = displ2.copy(); z = displ1.copy()
                    xRec=nameCh2;yRec=nameCh3;zRec=nameCh1
            elif "Up" in nameCh2 or "HNZ" in nameCh2:
                if "360" in nameCh1 or "180" in nameCh1:
                    xa = scaledAccel3.copy(); ya = scaledAccel1.copy(); za = scaledAccel2.copy()
                    xv = vel3.copy(); yv = vel1.copy(); zv = vel2.copy()
                    x = displ3.copy(); y = displ1.copy(); z = displ2.copy()
                    xRec=nameCh3;yRec=nameCh1;zRec=nameCh2
                else:
                    xa = scaledAccel1.copy(); ya = scaledAccel3.copy(); za = scaledAccel2.copy()
                    xv = vel1.copy(); yv = vel3.copy(); zv = vel2.copy()
                    x = displ1.copy(); y = displ3.copy(); z = displ2.copy()
                    xRec=nameCh1;yRec=nameCh3;zRec=nameCh2

            elif "Up" in nameCh3 or "HNZ" in nameCh3:
                if "360" in nameCh1 or "180" in nameCh1:
                    xa = scaledAccel2.copy(); ya = scaledAccel1.copy(); za = scaledAccel3.copy()
                    xv = vel2.copy(); yv = vel1.copy(); zv = vel3.copy()
                    x = displ2.copy(); y = displ1.copy(); z = displ3.copy()
                    xRec=nameCh2;yRec=nameCh1;zRec=nameCh3
                else:
                    xa = scaledAccel1.copy(); ya = scaledAccel2.copy(); za = scaledAccel3.copy()
                    xv = vel1.copy(); yv = vel2.copy(); zv = vel3.copy()
                    x = displ1.copy(); y = displ2.copy(); z = displ3.copy()
                    xRec=nameCh1;yRec=nameCh2;zRec=nameCh3
            
            if "360" in yRec:
                yRec = yRec.replace("360 Deg", "NS")
            elif "180" in yRec:
                ya[1,:]=[i*-1 for i in ya[1,:]]
                yv[1,:]=[i*-1 for i in yv[1,:]]
                y[1,:]=[i*-1 for i in y[1,:]]
                yRec = yRec.replace("180 Deg", "NS")
            
            if "90" in xRec:
                xRec = xRec.replace("90 Deg", "EW")
            elif "270" in xRec:
                xa[:]=[i*-1 for i in xa[:]]
                xv[:]=[i*-1 for i in xv[:]]
                x[:]=[i*-1 for i in x[:]]
                xRec = xRec.replace("270 Deg", "EW")
            
            if ooption =="Disp":
                rT='Disp (cm)'
                yV = displ1
                orx = x; ory = y
            elif ooption =="Vel":
                rT ='Vel (cm/sec)'
                yV = vel1
                orx = xv; ory = yv
            else:
                rT ='Accel (g)'
                yV = scaledAccel1
                orx = xa; ory = ya

            fig4, ax = plt.subplots(1,1,figsize=(width, width))
            ax.set_title("Orbit plot for " + ooption)
            orbitplotfn()
            st.pyplot(fig4)
            st.write('Note: Orbit plots can be misleading - points on the curves are resultant at a given time but are not the maximum resultant in that direction')



        st.subheader("Response Spectra")
        tab1, tab2, tab3 = st.tabs(["Type of Spectra", "Damping", "End Period"])
        with tab1:
            option = st.selectbox("Type of Spectra",("Accel", "Vel", "Disp"),)
        with tab2:
            xiStr = st.text_input("Damping values (separate with commas)",str("0.0, 0.02, 0.05"))
            xi =[float(i) for i in xiStr.split(",")]
        with tab3:
            endPeriod = float(st.text_input("End Period for Spectra",str("6.0")))

        tT = np.concatenate( (np.arange(0.05, 0.1, 0.005) , np.arange (0.1, 0.5, 0.01) , np.arange (0.5, 1, 0.02) , np.arange (1, endPeriod, 0.05) ) ) # Time vector for the spectral response
        freq = 1/tT # Frequenxy vector
        df = 1.0/dtAccel1
        
        respsec = st.checkbox("Create Response Spectra")
        if respsec:
            sAll = np.zeros((3,len(xi),len(tT)))

            yaxislimit = round(accelim(scaledAccel1, scaledAccel2, scaledAccel3)*1.1,2)
            nyaxislimit = 0.0 - yaxislimit
            fig2, ax = plt.subplots(3,1,sharex='col',sharey='all',figsize=(width, height))

            if option =="Disp":
                rT='SD'
                rL= 'SD (cm)'
                rU = 'cm'
            elif option =="Vel":
                rT ='SV'
                rL= 'SV (cm/sec)'
                rU = 'cm/sec'
            else:
                rT ='SA'
                rL ='SA (g)'
                rU ='g'

            ax[0].set_title(nameCh1)
            resSpectrafn(accel1,ax[0],rU,rT,0)
            ax[0].legend()
            
            ax[1].set_title(nameCh2)
            resSpectrafn(accel2,ax[1],rU,rT,1)
            ax[1].legend()

            ax[2].set_title(nameCh3)
            resSpectrafn(accel3,ax[2],rU,rT,2)
            ax[2].set_xlabel('Period (secs)')
            ax[2].legend()
            st.pyplot(fig2)

        respec3 = st.checkbox("Create PSA vs Disp Spectra")
        if respec3:
            deflt = int(len(xi)/2)
            dampoption = st.selectbox("Pick one damping ratio",xi,index=deflt)
            fig5, ax = plt.subplots(1,3,figsize=(width, height))
            ax[0].set_title(nameCh1)
            adrs(accel1, ax[0])
            ax[0].set_ylabel('Peak PSA (g)')
            ax[1].set_title(nameCh2)
            adrs(accel2, ax[1])
            ax[2].set_title(nameCh3)
            adrs(accel3, ax[2])
            st.pyplot(fig5)

        respsec2 = st.checkbox("Create Tripartite Spectra")
        if respsec2:
            fig3, ax = plt.subplots(3,1,sharex='col',sharey='all',figsize=(width, height*1.5))
            ax[0].set_title(nameCh1)
            resTripSpectrafn(accel1,ax[0])
            ax[0].legend(loc='upper right')
            ax[1].set_title(nameCh2)
            resTripSpectrafn(accel2,ax[1])
            ax[1].legend(loc='upper right')
            ax[2].set_title(nameCh3)
            resTripSpectrafn(accel3,ax[2]) 
            ax[2].legend(loc='upper right')
            st.pyplot(fig3)



    st.subheader("Download Accelerations")
    wch1 = st.checkbox("Download Acceleration " + nameCh1)
    if EOF == 1:
        wch2 = wch3 = False
    else:
        wch2 = st.checkbox("Download Accelertaion " + nameCh2)
        wch3 = st.checkbox("Download Acceleration " + nameCh3)
    if wch1 or wch2 or wch3:
        text_contents = saveFile()
        st.download_button("Save Acceleration file", text_contents, file_name="accelerations.csv",mime="text/csv",)

    if EOF != 1:
        if respsec:
            st.subheader("Download Response Spectra")
            rch1 = st.checkbox("Download Spectrum " + nameCh1)
            rch2 = st.checkbox("Download Spectrum" + nameCh2)
            rch3 = st.checkbox("Download Spectrum " + nameCh3)
            if rch1 or rch2 or rch3:
                text_contents = rsaveFile()
                st.download_button("Save Response Spectra file", text_contents, file_name="respspectra.csv",mime="text/csv",)
