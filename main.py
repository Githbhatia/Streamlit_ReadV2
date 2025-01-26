import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import zipfile, io
from itertools import islice
from RS_function import RS_function

def readFileV2c(f, f_name):
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

def startlimAccel():
    a1 = next(i for i, x in enumerate(accel1) if abs(x) >5)
    startTime  = max(a1*dtAccel1- 2, 0)
    if EOF==0:
        a2 = next(i for i, x in enumerate(accel2) if abs(x) >5)
        a3 = next(i for i, x in enumerate(accel3) if abs(x) >5)
        startTime  = max(min(a1*dtAccel1,a2*dtAccel2,a3*dtAccel3) - 2, 0)
    return round(startTime,2)

def endlimAccel():
    a1 = next(i for i, x in reversed(list(enumerate(accel1))) if abs(x) >5)
    endTime  = max(a1*dtAccel1+ 2, 0)
    if EOF==0:
        a2 = next(i for i, x in reversed(list(enumerate(accel2))) if abs(x) >5)
        a3 = next(i for i, x in reversed(list(enumerate(accel3))) if abs(x) >5)
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
    index = len(tT)
    while j < index:
        textstring += str(round(tT[j],3))
        if rch1:
            textstring += ", " + str(S1[j]) 
        if rch2:
            textstring += ", " + str(S2[j])
        if rch3:
            textstring += ", " + str(S3[j])
        textstring += "\n"
        j+= 1
    return (textstring)

def accelim(x,y,z):
    xmax = max([abs(i) for i in x])
    ymax = max([abs(i) for i in y])
    zmax = max([abs(i) for i in z])
    return max([xmax,ymax,zmax])

# Title
st.title("Read V2/V2c")

filenames=st.file_uploader("Upload V2/V2c zip file",type=[ "zip"])
V2c = V2 = False
f= None
f_all=[];f_name=[]

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

    if V2c:
        EOF = 0
        for index,vfl in enumerate(f_name):
            print("Reading V2c files")
            if "HNE.--.acc" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh1,dtAccel1,numofPointsAccel1,accel1 = readFileV2c(f,f_name[index])
            if "HNN.--.acc" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh2,dtAccel2,numofPointsAccel2,accel2 = readFileV2c(f,f_name[index])
            if "HNZ.--.acc" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh3,dtAccel3,numofPointsAccel3,accel3 = readFileV2c(f,f_name[index])
            if "HNE.--.vel" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh1,dtVel1,numofPointsVel1,vel1 = readFileV2c(f,f_name[index])
            if "HNN.--.vel" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh2,dtVel2,numofPointsVel2,vel2 = readFileV2c(f,f_name[index])
            if "HNZ.--.vel" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh3,dtVel3,numofPointsVel3,vel3 = readFileV2c(f,f_name[index])
            if "HNE.--.dis" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh1,dtDispl1,numofPointsDispl1,displ1 = readFileV2c(f,f_name[index])
            if "HNN.--.dis" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh2,dtDispl2,numofPointsDispl2,displ2 = readFileV2c(f,f_name[index])
            if "HNZ.--.dis" in vfl:
                f = f_all[index]
                recTime,latitude,longitude,nameCh3,dtDispl3,numofPointsDispl2,displ3 = readFileV2c(f,f_name[index])
        print("Completed reading V2c files")
        unitsAccel1 = unitsAccel2 = unitsAccel3 = "cm/sec2"
        unitsVel1 = unitsVel2 = unitsVel3 = "cm/sec"
        unitsDispl1 = unitsDispl2 = unitsDispl3 = "cm"
    elif V2:
        for line in islice(f, 2, 3):   
            recTime = line[:50].strip()
        # print(recTime)

        for line in islice(f, 2, 3):    
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
    starttime = float(st.text_input("Start Time",str(startlimAccel())))
    endtime = float(st.text_input("End Time",str(endlimAccel())))
    st.subheader("Recorded Accelerations")
    width = st.sidebar.slider("plot width", 1, 25, 10)
    height = st.sidebar.slider("plot height", 1, 10, 8)
    
    
    if EOF == 1:
        fig, ax = plt.subplots(1,1,sharex='col',sharey='all',figsize=(width, height))

        ax.set_xlabel('Time (secs)')
        ax.set_ylabel('Accel(g)')
        ax.set_xlim(starttime,endtime)

        ax.set_title(nameCh1)
        ax.plot(T1,scaledAccel1, label="Channel1", color= 'Red', linewidth=1.0)
        amax=maxaccel(scaledAccel1, T1); ax.annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(scaledAccel1, T1); ax.annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
        st.pyplot(fig)
    else:
        fig, ax = plt.subplots(3,1,sharex='col',sharey='all',figsize=(width, height))

        ax[2].set_xlabel('Time (secs)')
        ax[0].set_ylabel('Accel(g)')
        ax[1].set_ylabel('Accel(g)')
        ax[2].set_ylabel('Accel(g)')
        ax[0].set_xlim(starttime,endtime)

        ax[0].set_title(nameCh1)
        ax[0].plot(T1,scaledAccel1, label="Channel1", color= 'Red', linewidth=1.0)
        amax=maxaccel(scaledAccel1, T1); ax[0].annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(scaledAccel1, T1); ax[0].annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')

        ax[1].set_title(nameCh2)
        ax[1].plot(T1,scaledAccel2, label="Channel2", color= 'Red', linewidth=1.0)
        amax=maxaccel(scaledAccel2, T1); ax[1].annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(scaledAccel2, T1); ax[1].annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')

        ax[2].set_title(nameCh3)
        ax[2].plot(T1,scaledAccel3, label="Channel3", color= 'Red', linewidth=1.0)
        amax=maxaccel(scaledAccel3, T1); ax[2].annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(scaledAccel3, T1); ax[2].annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
        st.pyplot(fig)


        st.subheader("Response Spectra")
        respsec = st.checkbox("Create Response Spectra")
        if respsec:
            xi = float(st.text_input("Damping",str("0.05")))
            endPeriod = float(st.text_input("End Period for Spectra",str("6.0")))
            option = st.selectbox("Type of Spectra",("Accel", "Vel", "Disp"),)
            tT = np.concatenate( (np.arange(0.05, 0.1, 0.005) , np.arange (0.1, 0.5, 0.01) , np.arange (0.5, 1, 0.02) , np.arange (1, endPeriod, 0.05) ) ) # Time vector for the spectral response
            freq = 1/tT # Frequenxy vector
            df = 1.0/dtAccel1
            Sfin=[]

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
            ax[0].grid()
            Sfin= RS_function(accel1[int(starttime/dtAccel1):int(endtime/dtAccel1)], df, tT, xi, Resp_type = rT)
            if option =="Accel":
                S1=Sfin[0,:]*scaleValue(unitsAccel1)
            else:
                S1=Sfin[0,:]
            ax[0].set_ylabel(rL)
            ax[0].plot(tT,S1,color= 'Red', linewidth=1.0)
            amax=[tT[np.argmax(abs(S1))], max(abs(S1))]; ax[0].annotate(str(round(amax[0],3)) +"sec, "+str(round(amax[1],2)) + rU , xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            ax[0].text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=8, color ='Black',transform=ax[0].transAxes)
            
            ax[1].set_title(nameCh2)
            ax[1].grid()
            Sfin= RS_function(accel2[int(starttime/2):int(endtime/dtAccel2)], df, tT, xi, Resp_type = rT)
            if option =="Accel":
                S2=Sfin[0,:]*scaleValue(unitsAccel2)
            else:
                S2=Sfin[0,:]
            ax[1].set_ylabel(rL)
            ax[1].plot(tT,S2,color= 'Red', linewidth=1.0)
            amax=[tT[np.argmax(abs(S2))], max(abs(S2))]; ax[1].annotate(str(round(amax[0],3)) +"sec, "+str(round(amax[1],2)) + rU , xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            ax[1].text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=8, color ='Black',transform=ax[1].transAxes)

            ax[2].set_title(nameCh3)
            ax[2].grid()
            Sfin= RS_function(accel3[int(starttime/dtAccel3):int(endtime/dtAccel3)], df, tT, xi, Resp_type = rT)
            if option =="Accel":
                S3=Sfin[0,:]*scaleValue(unitsAccel1)
            else:
                S3=Sfin[0,:]
            ax[2].set_xlabel('Period (secs)')
            ax[2].set_ylabel(rL)
            ax[2].plot(tT,S3,color= 'Red', linewidth=1.0)
            amax=[tT[np.argmax(abs(S3))], max(abs(S3))]; ax[2].annotate(str(round(amax[0],3)) +"sec, "+str(round(amax[1],2)) + rU , xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            ax[2].text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=8, color ='Black',transform=ax[2].transAxes)
            st.pyplot(fig2)

 
    st.subheader("Download Accelerations")
    wch1 = st.checkbox("Download Acceleration " + nameCh1)
    wch2 = st.checkbox("Download Accelertaion " + nameCh2)
    wch3 = st.checkbox("Download Acceleration " + nameCh3)
    if wch1 or wch2 or wch3:
        text_contents = saveFile()
        st.download_button("Save Acceleration file", text_contents)

    if respsec:
        st.subheader("Download Response Spectra")
        rch1 = st.checkbox("Download Spectrum " + nameCh1)
        rch2 = st.checkbox("Download Spectrum" + nameCh2)
        rch3 = st.checkbox("Download Spectrum " + nameCh3)
        if rch1 or rch2 or rch3:
            text_contents = rsaveFile()
            st.download_button("Save Response Spectra file", text_contents)


