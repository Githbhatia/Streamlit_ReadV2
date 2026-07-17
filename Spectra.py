
import numpy as np
import certifi
import ssl
import geopy.geocoders
from geopy.geocoders import Nominatim
from geopy.geocoders import ArcGIS
from geopy.exc import GeocoderTimedOut
from geopy.extra.rate_limiter import RateLimiter
import urllib.request as ur
import json as js
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd
import math
import pydeck as pdk

def persistent_checkbox(label, key):
    state = st.checkbox(label, value=st.session_state.checklist_items.get(key, False), key=key)
    st.session_state.checklist_items[key] = state
    return state


@st.cache_resource
def myurlopen(url): 
    ctx = ssl.create_default_context(cafile=certifi.where())
    try:
        response = ur.urlopen(url)
    except ur.URLError as e:
        if hasattr(e, 'reason'):
            st.write('We failed to reach a server.')
            st.write('Reason: ', e.reason)
            return()
        elif hasattr(e, 'code'):
            st.write('The server couldn\'t fulfill the request.')
            st.write('Error code: ', e.code)
            return() 
    
    return(response.read())

@st.cache_resource
def mygeolocatorreverse(lat, longt):
    ctx = ssl.create_default_context(cafile=certifi.where())
    # ctx = ssl._create_unverified_context()
    geopy.geocoders.options.default_ssl_context = ctx
    if whichgeolocator == "OpenStreetMaps (Default)":
        geolocator = Nominatim(user_agent="STASCE722SpectraFp1", scheme='https')
    else:
        geolocator = ArcGIS(user_agent="STASCE722SpectraFp1", scheme='https')
    geocode = RateLimiter(geolocator.reverse, min_delay_seconds=1)
    try:
        location = geolocator.reverse(str(lat) + " ," + str(longt), timeout=2)
        if location != None:
            address = str(location.address)
            if whichgeolocator == "OpenStreetMaps (Default)":
                st.write("Using "+ address + " (Geocoding services provided by OpenStreetMaps)")
            else:
                st.write("Using "+ address + " (Geocoding services provided by ArcGIS)")
            return(location.address)
        else:
            st.write("Address not found: Continuing using "+ str(lat) + ", " + str(longt)) 
            return("")  
    except GeocoderTimedOut as e:
        st.write("Error: geocode failed on input %s with message %s"%(address, e.message))
        st.write("Continuing using "+ str(lat) + ", " + str(longt))
        return("")

@st.cache_resource
def mygeolocator(address):
    ctx = ssl.create_default_context(cafile=certifi.where())
    # ctx = ssl._create_unverified_context()
    geopy.geocoders.options.default_ssl_context = ctx
    # geolocator = Nominatim(user_agent="STASCE722SpectraFp2")
    if whichgeolocator == "OpenStreetMaps (Default)":
        geolocator = Nominatim(user_agent="STASCE722SpectraFp1", scheme='https')
    else:
        geolocator = ArcGIS(user_agent="STASCE722SpectraFp1", scheme='https')
    geocode = RateLimiter(geolocator.geocode, min_delay_seconds=1)
    try:
        location = geolocator.geocode(address, timeout=2)
        if (location != None):
            lat = str(location.latitude)
            longt = str(location.longitude)
            address = str(location.address)
            if whichgeolocator == "OpenStreetMaps (Default)":
                st.write("Using "+ str(lat) + ", " + str(longt) + " (Geocoding services provided by OpenStreetMaps)")
            else:
                st.write("Using "+ str(lat) + ", " + str(longt) + " (Geocoding services provided by ArcGIS)")
            return(location.latitude, location.longitude, location.address)
        else:
            st.write("Invalid Address:" + "Revise address and try again")
            return(0.0, 0.0, "")
    except GeocoderTimedOut as e:
        st.write("Error: geocode failed on input %s with message %s"%(address, e.message))
        return(0.0, 0.0, "")


    

def onclick():

    global address, lat,longt, textout, riskct, sitecl,sds, whichgeolocator
   
    if swv != 0.0:
        try:
            shearwavevel = float(swv)
        except ValueError:
            st.write("Invalid Shear Wave Velocity:"+ "Enter shear wave velocity in ft/sec and try again")
            return
        if shearwavevel==0:
            st.write("Invalid Shear Wave Velocity:"+ "Enter a non-zero shear wave velocity in ft/sec and try again")
            return
        shearwavevellimits = [('F',0.0),('E',500.0),('DE',700.0),('D',1000.0),('CD',1450.0),('C',2100.0),('BC',3000.0),('B',5000.0),('A',1000000.0)]
        centershearwave = [('E',500.0),('DE',600.0),('D',849.0),('CD',1200.0),('C',1732.0),('BC',2500.0),('B',3536.0),('A',1000000.0)]
        index = 0
        for a, b in shearwavevellimits:
            if shearwavevel <= b:
                sitecl = a
                break
        prev =0
        for a, b in centershearwave:
            if shearwavevel > b:
                sitecll = a
                prev = b
        if shearwavevel <= 500.0:
            sitecll = "E"
                
        for a, b in centershearwave:
            if shearwavevel <= b:
                siteclu = a
                siteclBMultp = (shearwavevel - prev)/(b- prev)
                break
                
        if estimatedswv==1:
            for a, b in shearwavevellimits:
                if shearwavevel/1.3 <= b:
                    sitecll = a
                    break
            for a, b in shearwavevellimits:
                if shearwavevel*1.3 <= b:
                    siteclu = a
                    break

        placeholder.selectbox("Site Class",siteClassList,index = siteClassList.index(sitecl)) 
            
        
    # elif siteclass=="Default":
    #     sitecl = "CD"
    #     siteclu = "C"
    #     sitecll = "D"
    else:
        sitecl = siteclass
    if sitecl == 'F': 
        st.write("Invalid Shear Wave Velocity:" + "Site Class F, Requires site response analysis studies")
        return(0)
    
    #print(sitecll+" "+siteclu)
    #print(st.session_state.siteclass)

    if siteclass=="Default" and swv == 0.0:
        st.write("Using site class " + siteclass)
    else:
        st.write("Using site class " + sitecl)

    sitetitle = mysite
    riskct = riskc
    address = addressg

    

    if address =="":
        lat = latitude
        longt = longitude
        try:
            address = mygeolocatorreverse(lat, longt)
        except Exception as e:
            st.write("Geolocator not available, try again " + str(e))
    else:
        try:
            lat, longt, address = mygeolocator(address)
        except Exception as e:
            st.write("Geolocator not available, try again " + str(e))
            st.stop()
            

    
    df = pd.DataFrame({"lat":[float(lat)], "lon":[float(longt)],"text": mysite})
    view = pdk.ViewState(
        latitude=float(lat),
        longitude=float(longt),
        zoom=11,)
    st.pydeck_chart(
    pdk.Deck(
        # map_style="mapbox://styles/mapbox/light-v9",
        initial_view_state=view,
        tooltip={"text": "{text}"},
        layers=[
            pdk.Layer(
                "ScatterplotLayer",
                data=df,
                get_position=["lon", "lat"],
                get_color=[255,0,0],
                get_radius=10,
                radiusMinPixels=5,
                radiusMaxPixels=50,
                pickable=True,
                ),
            pdk.Layer(
                "TextLayer",
                data=df,
                get_position=["lon", "lat"],
                get_color=[0, 0, 0],
                get_text=str("text"),
                get_size=11,
                get_text_anchor='"middle"',
                get_alignment_baseline='"top"',
                pickable=True,
                ),


            ],
        )
    )
  
  
    faultURL = "https://hcai.maps.arcgis.com/apps/webappviewer/index.html?id=783b89c5f62a4853a92f08a320d6d518&marker=" + \
        str(longt) + ";" + str(lat)+ ";;"+str(mysite)+\
            "&scale=200000&showLayers=Earthquake Faults and Folds in the USA - Qfaults_US_Database;Earthquake Faults and Folds in the USA - ca_offshore" 
    st.link_button("See proximity of site to seismic sources and faults", faultURL, type="primary")

    url = 'https://earthquake.usgs.gov/ws/building-codes/asce7-22/calculate?latitude='+ str(lat) + '&longitude=' + str(longt) +'&riskCategory='+ riskct +'&siteClass=' + sitecl + '&title=Example'
    
    # if  swv != 0.0 or siteclass=="Default":
    if swv != 0.0:
        urll = 'https://earthquake.usgs.gov/ws/building-codes/asce7-22/calculate?latitude='+ str(lat) + '&longitude=' + str(longt) +'&riskCategory='+ riskct +'&siteClass=' + sitecll + '&title=Example'
        urlu = 'https://earthquake.usgs.gov/ws/building-codes/asce7-22/calculate?latitude='+ str(lat) + '&longitude=' + str(longt) +'&riskCategory='+ riskct +'&siteClass=' + siteclu + '&title=Example'
        


    response = myurlopen(url)
    # if swv != 0.0 or siteclass=="Default":
    if swv != 0.0 :
        responsel = myurlopen(urll)
        responseu = myurlopen(urlu)

    
    rdata = js.loads(response)
    # if swv != 0.0 or siteclass=="Default":      
    if swv != 0.0:      
        rdatal = js.loads(responsel)
        rdatau = js.loads(responseu)

    # if self.SaveJson.get() == 1:
    #     with open("ASCE722.json", "w") as write_file:
    #         js.dump(rdata, write_file)
    #     if str(self.entry_SWVel.get()) != "" or str(self.SelectedSiteClass.get())=="Default":
    #         with open("ASCE722_lowerbound.json", "w") as write_file:
    #             js.dump(rdatal, write_file)
    #         with open("ASCE722_upperbound.json", "w") as write_file:
    #             js.dump(rdatau, write_file)

    output = 'Output for Latitude = ' + str(lat) + ' Longitude = ' + str(longt)
    t = rdata["response"]["data"]["multiPeriodDesignSpectrum"]["periods"]
    s = rdata["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"]
    # st.write (t, s )

    t2 = rdata["response"]["data"]["twoPeriodDesignSpectrum"]["periods"]
    s2 = rdata["response"]["data"]["twoPeriodDesignSpectrum"]["ordinates"]
        
    tmce = rdata["response"]["data"]["multiPeriodMCErSpectrum"]["periods"]
    smce = rdata["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]



    tmce2 = rdata["response"]["data"]["twoPeriodMCErSpectrum"]["periods"]
    smce2 = rdata["response"]["data"]["twoPeriodMCErSpectrum"]["ordinates"]

    # if swv != 0.0 or siteclass=="Default":    
    if swv != 0.0:
        tl = rdatal["response"]["data"]["multiPeriodDesignSpectrum"]["periods"]
        sl = rdatal["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"]
        
        tu = rdatau["response"]["data"]["multiPeriodDesignSpectrum"]["periods"]
        su = rdatau["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"]

        tmcel = rdatal["response"]["data"]["multiPeriodMCErSpectrum"]["periods"]
        smcel = rdatal["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]

        tmceu = rdatau["response"]["data"]["multiPeriodMCErSpectrum"]["periods"]
        smceu = rdatau["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]

    fig = plt.figure(figsize=(10, 10))
    fig2 = plt.figure(figsize=(10, 5))
    ax = fig.subplots(2,1)
    ax2 = fig2.subplots(1,1)
    ax[0].set_xlabel('Period')
    ax[0].set_title(sitetitle + " Design Spectrum")

    ax[1].set_xlabel('Period')
    ax[1].set_title(sitetitle + " MCE Spectrum")

    ax2.set_xlabel('Period')

    if (estimatedswv and swv != 0.0):

        sg = [max(sl,s,su) for sl,s,su in zip(sl,s,su)]
        ax[0].plot(t, sl, label="Multiperiod Des Spec lower bound SC= "+ sitecll, color='Red', linewidth=1.0)
        ax[0].plot(t, s, label="Multiperiod Des Spec SC= " + sitecl, color='Blue', linewidth=1.0)
        ax[0].plot(t, su, label="Multiperiod Des Spec upper bound SC= "+ siteclu, color='Green', linewidth=1.0)
        ax[0].plot(t, sg, label="Govering Multiperiod Des Spec", color='Black', linestyle='--', linewidth=2.0)
        ax[0].set_xlim([0, 5])
        ax[0].legend()  
        ax[0].grid()
        smcel = rdatal["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]
        smceu = rdatau["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]
        smceg = [max(smcel,smce,smceu) for smcel,smce,smceu in zip(smcel,smce,smceu)]
        ax[1].plot(tmce, smcel, label="MCE Multiperiod lower bound SC= "+ sitecll, color='Red', linewidth=1.0)
        ax[1].plot(tmce, smce, label="MCE Multiperiod Spec SC= " + sitecl, color='Blue', linewidth=1.0)
        ax[1].plot(tmce, smceu, label="MCE Multiperiod upper bound SC= "+ siteclu, color='Green', linewidth=1.0)
        ax[1].plot(tmce, smceg, label="Govering MCE Multiperiod", color='Black', linestyle='--', linewidth=2.0)
        ax[1].set_xlim([0, 5])
        ax[1].legend() 
        ax[1].grid()
        

        
        sds = 0.9 * max(sg[t.index(0.2):t.index(5.0)])
        st.session_state['sds'] = sds
        sd1min = sg[t.index(1.0)]
        sd1 = 0.0
        if shearwavevel > 1450:
            for i in range(t.index(1.0), t.index(2.0)+1):
                sd1 = max(0.9*sg[i]*t[i], sd1)
            sd1=max(sd1,sd1min)
        elif shearwavevel <= 1450:
            for i in range(t.index(1.0), t.index(5.0)+1):
                sd1 = max(0.9*sg[i]*t[i], sd1)
            sd1=max(sd1,sd1min)

        elfs = [min(x,sds) for x in sg ]
        focc = elfs.index(sds)
        elfs = [sds if ind <= focc else x for ind, x in enumerate(elfs) ]
        ax2.plot(t, elfs, label="Recommended ELF Design Spectrum", color='Purple', linewidth=1.0)
        ax2.set_xlim([0, 5])
        ax2.legend()
        ax2.grid()
        ax2.set_title(sitetitle + " Recommended ELF Design Spectrum (MPRS capped to SDS)")


        st.subheader("ASCE7-22 Seismic Parameter Output")
        st.write("Based on est. shear wave velocity per ASCE 7-22 Section 20.3 and 21.4")
        df = pd.DataFrame(
        {'Parameter':["sms","sm1","sds","sd1","pga"],'Values':[str(round(sds*1.5,3)),str(round(sd1*1.5,3)),str(round(sds,3)),str(round(sd1,3)),str(round(sg[0],3))]}
        )
        df.set_index('Parameter', inplace=True)
        st.dataframe(df)

        dfs=pd.DataFrame({"time period":t,"Governing Multiperiod Spec": sg, "Recommended ELF Design Spectrum": elfs,"Governing MCE Multiperiod":smceg})
        st.dataframe(dfs)
        textout = mywritefileEstSV(t, sg, tmce, smceg, sds, sd1, sitecl, elfs)



    elif  swv != 0.0:
        sexp = np.array(su)*siteclBMultp + np.array(sl)*(1-siteclBMultp)
        sexpmce = np.array(smceu)*siteclBMultp + np.array(smcel)*(1-siteclBMultp)
        ax[0].plot(t, s, label="Multiperiod Design Spectrum for " + sitecl, color='Red', linewidth=1.0)
        ax[0].plot(t2, s2, label="2-Period Design Spectrum for " + sitecl, color='Green', linewidth=1.0)
        #ax[0].plot(tl, sl, label="Lower Bound Design Spectrum for" + sitecll, color='black', linewidth=0.1)
        ax[0].plot(tl, sexp, label="Interpolated Spectrum for " + str(round(shearwavevel,0)) + " ft/s", color='black', linestyle='--', linewidth=1.0)
        ax[0].set_xlim([0, 5])
        ax[0].legend()
        ax[0].grid()
        ax[1].plot(tmce, smce, label="MCE Multiperiod Spectrum", color='Blue', linewidth=1.0)
        ax[1].plot(tmce2, smce2, label="MCE 2-Period  Spectrum", color='Green', linewidth=1.0)
        ax[1].plot(tmcel, sexpmce, label="Interpolated mCE Spectrum for " + str(round(shearwavevel,0)) + " ft/s", color='black', linestyle='--', linewidth=1.0)
        ax[1].set_xlim([0, 5])
        ax[1].legend()
        ax[1].grid()
        p = rdata["response"]["data"].items()
        sds = rdata["response"]["data"]["sds"]

        elfs = [min(x,sds) for x in s ]
        focc = elfs.index(sds)
        elfs = [sds if ind <= focc else x for ind, x in enumerate(elfs) ]
        ax2.plot(t, elfs, label="Recommended ELF Design Spectrum", color='Purple', linewidth=1.0)
        ax2.set_xlim([0, 5])
        ax2.legend()
        ax2.grid()
        ax2.set_title(sitetitle + " Recommended ELF Design Spectrum (MPRS capped to SDS)")

        st.subheader("ASCE7-22 Seismic Parameter Output")
        df = pd.DataFrame(p)
        df = df[0:11]
        df.columns = ['Parameter','Values']
        df['Values'] = df['Values'].astype(str)
        df.set_index('Parameter', inplace=True)
        st.dataframe(df)
        dfs=pd.DataFrame({"time period":t,"Multiperiod Spec": s, "Interpolated Spec": sexp, "Recommended ELF Design Spectrum": elfs, "MCE Multiperiod":smce, "Interpolated MCE spec": sexpmce })
        st.dataframe(dfs)
        sds =float(df.loc["sds"].values[0])
        st.session_state['sds'] = sds
        textout = mywritefileest(rdata, sitecl, sexp, elfs)
        


    # elif siteclass=="Default":   

    #     sg = [max(sl,s,su) for sl,s,su in zip(sl,s,su)]
    #     ax[0].plot(t, sl, label="Multiperiod Des Spec lower bound SC= "+ sitecll, color='Red', linewidth=1.0)
    #     ax[0].plot(t, s, label="Multiperiod Des Spec SC= " + sitecl, color='Blue', linewidth=1.0)
    #     ax[0].plot(t, su, label="Multiperiod Des Spec upper bound SC= "+ siteclu, color='Green', linewidth=1.0)
    #     ax[0].plot(t, sg, label="Govering Multiperiod Des Spec", color='Black', linestyle='--', linewidth=2.0)
    #     ax[0].set_xlim([0, 5])
    #     ax[0].legend()  
    #     ax[0].grid()
    #     smcel = rdatal["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]
    #     smceu = rdatau["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]
    #     smceg = [max(smcel,smce,smceu) for smcel,smce,smceu in zip(smcel,smce,smceu)]
    #     ax[1].plot(tmce, smcel, label="MCE Multiperiod lower bound SC= "+ sitecll, color='Red', linewidth=1.0)
    #     ax[1].plot(tmce, smce, label="MCE Multiperiod Spec SC= " + sitecl, color='Blue', linewidth=1.0)
    #     ax[1].plot(tmce, smceu, label="MCE Multiperiod upper bound SC= "+ siteclu, color='Green', linewidth=1.0)
    #     ax[1].plot(tmce, smceg, label="Govering MCE Multiperiod", color='Black', linestyle='--', linewidth=2.0)
    #     ax[1].set_xlim([0, 5])
    #     ax[1].legend() 
    #     ax[1].grid()
 

    #     sds = 0.9 * max(sg[t.index(0.2):t.index(5.0)])
    #     st.session_state['sds'] = sds
    #     sd1 = sg[t.index(1.0)]

    #     elfs = [min(x,sds) for x in sg ]
    #     focc = elfs.index(sds)
    #     elfs = [sds if ind <= focc else x for ind, x in enumerate(elfs) ]
    #     ax2.plot(t, elfs, label="Recommended ELF Design Spectrum", color='Purple', linewidth=1.0)
    #     ax2.set_xlim([0, 5])
    #     ax2.legend()
    #     ax2.grid()
    #     ax2.set_title(sitetitle + " Recommended ELF Design Spectrum (MPRS capped to SDS)")

    #     st.subheader("ASCE7-22 Seismic Parameter Output")
    #     st.write("Default = Max of Site Class C, CD, D")
    #     df = pd.DataFrame(
    #     {'Parameter':["sms","sm1","sds","sd1","pga"],'Values':[str(round(sds*1.5,3)),str(round(sd1*1.5,3)),str(round(sds,3)),str(round(sd1,3)),str(round(sg[0],3))]}
    #     )
    #     df.set_index('Parameter', inplace=True)
    #     st.dataframe(df)
    #     dfs=pd.DataFrame({"time period":t,"Governing Multiperiod Spec": sg, "Recommended ELF Design Spectrum": elfs, "Governing MCE Multiperiod":smceg })
    #     st.dataframe(dfs)
    #     textout = mywritefileEstSV(t, sg, tmce, smceg, sds, sd1, sitecl, elfs)


    else:
        ax[0].plot(t, s, label="Multiperiod Design Spectrum for " + sitecl, color='Red', linewidth=1.0)
        ax[0].plot(t2, s2, label="2-Period Design Spectrum for " + sitecl, color='Green', linewidth=1.0)
        ax[0].set_xlim([0, 5])
        ax[0].legend()
        ax[0].grid()
        ax[1].plot(tmce, smce, label="MCE Multiperiod Spectrum", color='Blue', linewidth=1.0)
        ax[1].plot(tmce2, smce2, label="MCE 2-Period  Spectrum", color='Green', linewidth=1.0)
        ax[1].set_xlim([0, 5])
        ax[1].legend()
        ax[1].grid()
        p = rdata["response"]["data"].items()
        sds = rdata["response"]["data"]["sds"]

        elfs = [min(x,sds) for x in s ]
        focc = elfs.index(sds)
        elfs = [sds if ind <= focc else x for ind, x in enumerate(elfs) ]
        ax2.plot(t, elfs, label="Recommended ELF Design Spectrum", color='Purple', linewidth=1.0)
        ax2.set_xlim([0, 5])
        ax2.legend()
        ax2.grid()
        ax2.set_title(sitetitle + " Recommended ELF Design Spectrum (MPRS capped to SDS)")

        st.subheader("ASCE7-22 Seismic Parameter Output")
        df = pd.DataFrame(p)
        df = df[0:11]
        df.columns = ['Parameter','Values']
        df['Values'] = df['Values'].astype(str)
        df.set_index('Parameter', inplace=True)
        st.dataframe(df)
        sds = float(df.loc["sds"].values[0])
        st.session_state['sds'] = sds
        dfs=pd.DataFrame({"time period":t,"Multiperiod Spec": s, "Recommended ELF Design Spectrum": elfs, "MCE Multiperiod":smce })
        st.dataframe(dfs)
        textout = mywritefile(rdata, sitecl, elfs)
        


    st.pyplot(fig)
    st.pyplot(fig2)

    return()

def contourf(lat, longt, riskct):
    sitecl = siteclass
    nlong = 7
    nlat= 7
    gridspacing = 0.5/60.0
    lat = float(lat)
    longt = float(longt)
    latgrid = np.arange(lat+(nlat//2)*gridspacing, lat-((nlat//2)+0.9)*gridspacing, -gridspacing)
    longgrid = np.arange(longt-(nlong//2)*gridspacing, longt+((nlong//2)+0.9)*gridspacing, gridspacing)
    xLong,xLat = np.meshgrid(longgrid,latgrid)
    ZSDS=np.zeros((nlong,nlat)); ZSD1=np.zeros((nlong,nlat))
    
    df = pd.DataFrame({"lat":xLat.flatten(), "lon":xLong.flatten()})
    mesg = st.empty()

    for i in range(nlong):
        for j in range(nlat):
            mesg.write("Getting gird " + str(i) + ", " + str(j))
            url = 'https://earthquake.usgs.gov/ws/building-codes/asce7-22/calculate?latitude='+ str(xLat[i,j]) + '&longitude=' + str(xLong[i,j]) +'&riskCategory='+ riskct +'&siteClass=' + sitecl + '&title=Example'
            response = myurlopen(url)
            rdata = js.loads(response)
            ZSDS[i,j] = rdata["response"]["data"]["sds"]
            ZSD1[i,j] = rdata["response"]["data"]["sd1"]
    mesg.write("Completed")
    #print(ZSDS, ZSD1)
    fig = plt.figure(figsize=(10, 20))
    ax = fig.add_subplot(211)
    CS = ax.contour(xLong,xLat,ZSDS) 
    ax.set_title('Local Variation of SDS around site')
    ax.text(longt,lat , '. Site '+ str(ZSDS[nlong//2, nlat//2]), fontsize = 10)
    ax.clabel(CS, inline=True, fontsize=10)
    ax = fig.add_subplot(212)
    CS2 = ax.contour(xLong,xLat,ZSD1) 
    ax.set_title('Variation of SD1 around site')
    ax.text(longt, lat, '. Site '+ str(ZSD1[nlong//2, nlat//2]), fontsize = 10)
    ax.clabel(CS2, inline=True, fontsize=10)
    textlist = [""]*49
    textlist[24] = "Site"

    df = pd.DataFrame({"lat":xLat.flatten(), "lon":xLong.flatten(), "weight":ZSDS.flatten(), "weight2":ZSD1.flatten(), "text": textlist})
    df.insert(1, "latlong", (df["lat"].round(3)).astype(str)+","+(df["lon"].round(3)).astype(str))
    view = pdk.data_utils.compute_view(df[["lon", "lat"]])
    st.write("Grid Used:(hover to see values)")
    st.pydeck_chart(
    pdk.Deck(
        #map_style="mapbox://styles/mapbox/light-v9",
        initial_view_state=view,
        tooltip={"text": "{latlong}, \n SDS={weight}, SD1={weight2}"},
        layers=[
            pdk.Layer(
                "ScatterplotLayer",
                data=df,
                get_position=["lon", "lat"],
                get_color=[255,0,0],
                get_radius=10,
                radiusMinPixels=5,
                radiusMaxPixels=50,
                pickable=True,
                ),
            pdk.Layer(
                "TextLayer",
                data=df,
                get_position=["lon", "lat"],
                get_color=[0, 0, 0],
                get_text=str("text"),
                get_size=11,
                get_text_anchor='"middle"',
                get_alignment_baseline='"top"',
                pickable=True,
                ),


            ],
        )
    )

    st.pyplot(fig)



def mywritefileEstSV(t, sg, tmce, smceg, sds, sd1, sitecl, elfs):
    sitetitle = mysite
    riskct = riskc

    textout = ""
    textout += versionstr + "\n"
    textout += "Data source is USGS (ASCE 722 Database) and OpenStreetMaps.\nAuthors do not assume any responsibility or liability for its accuracy.\n"
    textout += "Use of the output of this program does not imply approval by the governing building code bodies responsible for building code approval and interpretation for the building site described by latitude/longitude location.\n"
    textout += "\n \n"
    textout += sitetitle + "\n" + address + "\n"
    textout += "The location is " + str(lat) + ", " + str(longt) +  " and Risk Category "+ riskct + "\n"
    if (estimatedswv and swv == 0.0):
        textout += "Site Class based on an estimated shear wave velocity of " + str(swv) + "ft/s\n"
        textout += "Lower bound and upper bound site class considered in computation per ASCE 7-22 Section 20.3 and 21.4" + "\n"
    else:
        textout += "Default Site Class based on max of Site Class C, CD, D\n"
    textout += "sms from governing design spectra = " + str(round(sds*1.5, 3)) + "\n"
    textout += "sm1 from governing design spectra = " + str(round(sd1*1.5, 3)) + "\n"
    textout += "sds from governing design spectra = " + str(round(sds, 3)) + "\n"
    textout += "sd1 from governing design spectra = " + str(round(sd1, 3)) + "\n"
    textout += "pga from governing design spectra = " + str(round(sg[0], 3)) + "\n"
    textout += "Governing MultiPeriodDesignSpectrum\n"
    index = len(t)
    j = 0
    while j < index:
        textout += str(t[j])+ ", " + str(sg[j])+"\n"
        j+= 1
    j = 0
    textout += "Recommended ELF Design Spectrum\n"
    while j < index:
        textout += str(t[j])+ ", " + str(elfs[j])+"\n"
        j+= 1
    textout += "Governing MultiPeriodMCErSpectrum\n"
    index = len(tmce)
    j = 0
    while j < index:
        textout += str(tmce[j])+ ", " + str(smceg[j])+"\n"
        j+= 1
    return(textout)




def mywritefile( ldata, sitecl, elfs):
    sitetitle = mysite
    riskct = riskc

    textout = ""
    index = 0
    p = ldata["response"]["data"]
    t = ldata["response"]["data"]["multiPeriodDesignSpectrum"]["periods"]
    s = ldata["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"]
    tmce = ldata["response"]["data"]["multiPeriodMCErSpectrum"]["periods"]
    smce = ldata["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]
    textout += versionstr + "\n"
    textout += "Data source is USGS (ASCE 722 Database) and OpenStreetMaps.\nAuthors do not assume any responsibility or liability for its accuracy.\n"
    textout += "Use of the output of this program does not imply approval by the governing building code bodies responsible for building code approval and interpretation for the building site described by latitude/longitude location.\n"
    textout += "\n \n"
    textout += sitetitle + "\n" + address + "\n"
    textout += "The location is " + str(lat) + ", " + str(longt) + " with Site Class " + sitecl + " and Risk Category "+ riskct + "\n"
    if swv != 0.0:
        textout += "Site Class based on a shear wave velocity of " + str(swv) + "ft/s\n"
    textout += "pga from design spectra = " + str(round(s[0], 3)) + "\n"
    for key, value in p.items():
        if index <= 11:
            textout += str(key)+ ", " + str(value)+"\n"  
        index += 1
    
    textout += "MultiPeriodDesignSpectrum\n"
    index = len(t)
    j = 0
    while j < index:
        textout += str(t[j])+ ", " + str(s[j])+"\n"
        j+= 1
    j = 0
    textout += "Recommended ELF Design Spectrum\n"
    while j < index:
        textout += str(t[j])+ ", " + str(elfs[j])+"\n"
        j+= 1
    textout += "MultiPeriodMCErSpectrum\n"
    index = len(tmce)
    j = 0
    while j < index:
        textout += str(tmce[j])+ ", " + str(smce[j])+"\n"
        j+= 1
    return(textout)


def mywritefileest(ldata, sitecl, sexp, elfs):
    sitetitle = mysite
    riskct = riskc

    textout = ""
    index = 0
    p = ldata["response"]["data"]
    t = ldata["response"]["data"]["multiPeriodDesignSpectrum"]["periods"]
    s = ldata["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"]
    tmce = ldata["response"]["data"]["multiPeriodMCErSpectrum"]["periods"]
    smce = ldata["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]
    textout += versionstr + "\n"
    textout += "Data source is USGS (ASCE 722 Database) and OpenStreetMaps.\nAuthors do not assume any responsibility or liability for its accuracy.\n"
    textout += "Use of the output of this program does not imply approval by the governing building code bodies responsible for building code approval and interpretation for the building site described by latitude/longitude location.\n"
    textout += "\n \n"
    textout += sitetitle + "\n" + address + "\n"
    textout += "The location is " + str(lat) + ", " + str(longt) + " with Site Class " + sitecl + " and Risk Category "+ riskct + "\n"
    if swv != 0.0 :
        textout += "Site Class based on a shear wave velocity of " + str(swv) + "ft/s\n"
    textout += "pga from design spectra = " + str(round(s[0], 3)) + "\n"
    for key, value in p.items():
        if index <= 11:
            textout += str(key)+ ", " + str(value)+"\n"     
        index += 1
    
    textout += "MultiPeriodDesignSpectrum\n"
    index = len(t)
    j = 0
    while j < index:
        textout += str(t[j])+ ", " + str(s[j])+"\n"
        j+= 1
    j = 0
    textout += "Recommended ELF Design Spectrum\n"
    while j < index:
        textout += str(t[j])+ ", " + str(elfs[j])+"\n"
        j+= 1

    textout += "Interpolated MultiPeriodDesignSpectrum\n"
    index = len(t)
    j = 0
    while j < index:
        textout += str(t[j])+ ", " + str(sexp[j])+"\n"
        j+= 1
    textout += "MultiPeriodMCErSpectrum\n"
    index = len(tmce)
    j = 0
    while j < index:
        textout += str(tmce[j])+ ", " + str(smce[j])+"\n"
        j+= 1
    return(textout)


for k, v in st.session_state.items():
    st.session_state[k] = v
st.subheader(":blue[ASCE7-22 Seismic Parameter Input]")


# st.query_params.from_dict({"address": "elk grove, CA", "title": "Cool location", "long": -120, "lat": 39, "shearwavevelo": 1200})

if "title" in st.query_params:
    inTitle = st.query_params["title"]
    st.session_state['myTitle'] = inTitle

if "address" in st.query_params:
    inAdd = st.query_params["address"]
    st.session_state['myaddress'] = inAdd
elif 'myaddress' in st.session_state:
    inAdd = st.session_state['myaddress']
else:
    inAdd = ""

if "long" in st.query_params:
    inLong = float(st.query_params["long"])
    st.session_state['mylong'] = inLong
elif 'mylong' in st.session_state:
    inLong = st.session_state['mylong']
else:
    inLong = -121.0

if "lat" in st.query_params:
    inLat = float(st.query_params["lat"])
    st.session_state['mylat'] = inLat
elif 'mylat' in st.session_state:
    inLat = st.session_state['mylat']
else:
    inLat = 38.0

RiskCategoryList=["I","II","III","IV"]
if "riskcat" in st.query_params:
    if st.query_params["riskcat"] in RiskCategoryList:
        inRisk = st.query_params["riskcat"]
        st.session_state['myrisk'] = inRisk
    elif 'myrisk' in st.session_state:
        inRisk = st.session_state['myrisk']
    else:
        inRisk = "IV"
else:
    if 'myrisk' in st.session_state:
        inRisk = st.session_state['myrisk']
    else:
        inRisk = "IV"

if "shearwavevelo" in st.query_params:
    inSwv = float(st.query_params["shearwavevelo"])
    st.session_state['myswv'] = inSwv
elif 'myswv' in st.session_state:
    inSwv = st.session_state['myswv']
else:
    inSwv = 0.0
versionstr = "Version 2.0 (revision date 7/17/2026)"
st.badge(versionstr, color="green")
st.write("Data source is USGS (ASCE 722 Database) and OpenStreetMaps.\nAuthors do not assume any responsibility or liability for its accuracy.")
st.write("Use of the output of this program does not imply approval by the governing building code bodies responsible for building code approval and interpretation for the building site described by latitude/longitude location.")
st.divider()
if st.session_state["myTitle"] == "":
    mysite = st.text_input("Title for report", placeholder="Enter title for report", key="title")
else:
    mysite = st.text_input("Title for report", st.session_state["myTitle"] , key="title")
st.session_state['myTitle'] = mysite

st.write("Either enter Shear Wave Velocity or pick Site Class:" )
st.write("(Shear Wave Velocity will be used when entered)")


c1, c2 =st.columns(2)
with c1:
    t1, t2 = st.tabs(["Shear Wave Velocity", "Site Class"])
    with t1:
        swv = st.number_input("Shear Wave Velocity (ft/s)",value = inSwv, step = 100.0, min_value = 0.0, key="swvss")
        st.write("Note: Shear Wave Velocity of 0.0 will use Site Class selection to generate spectra")
        st.session_state['myswv'] = swv
        estimatedswv= persistent_checkbox("Estimated Shear Wave Velocity?" ,key="estswv")
    with t2:
        placeholder = st.empty()
        siteClassList=["A","B","BC","C","CD","D","DE","E", "Default"]
        siteclass = placeholder.selectbox("Site Class",siteClassList,index = 8, key="siteclass")
        if swv != 0:
            st.write("Note: Clear Shear Wave Velocity to 0.0 to use generate spectra via site class")
with c2:

    riskc = st.selectbox("Risk Category",RiskCategoryList, index = RiskCategoryList.index(inRisk))
    st.session_state['myrisk'] = riskc

st.divider()
st.write("Either provide Address or Lat/Long Pair (leave Address blank)")

tab1, tab2 = st.tabs(["Lat/Long", "Address"])

with tab2:
    addressg = st.text_input("Address", inAdd, placeholder="123, streat name, city, CA")
    st.session_state['myaddress'] = addressg

with tab1:
    latitude= st.number_input("Latitude",value=inLat, step = 0.1, min_value = -90.0, max_value= 90.0)
    st.session_state['mylat'] = latitude
    longitude= st.number_input("Longitude",value =inLong, step = 0.1, min_value =-180.0, max_value=180.0)
    st.session_state['mylong'] = longitude
    if addressg != "":
        st.write("Note: Clear address to generate spectra using lat/long pair")

whichgeolocator = st.radio("Select Geocoding Service (will use cache where it exists)", ("OpenStreetMaps (Default)", "ArcGIS"), index=0, key="geoloc")

st.write("Note: Geocoding services have usage limits and may be slow at times. If geocoding fails, try again or switch geocoding service.")
def click_button():
    st.session_state.clicked = True

st.button('Run', on_click=click_button)

# st.write(st.session_state)
if st.session_state.clicked:
    onclick()
    
    st.subheader("Download output file for Spectra")
    sfile= st.checkbox("Save output file")
    if sfile:
        st.download_button("Save output file", textout, file_name="respspectra.txt",)

sds_latex = "S_{DS}"
sd1_latex = "S_{D1}"

if st.session_state.clicked:
    st.subheader("ASCE7-22 Local Variation of Seismic Parameters")
    st.write("Computed for selected site class only,\n Will take some time depending on latency of USGS website,\n Select to start")
    locvart= st.checkbox(f"Check Local Variation of ${sds_latex}$ and ${sd1_latex}$")
    if locvart==1:
        contourf(lat, longt, riskct)

st.divider()


