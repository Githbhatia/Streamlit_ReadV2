import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd
import math

def persistent_toggle(label, key):
    state = st.toggle(label, value=st.session_state.checklist_items.get(key, False), key=key)
    st.session_state.checklist_items[key] = state
    return state

def getHf(zhratio):
    a1 = min(1/tA,2.5)
    a2 = max((1-(0.4/tA)**2),0.0)
    hF = 1+ a1*zhratio + a2*zhratio**10 
    # print ("Hf = " + str(hF))   
    return(hF)

def getaltHf(zhratio):
    hF = 1+ 2.5*zhratio 
    # print ("Hf = " + str(hF))   
    return(hF)

for k, v in st.session_state.items():
    st.session_state[k] = v
sds  = st.session_state['sds']
sds_latex = "S_{DS}"
sd1_latex = "S_{D1}"
st.subheader(":blue[ASCE7-22 Fp Calculation]")
versionstr = "Version 2.0 (revision date 7/17/2026)"
st.badge(versionstr, color="green")
st.write("USING THE DEFAULT OPTIONS WILL LEAD TO CONSERVATIVE RESULTS")
if not st.session_state.clicked:
    st.write (f"Please click on :red[RUN] button in the previous page to get Sds for this site")
    st.write (f"To proceed anyways, manually enter:red[${sds_latex}$]")

FP="F_{p}"
if st.session_state['myTitle'] == "":
    mysite = st.text_input("Title for report", placeholder="Enter title for report", key="title")
else:
    mysite = st.text_input("Title for report",st.session_state['myTitle'], key="title")
st.session_state['myTitle'] = st.session_state['title']
if sds <= 0.0:
    sds = st.number_input(f"Enter ${sds_latex}$", value= sds, format="%0.3f", min_value=0.0)
else:
    sds = st.number_input(f"${sds_latex}$, as obtained in previous page, can modify here", value= sds, format="%0.3f")
if sds <= 0.0:
    st.write(f":red[Enter a valid value for ${sds_latex}$]")
    st.stop()

df = pd.read_csv('ASCE722Ch13.csv')
df.set_index('Menuitems', inplace=True)
st.divider()

OPD = st.checkbox("Select for Fp for OPDs", value=True, key="check1")
if OPD:
    zOPD = [0.0, 1/3, 2/3, 1.0]
    nlist = ["Interior nonstructural walls and partitions;  Light frame > 9 ft (2.74 m) in height","Ceilings; All"]
    iP = 1.5
    rU = 1.3
    fPMax = 1.6*sds*iP
    fPMin = 0.3*sds*iP
    fPOPDalt=[ [0]*len(zOPD) for i in range(len(nlist))]
    omegaOPD= []
    for j in range(len(nlist)):
        selecteditem = nlist[j] 
        car0 = df.loc[selecteditem].values[0]
        car1 = df.loc[selecteditem].values[1]
        rPO = df.loc[selecteditem].values[2]
        omegaOPD.append(df.loc[selecteditem].values[3])
        for i in range(len(zOPD)):
            hFLalt = getaltHf(zOPD[i])
            if zOPD[i] == 0.0:
                fPOPDalt[j][i] = min(max(0.4*sds*iP*(hFLalt/1.0)*(car0/rPO),fPMin),fPMax)
            else:
                fPOPDalt[j][i] = min(max(0.4*sds*iP*(hFLalt/rU)*(car1/rPO),fPMin),fPMax)
    Type1 ='Fp/Wp for Partition Walls (OmegaOP =' + str(round(omegaOPD[0],2)) + ')'
    Type2 ='Fp/Wp for Ceilings (OmegaOP =' + str(round(omegaOPD[1],2)) + ')'
    st.write(":red[$F_{p}$ for OPDs:]")
    st.write("Assumes $I_{p}$ = 1.5, $R_{\\mu}$= 1.3, $H_{f}$ per ASCE 7-22 Eq 13.3-5")
    loc = ["At Grade","Lower Third","Middle Third","Upper Third (including roof)"]

    st.write("Governing $F_{p}$:")
    dfsfP=pd.DataFrame.from_dict({'Location, (Sds = ' + str(round(sds,3)) + ')': loc, Type1: fPOPDalt[0], Type2: fPOPDalt[1]})
    st.dataframe(dfsfP, hide_index=True)
st.divider()    


DfP = st.checkbox("Select for Detailed Fp Calulations for any nonstructural component", value=False, key="check2")
if DfP:
        
    st.write(":red[Detailed $F_{p}$ Calculations:]")
    if st.session_state.selecteditem != "":
        _selecteditem22 = st.session_state.selecteditem
        selecteditem22 = st.selectbox("Select Nonstructural item (ASCE 7-22 Tables 13.5-1 and 13.6-1)",df.index, index = list(df.index).index(_selecteditem22), key="nonstructural")
    else:
        selecteditem22= st.selectbox("Select Nonstructural item (ASCE 7-22 Tables 13.5-1 and 13.6-1)",df.index, index = 1, key="nonstructural")

    st.session_state.selecteditem22 = st.session_state.nonstructural

    car0 = df.loc[selecteditem22].values[0]
    car1 = df.loc[selecteditem22].values[1]
    rPO = df.loc[selecteditem22].values[2]
    omegaOP = df.loc[selecteditem22].values[3]

    sc1,sc2 =st.columns(2)
    with sc1:
        I_p = "I_{p}"
        if st.session_state.selectedIp != 0.0:
            iP = float(st.selectbox(f"${I_p}$, Component Importance Factor",(1.0,1.5), index = list((1.0,1.5)).index(st.session_state.selectedIp), key="Ip"))   
        else:
            iP = float(st.selectbox(f"${I_p}$, Component Importance Factor",(1.0,1.5), index = 1, key="Ip"))
        st.session_state.selectedIp = iP
    

    # sc3,sc4 =st.columns(2)
    # with sc3:
    #     Z = "Z"
    #     # z = st.number_input(f"${Z}$, height above base",value= 90.0)
    #     if st.session_state.UserZvalues != "":
    #         zStr = st.text_input(f"${Z}$, height above base (multiple ok,separate with commas)",value= st.session_state.UserZvalues,  key="zvalues")

    #     else:
    #         zStr = st.text_input(f"${Z}$, height above base (multiple ok,separate with commas)",str("0, 15, 30, 45, 60, 75, 90, 100"),key="zvalues")

    # st.session_state.UserZvalues = st.session_state.zvalues


    # if st.session_state.UserZlabels != "":
    #     zLbl = st.text_input("Labels corresponding to Z values (Separate with commas,Optional)",st.session_state.UserZlabels, key="zLables")
    # else:
    #     zLbl = st.text_input("Labels corresponding to Z values (Separate with commas,Optional)",str("Grnd Level, Level 2, Level 3, Level 4, Level 5, Level 6, Mech Level, Roof"),key="zLables")
    # zLblist = [i.strip() for i in zLbl.split(",")]
    # st.session_state.UserZlabels = st.session_state.zLables

    # try:
    #     z =[float(i) for i in zStr.split(",")]
    # except ValueError:
    #     st.write(":red[Invalid input, Please enter numbers separated by commas]")
    #     st.stop()

    # if len(z) > len(zLblist):
    #     for i in range(len(z)-len(zLblist)):
    #         zLblist.append("")
    # if len(z) < len(zLblist):
    #     for i in range(len(zLblist)-len(z)):
    #         zLblist.pop()

    z = [0,  100]
    zLblist = ["Grnd Level","Roof"]
   
    dftable=pd.DataFrame({"Location" :zLblist,"z":z})
    st.write("Enter z values and corresponding labels in the table below (can add/remove rows as needed):")
    st.write("Maximum z value will be taken as h (average roof height), and minimum z value should be >= 0.0 (grade plane)." )
    st.write("List can be out of order. \
    Navigating to another page will reset the table to default values")
    edited_df = st.data_editor(dftable, num_rows="dynamic",column_config={"Location": st.column_config.TextColumn("Location", help="Labels for Z values"), 
        "z": st.column_config.NumberColumn(
            "z (ft)",
            help="Height above grade plane, should be >=0.0",
            min_value=0.0 ),
    }, hide_index=True)
    edited_df = edited_df.dropna( subset=["z"])
    edited_df.fillna({"Location":""}, inplace = True)
    edited_df = edited_df.sort_values(by="z")
    
    z = edited_df["z"].tolist()
    h = max(z)
    zLblist = edited_df["Location"].tolist()


    # with sc4:
    #     H = "H"
    #     if st.session_state["UserHvalues"] != 0.0:
    #         h = st.number_input(f"${H}$, Average roof height of structure in ft",value= st.session_state["UserHvalues"], key="H")
    #     else:
    #         h = st.number_input(f"${H}$, Average roof height of structure in ft",value= 100.0, key="H")
    #     st.session_state.UserHvalues = st.session_state.H
    #     if h < max(z):
    #         st.write(":red[H is < highest value of z, Please correct]")
    #         st.stop()

    st.divider()
    knownstsys = persistent_toggle("Structural System Selection (Unknown system assumed if not enabled)", key="structuralselect")
    if knownstsys:
        dfs = pd.read_csv('ASCE722StructuralSystems.csv')
        dfs.set_index('StructuralSystem', inplace=True)
        if st.session_state.selecteditemStructSys != "":
            _selecteditem = st.session_state.selecteditemStructSys
            selecteditem = st.selectbox("Select Structural System of the Building (ASCE 7-22 Table 12.2-1):",dfs.index, index = list(dfs.index).index(_selecteditem), key="structural") 
            st.session_state.selecteditemStructSys = st.session_state.structural
        else:   
            selecteditem = st.selectbox("Select Structural System of the Building (ASCE 7-22 Table 12.2-1):",dfs.index, index = 49, key="structural")
            st.session_state.selecteditemStructSys = st.session_state.structural
        
        
        r = dfs.loc[selecteditem].values[0]
        oM = dfs.loc[selecteditem].values[1]
        cTs = dfs.loc[selecteditem].values[3]
        xs = dfs.loc[selecteditem].values[4]
    I_e = "I_{e}"
    if st.session_state["selectedIe"] != 0.0:
        ie = float(st.selectbox(f"${I_e}$, Importance Factor for Building",(1.0,1.25,1.5), index = list((1.0,1.25,1.5)).index(st.session_state["selectedIe"]), key="Ie"))
    else:
        ie = float(st.selectbox(f"${I_e}$, Importance Factor for Building",(1.0,1.25,1.5), index = 2, key="Ie"))
    st.session_state["selectedIe"] = st.session_state.Ie
    st.write(f"Selected Structural System: {selecteditem}")
    if knownstsys:
        c1,c2 = st.columns(2)
        with c1:
            R= "R"
            st.write(f"${R}$, Response modification Value = "+ str(round(r,2)))
        with c2:
            Om = "\\Omega_{o}"
            st.write(f"${Om}$ = " + str(round(oM,2)))

        rU = max((1.1*(r/(ie*oM)))**0.5, 1.3)
    else:
        rU = 1.3
    st.write(":red[ASCE 7-22 Equation 13.3-6:]")
    st.latex(r'''\color{red} R_{\mu} = \left[ \frac{1.1 R}{I_{e}\Omega_{o}} \right]^{1/2} \ge 1.3''')
    Ru = "R_{\\mu}"
    st.write(f"${Ru}$ = " +str(round(rU,3)) + " (1.0 used for z = 0.0 per ASCE 7-22 13.3.1.2)")

    st.divider()
    knownperiod = persistent_toggle("Period Known (if not enabled, period is calculated based on Height h)", key="periodselect")
    # st.write(st.session_state)
    if knownperiod:
        Ta = "T_{a}"
        if st.session_state.selecteditemTa != 0.0:
            tA = st.number_input(f"${Ta}$, Lowest fundamental period of structure:",value= st.session_state.selecteditemTa, key="Ta")
        else:
            tA = st.number_input(f"${Ta}$, Lowest fundamental period of structure:",value= 0.5, key="Ta")
        if tA <= 0:
            st.write("tA cannot be zero, revise")
            st.stop()
        st.session_state.selecteditemTa = st.session_state.Ta
    else:
        if knownstsys:
            st.write(f"Per ASCE 7-22 Table 12.2-1 for {selecteditem}:" )
        else:
            cTs = 0.02
            xs = 0.75
            st.write(f"Per ASCE 7-22 Eq 12.8-8 (for \"all other structural systems\"):" )
        tA = cTs*h**xs 
        Ta = "T_{a}"
        Ct= "C_{t}"
        Hx = "h^{x}"
        X= "x"
        st.write(f"${Ct}$ = " +str(round(cTs,3)))
        st.write(f"${X}$ = " +str(round(xs,3)))
        st.write(f"${Ta}$ = ${Ct}$ ${Hx}$ = "  +str(round(tA,3))+ " secs")
    st.divider()




    zhlist = np.concatenate((np.array([0.0,0.001]),np.arange(0.002, 1.001, 0.001)),axis=0)

    zh = [None]*len(z);hF = [None]*len(z);fP = [None]*len(z)
    fPMax = 1.6*sds*iP
    fPMin = 0.3*sds*iP
    fPMaxstr = "1.6 S_{DS} I_p W_p"
    fPMinstr = "0.3 S_{DS} I_p W_p"
    for i in range(len(z)):
        zh[i] =z[i]/h
        hF[i] = getHf(zh[i])
    # Hf = "H_{f}"
    # st.write(f"${Hf}$ = " +str(round(hF,3)))
        if z[i] == 0.0:
            fP[i] = 0.4*sds*iP*(hF[i]/1.0)*(car0/rPO)
        else:
            fP[i] = 0.4*sds*iP*(hF[i]/rU)*(car1/rPO)

        fP[i] = min(max(fP[i],fPMin),fPMax)


    Wp = "W_{p}"
    s1, s2 = st.columns(2)
    with s1:
        st.write(":red[ASCE 7-22 Equation 13.3-4:]")
        st.latex(r'''\color{red} H_{f} = 1 + a_{1} \left(\frac{z}{h} \right) + a_{2} \left(\frac{z}{h} \right)^{10}''')
        st.latex(r'''\color{red} a_{1} = 1/T_{a} \leq 2.5''')
        st.latex(r'''\color{red} a_{2} = [1-(0.4/T_{a})^{2} \geq 0''')

    with s2:
        st.write(":blue[ASCE 7-22 Equation 13.3-4:]   \n (Conservative for periods > 0.4 secs)")
        st.latex(r'''\color{blue} H_{f} = 1 + 2.5 \left(\frac{z}{h} \right) ''')
    st.write(":red[ASCE 7-22 Equation 13.3-1:]")
    st.latex(r'''\color{red} F_{p} = 0.4 S_{DS} I_{p} \left( \frac{H_{f}}{R_{\mu}} \right) \left( \frac{C_{AR}}{R_{po}} \right) W_{p}''')
    c1,c2 = st.columns(2)
    with c1:
        tfPmin=str(round(fPMin,3))
        st.write(f":red[Minimum ${FP}$ = ${fPMinstr}$ = {tfPmin} ${Wp}$]")
    with c2:
        tfmax=str(round(fPMax,3))
        st.write(f":red[Maximum ${FP}$ = ${fPMaxstr}$ = {tfmax} ${Wp}$]")
    c1,c2 = st.columns(2)
    with c1:
        CAR0 = "C_{AR}"
        if math.isnan(car0):
            st.write(f"${CAR0}$ supported at or below grade plane = N/A" )
        else:
            st.write(f"${CAR0}$ supported at or below grade plane = " + str(car0))
    with c2:
        CAR1 = "C_{AR}"
        st.write(f"${CAR1}$ above grade plane,supported by structure = " + str(car1))   
    Rpo = "R_{PO}"
    st.write(f"${Rpo}$ = " + str(rPO))
    Omop = "\\Omega_{op}"
    st.write(f" ${Omop}$ to be used for concrete or masonry post-installed anchors = " + str(round(omegaOP,3)) )

    fPlist = []; fPlistalt=[]
    for i in range(len(zhlist)):
        hFL = getHf(zhlist[i])
        hFLalt = getaltHf(zhlist[i])
        if zhlist[i] == 0.0:
            fPlist.append(min(max(0.4*sds*iP*(hFL/1.0)*(car0/rPO),fPMin),fPMax))
            fPlistalt.append(min(max(0.4*sds*iP*(hFLalt/1.0)*(car0/rPO),fPMin),fPMax))
        else:
            fPlist.append(min(max(0.4*sds*iP*(hFL/rU)*(car1/rPO),fPMin),fPMax))
            fPlistalt.append(min(max(0.4*sds*iP*(hFLalt/rU)*(car1/rPO),fPMin),fPMax))
            
    st.write(f":blue[Governing ${FP}$:]")
    dfsfP=pd.DataFrame({"Location, (Sds = " + str(round(sds,3)) + ")" :zLblist,"z":z,"z/h": zh,"Hf": hF, "Fp/Wp": fP})
    st.dataframe(dfsfP, hide_index=True)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.plot(fPlist, zhlist, label="Calculated Fp using Hf per Eq 13.3-4", color='Red', linewidth=1.0)
    ax.plot(fPlistalt, zhlist, label="Calculated Fp using Hf per Eq 13.3-5", color='Blue', linestyle='--',linewidth=1.0)
    ax.legend()
    for i in range(len(z)):
        ax.plot(fP[i], zh[i], marker='o', label="Governing Fp", color='Black', linestyle='--', linewidth=2.0)
        axmin,axmax = ax.get_xlim()
        arrowlength = (axmax - axmin)/20
        if z[i] >0.85* h:
            ax.annotate(f"{round(fP[i],3)} at " + str(z[i]) + " ("+ zLblist[i] + ")", ha = 'right', xy=(fP[i], zh[i]), xytext=(fP[i]-arrowlength, zh[i]+0.005), arrowprops=dict(facecolor='black', shrink=0.05))
        else:       
            ax.annotate(f"{round(fP[i],3)} at " + str(z[i]) + " ("+ zLblist[i] + ")", xy=(fP[i], zh[i]), xytext=(fP[i]+arrowlength, zh[i]+0.005), arrowprops=dict(facecolor='black', shrink=0.05))
    ax.grid()
    ax.set_xlabel("Fp/Wp")
    ax.set_ylabel("z/h")
    ax.set_title("Variation of Fp with z/h")
    info = (mysite[:100] + '..') if len(mysite) > 100 else mysite
    ax.text(0.99, 0.08, info, horizontalalignment='right', verticalalignment='top', fontsize=10, color ='Black',transform=ax.transAxes)
    ax.text(0.99, 0.05, "Sds = "+str(round(sds,3)), horizontalalignment='right', verticalalignment='top', fontsize=10, color ='Black',transform=ax.transAxes)
    info = (st.session_state['nonstructural'][:150] + '..') if len(st.session_state['nonstructural']) > 150 else st.session_state['nonstructural']
    ax.text(0.99, 0.02, info, horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
    st.pyplot(fig)
    st.divider()

    st.subheader(f":blue[Equivalent ${sds_latex}$ Calculation for OPDs/OPMs approved for previous versions of ASCE 7:]")
    # df16 = pd.read_csv('ASCE716Ch13.csv')
    # df16.set_index('Menuitems', inplace=True)

    # if st.session_state.selecteditem16 != "":
    #     _selecteditem16 = st.session_state.selecteditem16
    #     selecteditem16 = st.selectbox("Select Nonstructural item (ASCE 7-16 Tables 13.5-1 and 13.6-1)",df16.index, index = list(df16.index).index(_selecteditem16), key="nonstructural16")
    # else:
    #     selecteditem16 = st.selectbox("Select Nonstructural item (ASCE 7-16 Tables 13.5-1 and 13.6-1)",df16.index, index = 1, key="nonstructural16")

    # st.session_state.selecteditem16 = st.session_state.nonstructural16

    # ap = df16.loc[selecteditem16].values[0]
    # rp = df16.loc[selecteditem16].values[1]
    st.write("Equivalent ASCE7-16 Chapter 13 item  = " + df.loc[selecteditem22].values[4] )
    ap = df.loc[selecteditem22].values[5]
    rp = df.loc[selecteditem22].values[6]

    ccc = st.columns(2)
    with ccc[0]:
        apstr = "a_{p}"
        st.write(f"${apstr}$ = " + str(ap))
    with ccc[1]:
        rpstr = "R_{p}"
        st.write(f"${rpstr}$ = " + str(rp))
    equivalentSdsList = []
    for i in range(len(z)):
        equivalentSdsList.append((fP[i]/(0.4*ap))*(rp/iP)/3.0)
    
    st.write(f":red[Equivalent ${sds_latex}$ :]")
    dfssds=pd.DataFrame({"Location" :zLblist,"z":z,"z/h": zh,"Equivalent Sds for OPD/OPM": equivalentSdsList})
    st.dataframe(dfssds, hide_index=True)

    st.write(":red[Note: All older OPDs/OPMs are approved for z/h = 1.0]")

st.divider()



