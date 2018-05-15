import numpy as np
import mplstereonet.mplstereonet as st
import matplotlib.pyplot as plt
import envelopes as env



def line_in_cone(cplunge, cbearing, cangle, lplunge, lbearing):
    """
    Evaluates if line is within a small circle (cone). Angular difference 
    between the line, and the cone center line should be less than 
    the cone angle. 
    
    Parameters
    ----------
    cplunge, cbearing, cangle: integer or float
        plunge, bearing, and angle defining a cone
    lplunge, lbearing: integer or float
        plunge, bearing of a line
        
    Returns
    -------
    ang_diff<cangle: boolean
        True if line is within the cone 
    """
    ang_diff=np.degrees(st.angular_distance(st.line(cplunge,cbearing),st.line(lplunge,lbearing),False))
    return ang_diff<cangle
    
    

def line_above_plane(strike,dip,lplunge,lbearing):
    """
    Evaluates if line is above a plane (great circle). Line should be 
    above its vertical projection onto the plane (as rake), and should have
    the same bearing as the rake. 
    
    Parameters
    ----------
    strike : int or float
        The strike of the plane in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".
    dip : int or float
        The dip of the plane in degrees.
    lplunge, lbearing: integer or float
        plunge, bearing of a line
        
    Returns
    -------
    aboveplane*samebearing: boolean
        True if line is above the plane
    """
#    project line to plane as rake
    rake_angle=st.azimuth2rake(strike, dip, lbearing)
#    convert rake to line defined by plunge and bearing
    rlon,rlat=st.rake(strike, dip, rake_angle)
    rplunge,rbearing=st.geographic2plunge_bearing(rlon, rlat)
#    evaluate if line is above the rake
    aboveplane=rplunge>lplunge
#    evaluate if line has the same bearing as the rake
    samebearing=np.abs(rbearing-lbearing)<0.001
    return aboveplane*samebearing
    

def planarFailure(sstr,sdip,jfriction,jstr,jdip,to_plot=True):
    """
    Evaluates planar failure of joints vis-a-vis a slope face 
    with a given strike and dip.
    
    Parameters
    ----------
    sstr : int or float
        The strike of the slope face in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".
    sdip : int or float
        The dip of the slope face in degrees.
    jfriction : int, float, or array of int or float
        The friction angle of the joint plane in degrees.
    jstr : int, float, or array of int or float
        The strike of the joint plane in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".
    jdip : int, float, or array of int or float
        The dip of the joint plane in degrees.
        
    Returns
    -------
    planarFail: boolean array of size = len(np.atleast_1d(jstr)) 
        Indicates if corresponding joints will allow planar failure.
    """
#    ensure jstr, jdip, and jfriction are 1-d arrays
    jstr,jdip=np.atleast_1d(jstr,jdip)
    try:
        len(jfriction)
        uniformFriction=False
    except:
        jfriction=jfriction*(np.ones(len(jstr)))
        uniformFriction=True
#    determinde daylight and friction envelopes
    pde_plunge, pde_bearing, pde_angle=env.planar_daylight(sstr,sdip,False)
    pfe_plunge, pfe_bearing, pfe_angle=env.planar_friction(jfriction,False)
#    convert joint plane (strike-dip) to pole (plunge-bearing)
    jplunge,jbearing=st.pole2plunge_bearing(jstr,jdip)
#    evaluate if joint poles are contained within daylight and friction envelopes (cones)
    inDaylight=np.empty(len(jstr))
    outFriction=np.empty(len(jstr))
    for a in range(len(jstr)):
        inDaylight[a]=line_in_cone(pde_plunge, pde_bearing, pde_angle, jplunge[a], jbearing[a])
        outFriction[a]=~line_in_cone(pfe_plunge[a], pfe_bearing[a], pfe_angle[a], jplunge[a], jbearing[a])
    
    planarFail=((inDaylight==True)&(outFriction==True))
    
    if uniformFriction and to_plot:
        env.setup_axes(sstr,sdip,jfriction[a],failure='planar',to_plot=True)
        plt.gca().pole(jstr[~planarFail],jdip[~planarFail],color='0.5',marker='.')
        plt.gca().pole(jstr[planarFail],jdip[planarFail],color='r',marker='.')
    return planarFail



def topplingFailure(sstr,sdip,jfriction,jstr,jdip,to_plot=True):
    """
    Evaluates toppling failure of joints vis-a-vis a slope face 
    with a given strike and dip.
    
    Parameters
    ----------
    sstr : int or float
        The strike of the slope face in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".
    sdip : int or float
        The dip of the slope face in degrees.
    jfriction : int, float, or array of int or float
        The friction angle of the joint plane in degrees.
    jstr : int, float, or array of int or float
        The strike of the joint plane in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".
    jdip : int, float, or array of int or float
        The dip of the joint plane in degrees.
        
    Returns
    -------
    topplingFail: boolean array of size = len(np.atleast_1d(jstr)) 
        Indicates if corresponding joints will allow toppling failure.
    """
#    ensure jstr, jdip, and jfriction are 1-d arrays
    jstr,jdip=np.atleast_1d(jstr,jdip)
    try:
        len(jfriction)
        uniformFriction=False
    except:
        jfriction=jfriction*(np.ones(len(jstr)))
        uniformFriction=True
#    determine daylight and friction envelopes
    tsl_plunge,tsl_bearing,tsl_angle=env.toppling_slipLimits(sstr,sdip,False)
    tfe_strike, tfe_dip=env.toppling_friction(sstr,sdip,jfriction,False)
#    convert joint plane (strike-dip) to pole (plunge-bearing)
    jplunge,jbearing=st.pole2plunge_bearing(jstr,jdip)
#    evaluate if joint poles are contained within slip limits (cones) and friction envelope (great circles)
    inSlipLimit1=np.empty(len(jstr))
    inSlipLimit2=np.empty(len(jstr))
    aboveFriction=np.empty(len(jstr))
    for a in range(len(jstr)):
        inSlipLimit1[a]=~line_in_cone(tsl_plunge, tsl_bearing, tsl_angle, jplunge[a], jbearing[a])
        inSlipLimit2[a]=~line_in_cone(tsl_plunge, tsl_bearing+180, tsl_angle, jplunge[a], jbearing[a])
        aboveFriction[a]=line_above_plane(tfe_strike[0],tfe_dip[a],jplunge[a], jbearing[a])
        
    topplingFail=((inSlipLimit1==True)&(inSlipLimit2==True)&(aboveFriction==True))
    
    if uniformFriction and to_plot:
        env.setup_axes(sstr,sdip,jfriction[a],failure='toppling',to_plot=True)
        plt.gca().pole(jstr[~topplingFail],jdip[~topplingFail],color='0.5',marker='.')
        plt.gca().pole(jstr[topplingFail],jdip[topplingFail],color='r',marker='.')
    return topplingFail

plt.close('all')

sstr=0
sdip=70
jfriction=25    
n=200
jstr=np.random.randint(0,360,n)
jdip=np.random.randint(0,90,n)
planarFail=planarFailure(sstr,sdip,jfriction,jstr,jdip)
topplingFail=topplingFailure(sstr,sdip,jfriction,jstr,jdip)


