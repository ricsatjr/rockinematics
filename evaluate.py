import numpy as np
import mplstereonet.mplstereonet as st
import matplotlib.pyplot as plt
import envelopes as env
from itertools import combinations as itercomb

def line_in_cone(cplunge, cbearing, cangle, lplunge, lbearing):
    """
    Evaluates if line is within a small circle (cone), such that the acute 
    angle between the line, and the cone center line should be less than 
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
    ang_dist=np.degrees(st.angular_distance(st.line(cplunge,cbearing),st.line(lplunge,lbearing),False))
    return ang_dist<cangle

def line_above_plane(strike,dip,lplunge,lbearing):
    """
    Evaluates if line is above a plane (great circle), such that 
    1) the line has a smaller plunge compared to that of its vertical projection 
    onto the plane (as rake), and 2) it has the same bearing as that of the rake. 
    
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
    smallerplunge=rplunge>lplunge
#    evaluate if line has the same bearing as the rake
    samebearing=np.abs(rbearing-lbearing)<0.001
    return smallerplunge*samebearing
    

def planarFailure(sstr,sdip,jfriction,jstr,jdip,to_plot=True):
    """
    Evaluates planar failure of joints vis-a-vis a slope face 
    with a given strike and dip, such that a joint's pole plots 1) within the
    planar daylight envelope, and 2) outside the planar friction envelope
    
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
    planarFail=(inDaylight==True)&(outFriction==True)
#    plotting results
    if uniformFriction and to_plot:
        env.setup_axes(sstr,sdip,jfriction[0],failure='planar',to_plot=True)
        plt.gca().pole(jstr[~planarFail],jdip[~planarFail],color='0.5',marker='.')
        plt.gca().pole(jstr[planarFail],jdip[planarFail],color='r',marker='.')
    return planarFail

def wedgeFailure(sstr,sdip,jfriction,jstr,jdip,to_plot=True):
    """
    Evaluates wedge failure of joints vis-a-vis a slope face 
    with a given strike and dip, such that a line defined by the intersection
    of two joints plots 1) on the convex side of the wedge daylight envelope, 
    and 2) within the wedge friction envelope. For each line, it conservatively
    uses the smaller friction angle (jfriction)
    
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
    wedgeFail: boolean array of size = len(np.atleast_1d(jstr)) 
        Indicates if corresponding joints will allow wedge failure.
    """
#    ensure jstr, jdip, and jfriction are 1-d arrays
    jstr,jdip=np.atleast_1d(jstr,jdip)
    try:
        len(jfriction)
        uniformFriction=False
    except:
        jfriction=jfriction*(np.ones(len(jstr)))
        uniformFriction=True
#    get plunge and bearing of unique joint pair intersections
    c=np.array(list(itercomb(range(len(jstr)),2)))
    wl_plunge,wl_bearing=st.plane_intersection(jstr[c[:,0]], jdip[c[:,0]], jstr[c[:,1]], jdip[c[:,1]])
#    get minimum jfriction for each joint pair
    wl_friction=np.min((np.vstack([jfriction[c[:,0]],jfriction[c[:,1]]])),axis=0)
#    determinde daylight and friction envelopes
    wde_strike, wde_dip=env.wedge_daylight(sstr,sdip,False)
    wfe_plunge, wfe_bearing, wfe_angle=env.wedge_friction(wl_friction,False)
#    evaluate if wedge lines are within daylight and friction envelopes (cones)
    convexDaylight=np.empty(len(wl_plunge))
    inFriction=np.empty(len(wl_plunge))
    for a in range(len(wl_plunge)):
        convexDaylight[a]=line_above_plane(wde_strike, wde_dip, wl_plunge[a], wl_bearing[a])
        inFriction[a]=line_in_cone(wfe_plunge[a], wfe_bearing[a], wfe_angle[a], wl_plunge[a], wl_bearing[a])
    wedgeFail=((convexDaylight==True)&(inFriction==True))
#    plotting results
    if uniformFriction and to_plot:
        env.setup_axes(sstr,sdip,jfriction[0],failure='wedge',to_plot=True)
        plt.gca().line(wl_plunge[~wedgeFail],wl_bearing[~wedgeFail],color='0.5',marker='.')
        plt.gca().line(wl_plunge[wedgeFail],wl_bearing[wedgeFail],color='r',marker='.')
    return wedgeFail

def topplingFailure(sstr,sdip,jfriction,jstr,jdip,to_plot=True):
    """
    Evaluates toppling failure of joints vis-a-vis a slope face 
    with a given strike and dip, such that a joint's pole plots 1) within the
    toppling slip limits, and 2) on the convex side of the toppling friction 
    envelope
    
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
    convexFriction=np.empty(len(jstr))
    for a in range(len(jstr)):
        inSlipLimit1[a]=~line_in_cone(tsl_plunge, tsl_bearing, tsl_angle, jplunge[a], jbearing[a])
        inSlipLimit2[a]=~line_in_cone(tsl_plunge, tsl_bearing+180, tsl_angle, jplunge[a], jbearing[a])
        convexFriction[a]=line_above_plane(tfe_strike[0],tfe_dip[a],jplunge[a], jbearing[a])
    topplingFail=((inSlipLimit1==True)&(inSlipLimit2==True)&(convexFriction==True))
#    plotting results    
    if uniformFriction and to_plot:
        env.setup_axes(sstr,sdip,jfriction[a],failure='toppling',to_plot=True)
        plt.gca().pole(jstr[~topplingFail],jdip[~topplingFail],color='0.5',marker='.')
        plt.gca().pole(jstr[topplingFail],jdip[topplingFail],color='r',marker='.')
    return topplingFail

plt.close('all')

sstr=60
sdip=70
jfriction=25    
n=200
jstr=np.random.randint(0,360,n)
jdip=np.random.randint(30,90,n)
planarFail=planarFailure(sstr,sdip,jfriction,jstr,jdip)
topplingFail=topplingFailure(sstr,sdip,jfriction,jstr,jdip)
wedgeFailure(sstr,sdip,jfriction,jstr[:25],jdip[:25])

