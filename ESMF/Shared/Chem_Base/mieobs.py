"""
   Implements Python interface to Mie Calculator. Based on f2py extension Mie_.

"""

from numpy import array, isfortran, ones, float32, size, zeros
from pylab import plot, axis, title
from MAPL  import config

VNAMES = [ 'DU001', 'DU002', 'DU003', 'DU004', 'DU005',
           'SS001', 'SS002', 'SS003', 'SS004', 'SS005',
           'BCPHOBIC', 'BCPHILIC',
           'OCPHOBIC', 'OCPHILIC',
           'SO4' ]

import MieObs_ 

def getMieDims(rcfile='Aod_EOS.rc'):
    """
        Return dimensions of Mie-table like nPol and nMom.
    """
    cf = config.Config(rcfile)
    dutable = cf('filename_optical_properties_DU')
    nCh, nRh, nBin, nMom, nPol, rc = MieObs_.getmiedims(dutable) 
    if rc != 0:
       raise ValueError, "on return from getMieDims, rc = %d"%rc

    return (nMom, nPol)

#---
def _toArray(channels):
    """Return numpy array, if *channels* is a scalar."""

    if type(channels) is int:
        channels = (channels, )
    elif type(channels) is float:
        channels = (channels, )
    channels = array(channels)
    return channels
#---
def aerToUpper(aer):
    """
    Create upper case aliases for aer variables to cope with changes
    in filename.
    """
    vnames = aer.__dict__.keys()
    for v in vnames:
        V = v.upper()
        if v != V:
            aer.__dict__[V] = aer.__dict__[v] # make alias
#---
def getEdgeVars(aer,ptop=1.):
    """
    Given aer object with (airdens,delp) attributes
    returns
           pe --- layer edge pressure [Pa]
           ze --- layer edge height above sfc [m]
           te --- temperature [K] at layer edge

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    """
    if needs_transpose(aer):
        pe, ze, te = MieObs_.getedgevars(aer.AIRDENS.T,aer.DELP.T,ptop)
    else: 
        pe, ze, te = MieObs_.getedgevars(aer.AIRDENS,aer.DELP,ptop)

    return (pe,ze,te)

#--
def getAOPscalar(aer,channels,vnames=VNAMES,vtypes=None,Verbose=False,rcfile='Aod_EOS.rc'):
    """
    Compute (tau,ssa,g) given aer object.

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    """

    # Variable type for mie calculation
    # ---------------------------------
    if vtypes is None:
        vtypes = vname     # to be used in mie calculation

    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)

    # Pack inputs
    # -----------
    nq, nch = len(vnames), len(channels)
    nobs = size(aer.PS)
    if needs_transpose(aer): # aer is (nobs,km)
        km = aer.DELP.shape[1]
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH.T
#        rh = zeros((km,nobs))
        for n, v in zip(range(nq),vnames):
            V = v.upper()
            qm[:,n,:] = aer.__dict__[V].T * aer.DELP.T / 9.81 
    else:                    # aer is (km,nobs)
        km = aer.DELP.shape[0]
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH
#        rh = zeros((km,nobs))
        for n, v in zip(range(nq),vnames):
            V = v.upper()
            qm[:,n,:] = aer.__dict__[V] * aer.DELP / 9.81 

    # Do the Mie calculation
    # ----------------------
    tau,ssa,g,rc = MieObs_.getaopscalar(rcfile,channels,pad(vtypes),Verbose,qm,rh)

    if rc!=0:
        print "<<<ERROR>>> on return from MieObs_.getaopscalar, rc = ", rc
        raise ValueError, 'cannot get Aerosol Optical Properties (scalar version)'

    return (tau,ssa,g)
#--
def getAOPvector(aer,channels,I=None,vnames=VNAMES,vtypes=None,
                 Verbose=False,rcfile='Aod_EOS.rc',nMom=301):
    """
    Compute (tau,ssa,g,pmom) given aer object.

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    
    J --- index of subset of observations to process
    """

    # Variable type for mie calculation
    # ---------------------------------
    if vtypes is None:
        vtypes = vname     # to be used in mie calculation

    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)

    # Pack inputs
    # -----------
    nq, nch = len(vnames), len(channels) 
    if I is None:
        nobs = size(aer.PS)
        I = range(0,nobs)
    else :
        nobs = len(I)
    
    if needs_transpose(aer): # aer is (nobs,km)
        km = aer.DELP.shape[1]
        
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH[I].T
        for n, v in zip(range(nq),vnames):
            V = v.upper()
            qm[:,n,:] = aer.__dict__[V][I].T * aer.DELP[I].T / 9.81 
    else:                    # aer is (km,nobs)
        km = aer.DELP.shape[0]
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH[:,I]
        for n, v in zip(range(nq),vnames):
            V = v.upper()
            qm[:,n,:] = aer.__dict__[V][:,I] * aer.DELP[:,I] / 9.81 
    

    # Do the Mie calculation
    # ----------------------
    nMom_mieTable,nPol_mieTable = getMieDims(rcfile=rcfile)  # return nMom & nPol of the Mie tables
    nPol = 6                                                 # for dust non spherical
    tau,ssa,g,pmom,rc = MieObs_.getaopvector(rcfile,channels,pad(vtypes),Verbose,qm,rh,nMom,nPol)

    if rc!=0:
        print "<<<ERROR>>> on return from MieObs_.getaopvector, rc = ", rc
        raise ValueError, 'cannot get Aerosol Optical Properties (vector version)'

    return (tau,ssa,g,pmom)
#---
def getAOPext(aer,channels,I=None,vnames=VNAMES,vtypes=None,Verbose=False,rcfile='Aod3d_532nm.rc'):
    """
    Compute (ext,backscat,aback_sfc,aback_toa) given aer object.

    Input arrays can be (nobs,km) or (km,nobs).
    It always returns arrays that are (km,nobs).
    
    I --- index of subset of observations to process
    """

    # Variable type for mie calculation
    # ---------------------------------
    if vtypes is None:
        vtypes = vname     # to be used in mie calculation

    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)

    # Pack inputs
    # -----------
    nq, nch = len(vnames), len(channels) 

    if I is None:
        nobs = size(aer.PS)
        I = range(0,nobs)
    else :
        nobs = len(I)
    
    if needs_transpose(aer): # aer is (nobs,km)
        km = aer.DELP.shape[1]
        
        qc = ones((km,nq,nobs),dtype=float32)
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH[I].T
#        rh = zeros((km,nobs))
        for n, v in zip(range(nq),vnames):
            V = v.upper()
            qc[:,n,:] = aer.__dict__[V][I].T * aer.AIRDENS[I].T 
            qm[:,n,:] = aer.__dict__[V][I].T * aer.DELP[I].T / 9.81 
    else:                    # aer is (km,nobs)
        km = aer.DELP.shape[0]
        
        qc = ones((km,nq,nobs),dtype=float32)
        qm = ones((km,nq,nobs),dtype=float32)
        rh = aer.RH[I]
#        rh = zeros((km,nobs))
        for n, v in zip(range(nq),vnames):
            V = v.upper()
            qc[:,n,:] = aer.__dict__[V][:,I] * aer.AIRDENS[:,I] 
            qm[:,n,:] = aer.__dict__[V][:,I] * aer.DELP[:,I] / 9.81 

    # Do the Mie calculation
    # ----------------------
    ext,sca,backscat,aback_sfc,aback_toa,depol,rc = \
        MieObs_.getext(rcfile,channels,pad(vtypes),Verbose,qc,qm,rh)

    if rc!=0:
        print "<<<ERROR>>> on return from MieObs_.getaopvector, rc = ", rc
        raise ValueError, 'cannot get Aerosol Optical Properties (vector version)'

    return (ext,sca,backscat,aback_sfc,aback_toa,depol)

#........................................................................
def pad(names):
    """
    Make all strings in list *names* the same size for f2py's benefit.
    """
    return [ "%-16s"%v for v in names ]

#........................................................................
def needs_transpose(aer):
    """
    Returns
         True  if aer object has shapes (nobs,km)
         False if aer object has shapes (km,nobs)
    This is needed to cope with the way the fortran expects
    the data arrays.
    """
    if aer.PS.shape[0] == aer.DELP.shape[0]: # (nobs,km)
        return True
    else:
        return False
    
#........................................................................

if __name__ == "__main__":

    from pyobs import LIDAR_L2, NPZ
    from mie   import getTopo
    
    channels = array([532.,])

    c_dn = '/nobackup/2/vbuchard/LIDAR/dR_Fortuna-2-4-b4/' # Pete's Calipso interp
    a_dn = '/nobackup/2/vbuchard/LIDAR/dR_Fortuna-2-4-b4/'         # aer_Nv + interp

    c_fn = 'dR_Fortuna-2-4-b4.calipso_532nm.20090715.nc' # Pete's interp
    a_fn = 'dR_Fortuna-2-4-b4.calipso_aer.20090715.npz' # aer collocation

    # Read relevant data
    # ------------------
    c = LIDAR_L2(c_dn+c_fn)  # Pete's interp with CALPSO coordinates
    a = NPZ(a_dn+a_fn)        # aer_v interpolated to obs location
    aerToUpper(a)
#   ---------------------------------------------------------------------------    
#   NOTE:
#
#   For creating the npz file above with the aer_v data interpolated to obs ,
#   location you do this:
#      c.sampleFile(aer_v,npzFile=npzFile)   
#   where aer_v is the full, gridded, netcdf file.
#   In principle, you do not need to write the npzFile as the result of the
#   sampling is also available as attribute *sample*. So, instead of reading
#   the npzFile you could do
#     a = c.sample
#   But sampling takes time, so saving a "sampling  file" is convenient.
#   BTW, npzFile is too python specific. Soon, I'll implement a NetCDF option,
#   so you will be able to say
#     c.sample(aer_v,ncFile=NetCDF_filename)
#   See bottom of GMAO_pyobs/pyobs.lidar_l2.py file for an example.
#   ---------------------------------------------------------------------------    

    pe, ze, te = getEdgeVars(a)
    
    tau, ssa, g = getAOPscalar(a,channels,rcfile='Aod_EOS.rc')
    ext,sca,backscat,aback_sfc,aback_toa = getAOPext(a,channels)
