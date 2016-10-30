import numpy as np
from galpy.util import bovy_coords
# Position of stars:


def calc_actions(ra_deg, dec_deg, d_kpc, pm_ra_masyr, pm_dec_masyr,
                 v_los_kms):
    ra_rad = ra_deg * (np.pi / 180.)  # RA [rad]
    dec_rad = dec_deg * (np.pi / 180.)  # dec [rad]

    # Galactocentric position of the Sun:
    X_gc_sun_kpc = 8.  # [kpc]
    Z_gc_sun_kpc = 0.025  # [kpc]

    # Galactocentric velocity of the Sun:
    vX_gc_sun_kms = -9.58  # = -U              [kms]
    vY_gc_sun_kms = 10.52 + 220.  # = V+v_circ(R_Sun) [kms]
    vZ_gc_sun_kms = 7.01  # = W               [kms]

    # a. convert spatial coordinates (ra,dec,d) to (R,z,phi)

    # (ra,dec) --> Galactic coordinates (l,b):
    lb = bovy_coords.radec_to_lb(ra_rad, dec_rad, degree=False, epoch=2000.0)
    l_rad = lb[:, 0]
    b_rad = lb[:, 1]

    # (l,b,d) --> Galactocentric cartesian coordinates (x,y,z):
    xyz = bovy_coords.lbd_to_XYZ(l_rad, b_rad, d_kpc, degree=False)
    x_kpc = xyz[:, 0]
    y_kpc = xyz[:, 1]
    z_kpc = xyz[:, 2]

    # (x,y,z) --> Galactocentric cylindrical coordinates (R,z,phi):
    Rzphi = bovy_coords.XYZ_to_galcencyl(x_kpc, y_kpc, z_kpc,
                                         Xsun=X_gc_sun_kpc, Zsun=Z_gc_sun_kpc)
    R_kpc = Rzphi[:, 0]
    phi_rad = Rzphi[:, 1]
    z_kpc = Rzphi[:, 2]

    # b. convert velocities (pm_ra,pm_dec,vlos) to (vR,vz,vT)

    # (pm_ra,pm_dec) --> (pm_l,pm_b):
    pmlpmb = bovy_coords.pmrapmdec_to_pmllpmbb(pm_ra_masyr, pm_dec_masyr,
                                               ra_rad, dec_rad, degree=False,
                                               epoch=2000.0)
    pml_masyr = pmlpmb[:, 0]
    pmb_masyr = pmlpmb[:, 1]

    # (v_los,pm_l,pm_b) & (l,b,d) --> (vx,vy,vz):
    vxvyvz = bovy_coords.vrpmllpmbb_to_vxvyvz(v_los_kms, pml_masyr, pmb_masyr,
                                              l_rad, b_rad, d_kpc, XYZ=False,
                                              degree=False)
    vx_kms = vxvyvz[:, 0]
    vy_kms = vxvyvz[:, 1]
    vz_kms = vxvyvz[:, 2]

    # (vx,vy,vz) & (x,y,z) --> (vR,vT,vz):
    vRvTvZ = bovy_coords.vxvyvz_to_galcencyl(vx_kms, vy_kms, vz_kms, R_kpc,
                                             phi_rad, z_kpc,
                                             vsun=[vX_gc_sun_kms,
                                                   vY_gc_sun_kms,
                                                   vZ_gc_sun_kms],
                                             galcen=True)
    vR_kms = vRvTvZ[:, 0]
    vT_kms = vRvTvZ[:, 1]
    vz_kms = vRvTvZ[:, 2]

    print("R = ", R_kpc, "\t kpc")
    print("phi = ", phi_rad, "\t rad")
    print("z = ", z_kpc, "\t kpc")
    print("v_R = ", vR_kms, "\t km/s")
    print("v_T = ", vT_kms, "\t km/s")
    print("v_z = ", vz_kms, "\t km/s")
    return vz_kms

if __name__ == "__main__":

    ra_deg = np.array([7.7750132145])
    dec_deg = np.array([-26.8097293548])
    d_pc = np.array([890.547792917])

    # Velocity of stars:
    pm_ra_masyr = np.array([24.965])  # pm in direction of RA [mas/yr]
    pm_dec_masyr = np.array([-9.683])  # pm in direction of dec [mas/yr]
    v_los_kms = np.array([-4.351])  # line-of-sight velocity [km/s]
    calc_actions(ra_deg, dec_deg, d_pc, pm_ra_masyr, pm_dec_masyr, v_los_kms)
