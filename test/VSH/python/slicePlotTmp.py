import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from pytools.vtk_io import readVTK
import inifix
import os

rc('text.latex', preamble = r'\usepackage{txfonts}')
rc('text', usetex=True, color='black') ; rc('font', family='serif', serif='Times')

conf = inifix.load('../idefix.ini')
gamma = conf["Hydro"]["gamma"]
rmin = conf["Grid"]["X1-grid"][1]
rmax = conf["Grid"]["X1-grid"][4]
nphi = conf["Grid"]["X3-grid"][2]

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-fcb", "--fixed_colorbar", default=-1, help="use a fixed colorbar, set on the rms from the given file")
parser.add_argument("-all_slice", "--make_all_slice", default=False, help="Make all slices", action='store_true')
parser.add_argument("-r", "--radial_slice", default=False, help="Make radial slice plots for all vtk files", action='store_true')
parser.add_argument("-th", "--theta_slice", default=False, help="Make meridional slice plots for all vtk files", action='store_true')
parser.add_argument("-phi", "--phi_slice", default=False, help="Make equatorial slice plots for all vtk files", action='store_true')
parser.add_argument("-mhd", "--make_mhd_plot", default=False, help="Make magnetic field plot as well", action='store_true')
parser.add_argument("-mhd_only", "--make_only_mhd_plot", default=False, help="Make magnetic field plot only", action='store_true')
parser.add_argument("-all_vtk", "--make_all_vtk", default=False, help="Make plot for all vtk", action='store_true')
parser.add_argument("-all_var", "--make_all_var", default=False, help="Make plot for all variables", action='store_true')
parser.add_argument("-all_proj", "--make_all_proj", default=False, help="Make plot with all projections", action='store_true')
parser.add_argument("-vtk", "--vtk_index", nargs='+', default=-1, help="Vtk indexes to be plotted")
parser.add_argument("-var", "--var_name", nargs='+', default='', help="Physical variables to be plotted")
parser.add_argument("-proj", "--proj_name", nargs='+', default='', help="Projections to be used in radial slices")
parser.add_argument("-all", "--make_all", default=False, help="Make all previous parser to true", action='store_true')
parser.add_argument("-show", "--show_plot", default=False, help="Show the plots", action='store_true')
parser.add_argument("-save", "--save_plot", default=False, help="Save the plots", action='store_true')
parser.add_argument("-gif", "--make_gif", default=False, help="Make a gif from plots", action='store_true')
args = parser.parse_args()

variables = {"VX1": 0, "VX2": 1, "VX3": 2, "DTT": 3, "DPP": 4, "DRR": 4, "VEL": 5, "MACH": 6, "BX1": 7, "BX2": 8, "BX3": 9, "B": 10}
baseclb = {"VX1": 1., "VX2": 1., "VX3": 1., "DTT": 1., "DPP": 1., "DRR": 1., "VEL":0., "MACH": 0., "BX1": 1., "BX2": 1., "BX3": 1., "B": 0.}
cmapclb = {"VX1": 'PuOr', "VX2": 'PRGn', "VX3": 'PiYG', "DTT": 'coolwarm', "DPP": 'Spectral', "DRR": 'RdBu', "VEL": 'PuRd', "MACH": 'Reds', "BX1": 'BrBG', "BX2": 'RdGy', "BX3": 'RdYlBu', "B": 'magma'}
what_title = {"VX1": 'Radial velocity', "VX2": 'Polar velocity', "VX3": 'Azimuthal velocity', "DTT": 'Temperature fluctuation', "DPP": 'Pressure fluctuation', "DRR": 'Density fluctuation', "VEL": 'Velocity norm', "MACH": 'Mach number', "BX1": 'Radial magnetic field', "BX2": 'Polar magnetic field', "BX3": 'Azimuthal magnetic field', "B": 'Magnetic field strength'}
symbol = {"VX1": r'\varv_r', "VX2": r'\varv_\theta', "VX3": r'\varv_\varphi', "DTT": r'\delta T', "DPP": r'\delta P', "DRR": r'\delta\rho', "VEL": r'\varv', "MACH": r'\mathcal{M}', "BX1": r'B_r', "BX2": r'B_\theta', "BX3": r'B_\varphi', "B": r'B'}
#variables = {"VX1": 0, "VX2": 1, "VX3": 2, "deltaT_T": 3, "deltaP_P": 4, "deltaRho_Rho": 4, "VEL": 5, "MACH": 6, "BX1": 7, "BX2": 8, "BX3": 9, "B": 10}
#baseclb = {"VX1": 1., "VX2": 1., "VX3": 1., "deltaT_T": 1., "deltaP_P": 1., "deltaRho_Rho": 1., "VEL":0., "MACH": 0., "BX1": 1., "BX2": 1., "BX3": 1., "B": 0.}
#cmapclb = {"VX1": 'PuOr', "VX2": 'PRGn', "VX3": 'PiYG', "deltaT_T": 'coolwarm', "deltaP_P": 'Spectral', "deltaRho_Rho": 'RdBu', "VEL": 'PuRd', "MACH": 'Reds', "BX1": 'undefined', "BX2": 'undefined', "BX3": 'undefined', "B": 'magma'}
#what_title = {"VX1": 'Radial velocity', "VX2": 'Polar velocity', "VX3": 'Azimuthal velocity', "deltaT_T": 'Temperature fluctuation', "deltaP_P": 'Pressure fluctuation', "deltaRho_Rho": 'Density fluctuation', "VEL": 'Velocity norm', "MACH": 'Mach number', "BX1": 'Radial magnetic field', "BX2": 'Polar magnetic field', "BX3": 'Azimuthal magnetic field', "B": 'Magnetic field strength'}
#symbol = {"VX1": r'\varv_r', "VX2": r'\varv_\theta', "VX3": r'\varv_\varphi', "deltaT_T": r'\delta T', "deltaP_P": r'\delta P', "deltaRho_Rho": r'\delta\rho', "VEL": r'\varv', "MACH": r'\mathcal{M}', "BX1": r'B_r', "BX2": r'B_\theta', "BX3": r'B_\varphi', "B": r'B'}
nvar = len(variables)

for varname in list(variables.keys()):
        exec('{0} = {1}'.format(varname, variables[varname]))

fcb = np.int32(args.fixed_colorbar)
if args.make_all:
    make_all_slice = True
    make_all_vtk = True
    make_all_var = True
    make_all_proj = True
else:
    make_all_slice = args.make_all_slice
    make_all_vtk = args.make_all_vtk
    make_all_var = args.make_all_var
    make_all_proj = args.make_all_proj
if make_all_slice:
    radial_slice = True
    theta_slice = True
    phi_slice = True
else:
    radial_slice = args.radial_slice
    theta_slice = args.theta_slice
    phi_slice = args.phi_slice
if radial_slice:
    import cartopy.crs as ccrs
    if make_all_proj:
        projections = [ccrs.Mollweide(), ccrs.Orthographic(), ccrs.Orthographic(central_latitude=90), ccrs.Orthographic(central_latitude=-90)]
        projections_name = ['MOLL', 'ORTHO', 'ORTHONP', 'ORTHOSP']
    else: 
        projections = list()
        projections_name = list()
        if 'moll' in args.proj_name:
            projections.append(ccrs.Mollweide())
            projections_name.append('MOLL')
        if 'ortho' in args.proj_name:
            projections.append(ccrs.Orthographic())
            projections_name.append('ORTHO')
        if 'orthonp' in args.proj_name:
            projections.append(ccrs.Orthographic(central_latitude=90))
            projections_name.append('ORTHONP')
        if 'orthosp' in args.proj_name:
            projections.append(ccrs.Orthographic(central_latitude=-90))
            projections_name.append('ORTHOSP')
        if len(projections) == 0:
            raise Exception("Need to specify at least one projection for radial slices")    
else:
    radial_slice = args.radial_slice
    theta_slice = args.theta_slice
    phi_slice = args.phi_slice
if make_all_vtk:
    tstop = conf["TimeIntegrator"]["tstop"]
    time_step = conf["Output"]["vtk"]
    nfile = np.int32(tstop//time_step + 1)
    list_time = [i*time_step for i in range(nfile)]
    vtk_beg_idx = 0
    vtk_end_idx = nfile
    step = 1
    vtk_index = np.arange(vtk_beg_idx, vtk_end_idx, step)
elif args.vtk_index==-1:
    raise Exception("VTK indexes should be specified when not all vtk need to be plotted")
else:
    try:
        try:
            vtk_beg_idx, vtk_end_idx, step = np.int32(args.vtk_index[0]), np.int32(args.vtk_index[1]), np.int32(args.vtk_index[2])
        except:
            vtk_beg_idx, vtk_end_idx = np.int32(args.vtk_index[0]), np.int32(args.vtk_index[1])
            step = 1
    except:
        vtk_beg_idx = np.int32(args.vtk_index[0])
        vtk_end_idx = vtk_beg_idx + 1
        step = 1
    vtk_index = np.arange(vtk_beg_idx, vtk_end_idx, step)
    nfile = (vtk_end_idx - vtk_beg_idx)//step
    time_step = conf["Output"]["vtk"]*step
    list_time = [(i+vtk_beg_idx)*time_step for i in range(nfile)]
if make_all_var and args.make_mhd_plot:
    varvalues = list(variables.values())
    varnames = list(variables.keys())
elif args.make_only_mhd_plot:
    varvalues = list(variables.values())[-4:]
    varnames = list(variables.keys())[-4:]
elif make_all_var:
    varvalues = list(variables.values())[:-4]
    varnames = list(variables.keys())[:-4]
else:
    new_varnames = list()
    new_varvalues = list()
    for varname in args.var_name:
        if varname in args.var_name:
            new_varnames.append(varname)
            new_varvalues.append(variables[varname])
    varnames = new_varnames
    varvalues = new_varvalues
print("Plot information")
print("List of variables: {0}".format(varnames))
print("List of VTK files: {0}".format(vtk_index))
print("Corresponding simulation time: {0}".format(list_time))
if args.make_gif:
    print("Creating animated gifs as well")

if fcb == -2:
    fcb = np.int32(vtk_end_idx-1)
    prefactor = 1.3
    folder = "sliceplot_{0}fcb".format(fcb)
if fcb == -1:
    prefactor = 2.
    folder = "sliceplot_nofcb"
if not os.path.exists(folder):
    os.mkdir(folder)

fgs = 3.5

VTK = readVTK('../data.{0:04d}.vtk'.format(0), 'spherical')
R = VTK.r
TH = VTK.theta
PHI = VTK.phi
RR, THTH, PHIPHI = np.meshgrid(R,TH,PHI, indexing='ij')
RR = np.moveaxis(RR, (2,0), (0,2))
THTH = np.moveaxis(THTH, (2,0), (0,2))
PHIPHI = np.moveaxis(PHIPHI, (2,0), (0,2))
XX = RR*np.sin(THTH)*np.cos(PHIPHI)
YY = RR*np.sin(THTH)*np.sin(PHIPHI)
ZZ = RR*np.cos(THTH)
nr = R.size
ntheta = TH.size
nphi = PHI.size
PHY = np.zeros((nvar, nfile, nphi, ntheta, nr))
for vtk_idx in vtk_index:
    VTK = readVTK('../data.{0:04d}.vtk'.format(vtk_idx), 'spherical')
    for var_idx, var in enumerate(varvalues):
        varname = varnames[var_idx]
        if varname not in ['VEL', 'MACH', 'B']:
            PHY[var,vtk_idx-vtk_beg_idx] = np.float64(np.moveaxis(VTK.data[varname], (0,2), (2,0)))
        elif varname == 'VEL':
            PHY[var,vtk_idx-vtk_beg_idx] = np.float64(np.moveaxis(VTK.data['VX1'], (0,2), (2,0))**2)
            PHY[var,vtk_idx-vtk_beg_idx] += np.float64(np.moveaxis(VTK.data['VX2'], (0,2), (2,0))**2)
            PHY[var,vtk_idx-vtk_beg_idx] += np.float64(np.moveaxis(VTK.data['VX3'], (0,2), (2,0))**2)
            PHY[var,vtk_idx-vtk_beg_idx] = np.float64(np.sqrt(PHY[var,vtk_idx-vtk_beg_idx]))
        elif varname == 'B':
            PHY[var,vtk_idx-vtk_beg_idx] = np.float64(np.moveaxis(VTK.data['BX1'], (0,2), (2,0))**2)
            PHY[var,vtk_idx-vtk_beg_idx] += np.float64(np.moveaxis(VTK.data['BX2'], (0,2), (2,0))**2)
            PHY[var,vtk_idx-vtk_beg_idx] += np.float64(np.moveaxis(VTK.data['BX3'], (0,2), (2,0))**2)
            PHY[var,vtk_idx-vtk_beg_idx] = np.float64(np.sqrt(PHY[var,vtk_idx-vtk_beg_idx]))
        elif varname == 'MACH':
            PHY[var,vtk_idx-vtk_beg_idx] = PHY[VEL,vtk_idx-vtk_beg_idx]/np.float64(np.sqrt(np.moveaxis(VTK.data['PRS']/VTK.data['RHO']**gamma, (0,2), (2,0))))
    
def make_slice(vtk_idx):
    current_time = list_time[vtk_idx-vtk_beg_idx]
    if theta_slice:
        ith_ref = ntheta//2
        where = 'in the equatorial plane'
        current_XX = XX[:,ith_ref,:]
        current_YY = YY[:,ith_ref,:]
        xy_anno = (0.9*current_XX.min(), 0.9*current_YY.max())
        for idx, var in enumerate(varvalues):
            if fcb == -1:
                rms_value = np.sqrt(np.mean(PHY[var,vtk_idx-vtk_beg_idx,:,ith_ref,:]**2.))*prefactor
            else:
                rms_value = np.sqrt(np.mean(PHY[var,fcb-vtk_beg_idx,:,ith_ref,:]**2.))*prefactor
            varname = varnames[idx]
            what = what_title[varname]
            fig = plt.figure()
            ax = fig.add_subplot()
            p = ax.pcolor(current_XX, current_YY, PHY[var,vtk_idx-vtk_beg_idx,:,ith_ref,:], vmin=-rms_value*baseclb[varname], vmax=rms_value, cmap=cmapclb[varname])
            ax.set_title('{0} {1}'.format(what, where))
            ax.annotate(r't = {0:.1f}'.format(current_time), xy_anno)
            cbar = fig.colorbar(p)
            cbar.set_label(r'${0}$'.format(symbol[varname]))
            if args.show_plot:
                plt.show()
            if args.save_plot:
                fig.savefig('{0}/THETA_{1}_{2:03d}.png'.format(folder,varname,vtk_idx))
            plt.close(fig)

    if phi_slice:
#        iphi_ref = 0
#        current_XX = XX[iphi_ref,:,:]
#        current_ZZ = ZZ[iphi_ref,:,:]
        iphi_ref = nphi//2 ### to get a yt-like plot
        current_XX = XX[0,:,:]
        current_ZZ = ZZ[0,:,:]
        where = 'in the polar plane'
        xy_anno = (0.8*current_XX.max(), 0.9*current_ZZ.max())
        for idx, var in enumerate(varvalues):
            if fcb == -1:
                rms_value = np.sqrt(np.mean(PHY[var,vtk_idx-vtk_beg_idx,iphi_ref,:,:]**2.))*prefactor
            else:
                rms_value = np.sqrt(np.mean(PHY[var,fcb-vtk_beg_idx,iphi_ref,:,:]**2.))*prefactor
            varname = varnames[idx]
            what = what_title[varname]
            fig = plt.figure(figsize=(1.3*fgs, 2.*fgs))
            ax = fig.add_subplot()
            p = ax.pcolor(current_XX, current_ZZ, PHY[var,vtk_idx-vtk_beg_idx,iphi_ref,:,:], vmin=-rms_value*baseclb[varname], vmax=rms_value, cmap=cmapclb[varname])
            ax.set_title('{0} {1}'.format(what, where))
            ax.annotate(r't = {0:.1f}'.format(current_time), xy_anno)
            ax.set_xlim(0., rmax)
            ax.set_ylim(-rmax, rmax)
            cbar = fig.colorbar(p)
            cbar.set_label(r'${0}$'.format(symbol[varname]))
            if args.show_plot:
                plt.show()
            if args.save_plot:
                fig.savefig('{0}/PHI_{1}_{2:03d}.png'.format(folder,varname,vtk_idx))
            plt.close(fig)

    if radial_slice:
        ir_ref = nr//2
        where = 'in the middle shell'
        current_TH = np.linspace(np.pi/2., -np.pi/2., ntheta)
        current_PHI = np.linspace(0., 2.*np.pi, nphi)
        current_PHIPHI, current_THTH = np.meshgrid(PHI, current_TH)
        current_THTH = np.rad2deg(current_THTH)
        current_PHIPHI = np.rad2deg(current_PHIPHI)
        for idx_proj, projection in enumerate(projections):
            proj_name = projections_name[idx_proj]
            xy_anno = (-10., 0.)
            for idx, var in enumerate(varvalues):
                if fcb == -1:
                    rms_value = np.sqrt(np.mean(PHY[var,vtk_idx-vtk_beg_idx,:,:,ir_ref]**2.))*prefactor
                else:
                    rms_value = np.sqrt(np.mean(PHY[var,fcb-vtk_beg_idx,:,:,ir_ref]**2.))*prefactor
                varname = varnames[idx]
                what = what_title[varname]
                fig = plt.figure()
                ax = fig.add_subplot(projection=projection)
    #            p = ax.pcolor(current_PHIPHI, current_THTH, PHY[var,vtk_idx-vtk_beg_idx,:,:,ir_ref].T, vmin=-rms_value*baseclb[varname], vmax=rms_value, cmap=cmapclb[varname], transform=ccrs.PlateCarree())
                p = ax.pcolor(current_PHIPHI, current_THTH, np.roll(PHY[var,vtk_idx-vtk_beg_idx,:,:,ir_ref].T, nphi//2, axis=1), vmin=-rms_value*baseclb[varname], vmax=rms_value, cmap=cmapclb[varname], transform=ccrs.PlateCarree()) ### to get a yt-like plot
                ax.set_title('{0} {1}'.format(what, where))
    #            ax.annotate(r't = {0:.1f}'.format(current_time), xy_anno)
                cbar = fig.colorbar(p)
                cbar.set_label(r'${0}$'.format(symbol[varname]))
                if args.show_plot:
                    plt.show()
                if args.save_plot:
                    fig.savefig('{0}/RAD_{1}_{2}_{3:03d}.png'.format(folder,proj_name, varname,vtk_idx))
                plt.close(fig)

if args.save_plot and not args.show_plot:
    for vtk_idx in vtk_index:
        make_slice(vtk_idx)

rs = [5,100]
normals = list()
if radial_slice:
    normals.append('RAD')
if theta_slice:
    normals.append('THETA')
if phi_slice:
    normals.append('PHI')

if args.make_gif:
    for r in rs:
        for varname in varnames:
            for normal in normals:
                if normal == 'RAD':
                    for proj_name in projections_name:
                        os.system("ffmpeg -r {1} -i {0}/{2}_{4}_{3}_%3d.png -loop 0 {0}/{2}_{4}_{3}_r{1}.gif -y".format(folder, r, normal, varname, proj_name))
                else:
                    os.system("ffmpeg -r {1} -i {0}/{2}_{3}_%3d.png -loop 0 {0}/{2}_{3}_r{1}.gif -y".format(folder, r, normal, varname))

