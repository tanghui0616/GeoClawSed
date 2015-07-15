Document for GeoClawSed variable
dbugfile: debug file
parmfile: parameter file
maxthreads: maximum number of thread
amrfile: AMRClaw parameter file
clawfile: Clawpack paramter file
ndim: Number of dimension 
xlower: 
ylower:
xupper:
yupper:
nx: number of cell in x direction
ny: number of cell in y direction 
nvar or meqn: number of variable
mwaves: number of waves
naux: number of auxiliary array
t0: start run time
output_style: output data style 
nout/nstop: number of output times
tfinal: the final output time
output_t0: the first output time
tout: output time
rinfinity: infinity 
output_format: the way to output
output_q_components: output q components
output_aux_components: output auxilary array components 
output_aux_onlyonce: only output auxilary once
possk: dt initial 
dt_max: largest allowable dt
cflv1: max cfl number
cfl: clf number
nv1: max number of steps
vtime: dt variable or not
method(1):dt variable or not
method(2): order of accuracy 
method(3): transverse waves
method(4): verbosity 
method(5): source term split
dimensional_split: can split dimension or not
mcapa1: capa index
use_fwaves: use fwave or not
mthlim: flux limiter type 
nghost: number of ghost cell
mthbc(1:4): boundary type left:1, right:2, bottom:3, top:4; 1: 0 order; 2: periodic boudary; 3:solid; 4:?; 5: sphere boundary
xperdom/yperdpm: periodic x/y domain
spheredom: sphere domain
rest: restart or not
rstfile: restart file 
checkpt_style: check style: 0: never check, 2: given number of check times; 3: given check time interval
checkpt_interval: check time interval
nchkpt: number of check times
mxnest: the maximum refinement level 
intratx: refinement rate in x direction 
intraty: refinement rate in y direction
kratio: refinement rate for time 
max1d: maximum possible refinement level
auxtype: auxiary array type
flag_richardson: use richardson method or not
tol: error tol
flag_gradient: refine or not
tolsp: spatial error tol
iorder: order of integrator
kcheck: error checking interval
ibuff: buffer zone size
cut:volume ratio cutoff
verbosity_regrid:verbosity for regrid
dprint
eprint
edebug
gprint
nprint
pprint
rprint
sprint
tprint
uprint
mcapa: number of capacity 
hxposs: dx from coarse to fine
hyposs: dy from coarse to fine 
matlabu: frame number for output
tstart_thisrun: restart run time
outfile: output file name
cflmax: max possible cfl number
dxmin: smallest dx for refine level one
dymin: smallest dy for refinement level one
num_gauge_SAVE: store the number of gauge
num_gauges: number of gauge
    
    
    
    
    
