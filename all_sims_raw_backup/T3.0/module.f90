module md_vars
! Bk: Boltzmann Constant,  sigma & e: coefficients in L-J potential (before normalization) of Ar atoms
   real(8), parameter :: Bk = 1.38d-23, sigma = 3.405d-10, e = 119.8d0*Bk
! dmass: mass of an Ar atom in kilogram
   real(8), parameter :: dmass = 4d-25/6.02d0, Pi = 3.1415926d0
! NMAX: max number of atoms, NBMAX: max number of atoms in a reduced cell and its neighboring cells
! NTMAX: max number of points in the force table, LMAX: max number of cells in one direction
! NSMAX: max number of steps, NCMAX: max number of unit cells in one direction
   integer,parameter :: NMAX=2048, NBMAX=4000, NTMAX=10001
   integer,parameter :: NTMAX1=NTMAX - 1, LMAX=20, NGMAX=1000
   integer,parameter :: NSMAX=100000, NCMAX=8
! NMAX3 is the max number of coordinates = 3 * max number of particles (for positions, velocities cand forces) because there are three dimensions.
   integer,parameter :: NMAX3=NMAX*3, NBMAX3=3*NBMAX
! tau=sqrt(dmass/e)*sigma: time scaling constant used for normalization
   real(8),parameter :: tau=2.355d-11, taui=1d0/tau
! dt: time step length in normalized units, dtsq=half of dt square, dti=inverse of dt, dthalf=half dt
   real(8),parameter :: dt=2d-4, dtsq=0.5d0*dt*dt, dti=1d0/dt
   real(8),parameter :: dthalf=0.5d0*dt
! cutoff: rc in the notes,   cutoffsq is its square
   real(8),parameter :: cutoff=2.55d0, cutoffsq=cutoff*cutoff
! dr is the distance between neighboring points in the force table (in normalized unit)
   real(8),parameter :: dr=cutoff/NTMAX1, dri=1d0/dr
! xlc defines the shape of the unit cell, In the case of Ar, it's a cube so all xlc = 1
! v1 v2 v3 are the three unit vectors that define the lattice
! mx my mz form index vectors pointing to the neighboring 13 reduced cells in single direction
! For instance, the first neighboring cell will have index offset of (0, 0, 1) and the 13th will
! have index offset of (-1, 1, 1).
   real(8) :: xlc(3) = (/1.d0, 1.d0, 1.d0/)
   real(8) :: v1(3) = (/5d-1, 5d-1, 0d0/)
   real(8) :: v2(3) = (/0d0, 5d-1, 5d-1/)
   real(8) :: v3(3) = (/5d-1, 0d0, 5d-1/)
   integer :: mx(14) = (/0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, -1/)
   integer :: my(14) = (/0, 0, 1, 0, 1, 0, 1, 1, 0, -1, 1, 1, -1, 1/)
   integer :: mz(14) = (/0, 1, 0, 0, 1, 1, 0, -1, -1, 0, 1, -1, 1, 1/)

! istep is the number of time step, istable is from which step to start calculating Cv and
! diffusion.
   integer :: istep, istable, iprint
! req_temp: requested temperature in normalized units. coeff_quench is for quenching temperature
   real(8) :: req_temp, coeff_quench
! xlcd: Ar lattice constant in scaled coordinates. ncell are the numbers of unit cells in x, y and z direction
! ini_new: if we wanna generate new system (1 generates new configurations and 0 means read in).
   real(8) :: xlcd
   integer :: ncell(3), ini_new, iscale_x, iscale_v, isquench, iquench
! x are the coordinates of atoms saved in a one dimensional array. v are the velocities and f the forces
   real(8) :: x(NMAX3), f(NMAX3), v(NMAX3), originx(NMAX3)
! epot is the potential energy and ekin is the kinetic energy, with akin and akinsq for calculating Cv
! density, volume, pressure are obvious
   real(8) :: epot, ekin, etot, diffusion(5000)
   real(8) :: akin, akinsq, volume, pressure, density, energy(2000, 3)
! boxmd are the size of the total md box, cellsize are the size of a reduced cell, mc are the numbers of
! reduced cells. They are all vectors in 3 dimensions. halfboxmd = 0.5 * boxmd
   real(8) :: boxmd(3), cellsize(3), halfboxmd(3)
   integer :: mc(3)
! nheader is the header of linked list. lsize stores the number of atoms in each reduced cell
! linklist is the body of the linked list
   integer :: nheader(LMAX, LMAX, LMAX), lsize(LMAX, LMAX, LMAX), linklist(NMAX)
! ftable is the force table and vtable is the potential table. They are constructed and used through
! interpolation to save number of calculations and computing time.
   real(8) :: ftable(NTMAX), vtable(NTMAX)
! localN is the number of atoms in the md box, and volume, pressure are the system's properties
   integer :: localN
! This is the seed used to generate random numbers
   real(8) :: dseed
! These are the offsets applied to ftable and vtable to make them continuous and smooth at the cutoff
   real(8) :: voffset, foffset
! g(r) and coordination number array
   real(8) :: g(NGMAX), cnum(NGMAX), gdr, grcutoff, grcutoffsq

! # of unit cells used in the g(r) cutoff
   integer,parameter :: NUMXLCDFORGR=3

! freeze atoms from bottom in z
   real(8) :: lfreeze=0.0d0
   real(8),allocatable :: msd_by_layer(:,:)
   real(8) :: msd_layer_thickness=2.d0
   integer :: msd_num_layers

end module
