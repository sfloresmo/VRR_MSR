# Constants of the code
const G = 1.0 # Gravitational constant
const MBH = 1.0 # Mass of the central BH
const MTOT = 1.0 # Total stellar mass
const NPART = 1 # Total number of particles
const M_STAR = MTOT / NPART # Individual stellar mass
const A_STAR = 1.0 # Semi-major axis of the particles
const ECC_STAR = 0.0 # Eccentricity of the particles
const L_STAR = M_STAR * sqrt(G * MBH * A_STAR * (1 - ECC_STAR^2)) # Norm of the angular momentum of every star
const H2_COUPLING = (pi/5) * (G * M_STAR^(2))/(A_STAR) # \ell=2 Hamiltonian coupling coefficient
const J2_COUPLING = H2_COUPLING / L_STAR # \ell=2 coupling coefficient. ATTENTION, to the rescaling by 1/L
const B2_COUPLING = 30 / (8 * pi) # Constant \ell=2 coefficient B_l = l*(l+1)(2*l+1)/(8*pi)
const C0 = NPART / (4 * pi) # Normalisation constant for the correlation function
const TC = sqrt((4 * pi) / (NPART * B2_COUPLING * J2_COUPLING^2)) # Typical coherence time
