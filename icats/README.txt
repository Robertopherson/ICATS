# just add the icats directory to some PYTHONPATH path directory
# there are example and test script that need the following to work

#to run rundyn.py (outputs directory) and this scripts you need something like:
#pip install pyscf
#pip install pyscf-semiempirical
#pip install -U pyberny

# mindo_pyscf_hess.py* can be used to generate hessians using semiempirical method by providing some molecular xyz...

#to generate the distribuitions one runs:
./test.py input

# the input file reads the molecule input files h2om_dat.txt and ohm_dat.txt, which themselves require the reference geometry (*_geom.xyz)
# and the (if available) mass-scaled hessian (*_hess.txt), there is documentaiton on the input files.. 

#as well as generating the initial conditions files, it generates a "out_dist.log" file which contains, amonst other initial condition information, the distribution 
#of the rotational/vibrational/velocity information. It also generates a logfile with the same information (can be overwritten).

# the outputs/*xyz and *vel contain the initial position and velocities in Angstrom and femtoseconds units. the *info file first prints 
# the Sampled  vib/rot/trans/orient/impactparam  from the MC distributiions and prints their information (as well as the energy breakdown). It then
#subsequently analyses/recalculates all this same information from the *xyz and *vel data; for harmonic approximations and accurate hessians, the analyed
# results should closely match the Sampled results. 

#that generates initial conditions into the "outputs" directory
# the script "runem" then uses rundyn.py on all outputted initial conditions (uses en0 file to rescale the zero of the potential energy)
cd outputs
./runem
cd ..
# then run icats.analyse to analyse concordant *xyz/*vel pairs from dynamics output
# the dynamics*analinfo 
./icats.analyse input --dir outputs --prefix out

#the outputs/*analinfo files runs the analysis vib/rot/trans/orient/space momentum breakdown of each time step available from the dynamics outputted *.md.vel and *.md.xyz data

#you can run this script if the python libraries are present and it should work.. fingers crossed


