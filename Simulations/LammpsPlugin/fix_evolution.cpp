#include "fix_evolution.h"
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"

#define PI 3.14159265


using namespace LAMMPS_NS;
using namespace FixConst;


Evolution2D::Evolution2D(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
  
{
  
/*
#########################################################################

    Parameters definition
    
#########################################################################
*/
 
  if (narg != 15) error->all(FLERR,"Illegal evolution2D command. \n Expected arguments : eta - D - tau_v - tau_n - J - epsilon - alpha - alphaq - comm_radius - region - controller - seed");
  eta = utils::numeric(FLERR,arg[3],false,lmp);
  D = utils::numeric(FLERR,arg[4],false,lmp);
  if (eta < 0.0) error->all(FLERR,"Diffusion coefficient eta must be >= 0.0");
  if (D < 0.0) error->all(FLERR,"Diffusion coefficient D must be >= 0.0");
  
  tau_v = utils::numeric(FLERR,arg[5],false,lmp);
  tau_n = utils::numeric(FLERR,arg[6],false,lmp);
  inertiaJ = utils::numeric(FLERR,arg[7],false,lmp);
  epsilon = utils::numeric(FLERR,arg[8],false,lmp);
  alpha = utils::numeric(FLERR,arg[9],false,lmp);
  alphaq = utils::numeric(FLERR,arg[10],false,lmp);
  comm_radius = utils::numeric(FLERR,arg[11],false,lmp);
  region = domain->get_region_by_id(arg[12]);
  idregion = utils::strdup(arg[12]);
  controller_flag = utils::numeric(FLERR,arg[13],false,lmp);
  seed = utils::numeric(FLERR,arg[14],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal evolution2D command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
}

Evolution2D::~Evolution2D()
{
  
  delete random;

}



int Evolution2D::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}


void Evolution2D::init()
{

}


void Evolution2D::initial_integrate(int vflag)
{
 
  // function to update x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **mu = atom->mu;
  double **omega = atom->omega;
  double **phi = atom->phi;
  double *q_reward = atom->q_reward;
  double mag;

  double F_min = 0;
  double D_min = 0;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int step = update->ntimestep;
  
  double q_received = 0.0;
  double comm_sq = comm_radius * comm_radius;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dt = update->dt;
  sqrtdt = sqrt(dt);
  
  // Set initial particle orientation
  if (step <= 1) {
    for (int i = 0; i < nlocal; i++){       
    
/*
#########################################################################

    Initial conditions
    
#########################################################################
*/
       // Initialize Active vector with random direction at beginining of simulation
       mu[i][0] = random->gaussian();       
       mu[i][1] = random->gaussian(); 
       mu[i][2] = 0; // 2D

       // Normalize active vector
       mag = sqrt(pow(mu[i][0],2)+pow(mu[i][1],2));
       double mag_inv = 1.0 / mag;
       mu[i][0] *= mag_inv;
       mu[i][1] *= mag_inv;
    
       // Initialize parameters with random direction at beginining of simulation
       phi[i][0] = random->uniform();
       phi[i][1] = random->uniform();
       phi[i][2] = random->uniform();
       q_reward[i] = 0.0;
       }
  }
  int *ilist,*jlist,*numneigh,**firstneigh;
  int i, j, ii, jj, inum, jnum;
  double rsq, xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  NeighList *list = neighbor->lists[0];  // Works but might break
 
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  // Update parameters and reward
  if(step > 1){
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely;
      if(rsq == 0) continue;

      if (rsq < comm_sq) {

/*
#########################################################################

    Communication between j and i
    
#########################################################################
*/

        if (q_reward[j]  >= q_reward[i]){
	        phi[i][0] += alpha* (phi[j][0] - phi[i][0]) * dt;
	        phi[i][1] += alpha* (phi[j][1] - phi[i][1]) * dt;
	        phi[i][2] += alpha* (phi[j][2] - phi[i][2]) * dt;
	        q_reward[i] += alpha* (q_reward[j] - q_reward[i]) * dt;
	      }
	      if (q_reward[i]  >= q_reward[j]){
	        phi[j][0] += alpha* (phi[i][0] - phi[j][0]) * dt;
	        phi[j][1] += alpha* (phi[i][1] - phi[j][1]) * dt;
	        phi[j][2] += alpha* (phi[i][2] - phi[j][2]) * dt;
		      q_reward[j] += alpha* (q_reward[i] - q_reward[j]) * dt;
	      }
      }
    }

  }
  }
  // Integrator 
  for (int i = 0; i < nlocal; i++)
	  
    if (mask[i] & groupbit) {


/*
#########################################################################

    Model definition
    
#########################################################################
*/
      // Mutation noise
      phi[i][0] += random->gaussian() * sqrt(2*dt*eta); 
      phi[i][1] += random->gaussian() * sqrt(2*dt*eta); 
      phi[i][2] += random->gaussian() * sqrt(2*dt*eta);
      
      // Clamp values with hard collisions to 0 and 1
      if(phi[i][0] > 1.0) phi[i][0] = 2 - phi[i][0];
      if(phi[i][0] < 0.0) phi[i][0] = - phi[i][0];
      
      if(phi[i][1] > 1.0) phi[i][1] = 2 - phi[i][1];
      if(phi[i][1] < 0.0) phi[i][1] = - phi[i][1];
      
      if(phi[i][2] > 1.0) phi[i][2] = 2 - phi[i][2];
      if(phi[i][2] < 0.0) phi[i][2] = - phi[i][2];
        
      if(region->match(x[i][0], x[i][1], x[i][2])){
        q_received = 0.75;
      } else {
        q_received = 0.25;
      }
      // Update reward
      q_reward[i] += alphaq*(q_received - q_reward[i]) *dt;

      
      double Fa = 0;
      double th = tanh((phi[i][0]- q_received)/phi[i][1]);
      switch(controller_flag){
      
        // ###### AGNOSTIC CONTROLLERS ###### No response to light
        case 0:
        Fa = phi[i][0];
        break;
        
        case 1:
        Fa = 0.5*(1+cos(PI * (phi[i][0]+1)));
        break;
        
        case 2:
        Fa = 0.5*(1+cos(2*PI*phi[i][0]));
        break;
        
        case 3:
        Fa = -2*phi[i][0] + 1;
        if(phi[i][0] >= 0.5) Fa = 2*phi[i][0] - 1;
        break;
        
        case 4:
        Fa = 0.5*(1+cos(2*PI*phi[i][0] + PI));
        break;
        
        case 5:
        Fa = 2*phi[i][0];
        if(phi[i][0] >= 0.5) Fa = -2*phi[i][0] + 2;
        break;
        
        // ###### PRESETS ######
        case 6:
        // Preset still
        if(q_received >= 0.5) Fa = 0.;
        else Fa = 1.;
        break;
        
        case 7:
        // Preset move
        if(q_received >= 0.5) Fa = 0.33;
        else Fa = 1.;
        break;
        
        // ###### LEARNING ######
        case 8:
        // Bi-modal
        if(q_received >= 0.5) Fa = phi[i][1];
        else Fa = phi[i][0];
        break;
        
        case 9:
        // Threshold still
        if(q_received >= phi[i][0]) Fa = 0.0;
        else Fa = 1.;
        break;
        
        case 10:
        // Threshold move
        if(q_received >= phi[i][0]) Fa = 0.33;
        else Fa = 1.;
        break;
        
        case 11:
        // TanH
        Fa = 0.5 * (1+th) * phi[i][2] + 0.5 * (1-th) * (1-phi[i][2]);
        break;
        
        default:
        Fa = 1.;
        break;
      }
    

      // Update velocities
      v[i][0] += dt*(Fa*mu[i][0] - v[i][0] + f[i][0])/tau_v;
      v[i][1] += dt*(Fa*mu[i][1] - v[i][1] + f[i][1])/tau_v;
      v[i][2] = 0;

      // overdamped
      /*v[i][0] = mu[i][0] + f[i][0];
      v[i][1] = mu[i][1] + f[i][1];
      v[i][2] = 0;*/

       // Update Active Vector
      double mux, muy;
      if(inertiaJ == 0){
        mux = epsilon * (mu[i][1]*mu[i][1]*v[i][0] - mu[i][0]*mu[i][1]*v[i][1]) / tau_n;
        muy = epsilon * (mu[i][0]*mu[i][0]*v[i][1] - mu[i][0]*mu[i][1]*v[i][0]) / tau_n;
        mu[i][0] += dt * mux;
        mu[i][1] += dt * muy;
      } else {
        omega[i][2] += dt * (epsilon * (mu[i][0]* v[i][1] - mu[i][1]*v[i][0]) - tau_n*omega[i][2]) / inertiaJ;
        //omega[i][2] += dt /inertiaJ;
        mu[i][0] += -dt * mu[i][1] * omega[i][2];
        mu[i][1] += dt * mu[i][0] * omega[i][2] ;
      }
      
      // Initialise angular noise
      double ang_noise = random->gaussian();
      ang_noise *= sqrt(2*D*dt);
      
      mux = mu[i][0];
      muy = mu[i][1];

      mu[i][0] = cos(ang_noise) * mux - sin(ang_noise) * muy;
      mu[i][1] = sin(ang_noise) * mux + cos(ang_noise) * muy;

      /*double vx, vy;
      vx = v[i][0];
      vy = v[i][1];

      v[i][0] = cos(ang_noise) * vx - sin(ang_noise) * vy;
      v[i][1] = sin(ang_noise) * vx + cos(ang_noise) * vy;*/
      
      
      // Normalise updated Active vector
      mag = sqrt(pow(mu[i][0],2)+pow(mu[i][1],2)+pow(mu[i][2],2));
      double mu_mag_inv = 1.0/mag;
      mu[i][0] *= mu_mag_inv;
      mu[i][1] *= mu_mag_inv;
      
      // Update positions
      x[i][0] +=  v[i][0]*dt;
      x[i][1] +=  v[i][1]*dt;
      x[i][2] = 0;    
    }
}

/* ---------------------------------------------------------------------- */
