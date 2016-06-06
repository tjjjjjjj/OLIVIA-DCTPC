#include "SimRings.hh"
#include <cmath>
#include <iostream>


ClassImp(SimRings); 

SimRings::SimRings(int n, double x0, double y0, double rmin, double rmax, double zmin, double zmax, double dz, double ztop, double surftop, double zbottom, double surfbottom)
{
  nrings = n; 
  this->x0 = x0;
  this->y0 = y0;
 
  this->dz = dz; 
  this->rin = rmin; 
  this->rout = rmax; 

  this->ztop = ztop; 
  this->zbottom = zbottom ; 
  this->surftop = surftop; 
  this->surfbottom = surfbottom; 

  double currentz = 0; 
  zgap = ((zmax-zmin) - n * dz) / (n-1.); 
  double ringarea = M_PI * (rmax *rmax - rmin*rmin) * n; 
  double tubearea = 2*M_PI*rmin * dz * n; 

  double top_area = M_PI*rmin*rmin * surftop; 
  double bottom_area = M_PI*rmin*rmin * surfbottom; 

  dr = rout - rin; 

  ring_fraction = ringarea / (tubearea + ringarea + top_area + bottom_area); 
  top_fraction = top_area / (tubearea + ringarea + top_area + bottom_area); 
  bottom_fraction = bottom_area / (tubearea + ringarea + top_area + bottom_area); 


  std::cout << "Ring Fraction: " << ring_fraction << std::endl; 
  std::cout << "Top Fraction: " << top_fraction << std::endl; 
  std::cout << "Bottom Fraction: " << bottom_fraction << std::endl; 

  for (int i = 0; i < n; i++)
  {
    zs.push_back(currentz);                
    zs.push_back(currentz+dz);                
    currentz+= zgap; 
  }

}

void SimRings::generateDecay(double * pos, double * dir)
{
  //Generate until we reach something valid. 
  while(true)
  {
    //Pick horizontal or vertical     

    double surface = rand.Rndm(); 

    if (surface < ring_fraction)
    {
      //on ring. pick which surface: 
      int i = (int)  (rand.Rndm() * 2 * nrings); 
      bool top = i % 2; 

      double z = zs[i]; 
      double r  = sqrt(rand.Rndm()*(rout*rout - rin*rin) + rin*rin);  
      
      std::cout<<"generated r "<<r<<std::endl;

      double p_theta = acos(rand.Rndm() * (top ? 1 : -1)); 
      double phi = rand.Rndm() * 2 * M_PI; 
      double p_phi = rand.Rndm() * M_PI + M_PI/2 + phi; //only consider directions towards inside. 
      double px = cos(p_phi) * sin(p_theta);
      double py = sin(p_phi) * sin(p_theta);
      double pz = cos(p_theta); 

      if (pz == 0) continue; //No pure horizontal possible

      double x = r * cos(phi); 
      double y = r * sin(phi); 

      // Solve for when we hit the inner radius
     
      double A = px*px + py*py; 
      double B = 2*px*x + 2*py*y; 
      double C = x*x + y*y - rin*rin; 
     
      double discr = B*B - 4*A*C; 
      if (discr <= 0) continue;  //Doesn't go inside rings or just barely grazes

      double tplus = (-B + sqrt(discr)) / (2*A); 
      double tminus = (-B - sqrt(discr)) / (2*A); 

      if (tplus < 0 || tminus < 0)
      {
        std::cerr << "ERROR: PARTICLE GOING WRONG WAY!!!" << std::endl; 
        continue; 
      }

      if (top != pz > 0)  
      {
        std::cerr << "ERROR: PARTICLE INTO OWN RING!!!" << std::endl; 
        continue; 
      }

      //Pick minimum time
      double t = tplus < tminus ? tplus : tminus; 
  
      
      //std::cout << pz*t<< " " << top << " " << zgap <<std::endl; 

      //make sure we don't hit the other ring
      if (fabs(pz*t) > zgap) continue; 

      //Yay, it worked! 
      pos[0] = x + x0;  
      pos[1] = y + y0;  
      pos[2] = z;  
      dir[0] = px; 
      dir[1] = py; 
      dir[2] = pz; 
      //std::cout << "outside: " << p_theta << std::endl; 
      break; 

    }
    else if ( surface > ring_fraction && surface < ring_fraction + top_fraction)
    {
      double r = sqrt(rand.Rndm()) * rin; 
      double phi = rand.Rndm() * 2*M_PI; 
      double z = ztop; 
      double p_theta = acos(rand.Rndm()); 
      double p_phi = rand.Rndm() * 2 * M_PI; 

      double x= r * cos(phi); 
      double y= r * sin(phi); 
      double px = cos(p_phi) * sin(p_theta);
      double py = sin(p_phi) * sin(p_theta); 
      double pz = cos(p_theta); 

      pos[0] = x + x0;  
      pos[1] = y + y0;  
      pos[2] = z;  
      dir[0] = px; 
      dir[1] = py; 
      dir[2] = pz; 
      break;
    }
    else if ( surface > ring_fraction + top_fraction && surface < ring_fraction + top_fraction + bottom_fraction)
    {
      double r = sqrt(rand.Rndm()) * rin; 
      double phi = rand.Rndm() * 2*M_PI; 
      double z = zbottom; 
      double p_theta = acos(-rand.Rndm()); 
      double p_phi = rand.Rndm() * 2 * M_PI; 
      double x= r * cos(phi); 
      double y= r * sin(phi); 
      double px = cos(p_phi) * sin(p_theta);
      double py = sin(p_phi) * sin(p_theta); 
      double pz = cos(p_theta); 
      pos[0] = x + x0;  
      pos[1] = y + y0;  
      pos[2] = z;  
      dir[0] = px; 
      dir[1] = py; 
      dir[2] = pz; 
      break;
    }
    else
    {
      //On inside
      int i = (int)  (rand.Rndm() * nrings); 
      double phi = rand.Rndm() * 2*M_PI; 
      double z = rand.Rndm()*dz + zs[2*i]; 
      double p_theta = acos(rand.Rndm()  * 2 - 1); 
      double p_phi = rand.Rndm() * M_PI + M_PI/2 + phi; //only consider directions towards inside. 
      double x = rin * cos(phi); 
      double y = rin * sin(phi); 
      double px = cos(p_phi) * sin(p_theta);
      double py = sin(p_phi) * sin(p_theta); 
      double pz = cos(p_theta); 

      pos[0] = x + x0;  
      pos[1] = y + y0;  
      pos[2] = z;  
      dir[0] = px; 
      dir[1] = py; 
      dir[2] = pz; 
      std::cout <<"inside: "<< p_theta <<std::endl;; 
      break;
    }
  }

}
