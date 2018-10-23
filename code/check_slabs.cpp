

#include"check_slabs.hpp"

/* ==========================================================================

                    FUNCTION  :  check_slab

=============================================================================

this function checks how many neighbor slabs in x,y,z directions need to be included
to fully cover the interaction in the FCs_cut_sym

input: FCs_cut_sym, LMs_ord, LMs_uc_ord, Lxyz_uc, NFold, comp_diff

output : (3,1)matrix nneigh:(num. slabs in x-dir, num. slabs in y-dir, num. slabs in z-dir, )
*/

Eigen::MatrixXd check_slabs(Eigen::MatrixXd FCs_cut_sym, Eigen::MatrixXd LMs_ord,
                            Eigen::MatrixXd LMs_uc_ord, Eigen::MatrixXd Lxyz_uc,
                            Eigen::MatrixXd NFold,double comp_diff )
{

  Eigen::MatrixXd nneigh(3,1); /* matrix to return
                                  Defines how many slabs need to be included 
                                  in the corresponding direction to cover 
                                  all interactions. dim(3,1) 
                                  (num. slabs in x-dir, num. slabs in y-dir, 
                                  num. slabs in z-dir, )*/
  
  // nneigh inizialization 
  nneigh(0,0)= 1; // 1 slab along x
  nneigh(1,0)= 1; // 1 slab along y
  nneigh(2,0)= 1; // 1 slab along z

// -------------- local variables --------------------

  int Na_uc = LMs_uc_ord.rows();

  int Na = LMs_ord.rows();

  int Hs_row = Na_uc*3*Na*3;


  int num_CELL = NFold(0)*NFold(1)*NFold(2); //number of cells repeated (which compose the supercell)


  Eigen::MatrixXd HS(Na_uc*3,Na_uc*3*num_CELL); /*matrix filled with FCs_cut_sym, following 
                                                  the repeated cell order
                                                  dim(3*Na_uc,3*Na_uc*num_CELL)
                                                  it is composed by num_CELL sub-matrices
                                                  each sub matrix is a (3*Na_uc,3*Na_uc) matrix
                                                  filled with the FCs_cut_sym elements corresponding
                                                  to the atoms in the unit cell and the n-th cell
                                                  remember that since FCs_cut_sym is the FC matrix cut after
                                                  a cut-off radius, some elements are 0. Thus some of the 
                                                  sub-matrices can be zeros matrix. in this case the number of slab
                                                  in that direction must be increased in order to fully cover the
                                                  interaction.
                                                The first 3*Na_uc columns define the sub-matrix of cell 1,
                                                The second 3*Na_uc columns define the sub-matrix of cell 2,
                                                and so on..
                                                since the loop is 
                                                for(ix=0,..,NFold-1){for(iy=0,..,NFoldy-1){for(iz=0,..,NFoldz-1)}}
                                                if NFold = (4,1,2), we have 
                                                (0,	0,	0) = cell 1
                                                (0,	0,	1) = cell 2
                                                (1,	0,	0) = cell 3
                                                (1,	0,	1) = cell 4
                                                (2,	0,	0) = cell 5
                                                (2,	0,	1) = cell 6
                                                (3,	0,	0) = cell 7
                                                (3,	0,	1) = cell 8
                                                */

  int IA1, IA2;

  int I =0;


//--- temporary LM matrixes are created --------------------------------------- 

/* LM_coord is  LMs_ord matrix, reduced to only coordinates 
   (rid off col 4 = atomic species, col 5 = original order) */    
 Eigen::MatrixXd LM_coord(Na,3); 

for(int i = 0; i < 3 ; i++)
  {
LM_coord.col(i) = LMs_ord.col(i);
  } 
            
/* LM_uc_coord is  LMs_uc_ord matrix, reduced to only coordinates 
   (rid off col 4 = atomic species, col 5 = original order)*/
Eigen::MatrixXd LM_uc_coord(Na_uc,3); 

for(int i = 0; i < 3 ; i++)
  {
LM_uc_coord.col(i) = LMs_uc_ord.col(i);
  }

//std::cout<<"\n\n"<< LM_coord << std::endl;
// std::cout<<"\n\n"<< LM_uc_coord << std::endl;


 Eigen::MatrixXd xyz2(1,3);
//-------------------------------------------------------------------------

  for(int IAuc1 = 0; IAuc1 < Na_uc; IAuc1++)
    {

      IA1 = find_vec(LM_uc_coord.row(IAuc1), LM_coord, comp_diff);


      //--building the supercell--------------     
      for(int IAuc2 = 0; IAuc2 < Na_uc; IAuc2++)
        { 
          for(int ix = 0; ix < NFold(0); ix++)
            {
              for(int iy = 0; iy < NFold(1); iy++)
                {
                  for(int iz =0; iz < NFold(2); iz++)
                    {                      

                      xyz2(0,0) = LM_uc_coord(IAuc2,0) + ix*Lxyz_uc(0);
                      xyz2(0,1) = LM_uc_coord(IAuc2,1) + iy*Lxyz_uc(1);
                      xyz2(0,2) = LM_uc_coord(IAuc2,2) + iz*Lxyz_uc(2);

                      IA2 = find_vec(xyz2,LM_coord, comp_diff);
                      //--- cell index I -----
                      if(I==num_CELL)
                        {
                          I=1;
                        }
                      else
                        {
                          I++; 
                        }
                      //--------Fill the HS matrix --------------------
                      
                      
                      for(int coord1=0; coord1 < 3; coord1++)
                        {
                          for(int coord2=0; coord2 < 3; coord2++)
                            {
                              HS(3*IAuc1+coord1, 3*IAuc2+coord2 + (Na_uc*3*(I-1)))= FCs_cut_sym(3*IA1+coord1,3*IA2+coord2);

                            }
                        }

                    }
                }
            }
        }
    }

  // ---- check slabs along z -------------------
  if(NFold(2)>1)
    {
      int startz=Na_uc*3;
      int endz= startz + Na_uc*3;
      int n=1;
      double tmpz=0.0;
      while(n < NFold(2))
        {
          for(int icz=startz; icz < endz; icz ++ )
            {
              tmpz+=HS.col(icz).sum();
            }
          startz = endz+1;
          endz = startz + Na_uc*3;
          n++;
          /* if all the elements of the sub-matrix are zeros 
             increase the NFold_z by 1 */
          if(tmpz==0.0) nneigh(2,0)+=1; 
            tmpz=0.0;// re-inizialize tmpz
        }
    }

  // ---- check slabs along y -------------------
  if(NFold(1)>1)  
    {
      int starty=Na_uc*3*NFold(2);
      int endy= starty + Na_uc*3;
      int ny=1;
      double tmpy=0.0;
      while(ny < NFold(1))
        {
          for(int icy=starty; icy < endy; icy ++ )
            {
              tmpy+=HS.col(icy).sum();
            }
          ny++;
          starty = Na_uc*3*NFold(2)*ny;
          endy = starty + Na_uc*3;
          /* if all the elements of the sub-matrix are zeros 
             increase the NFold_y by 1 */
          if(tmpy==0.0) nneigh(1,0)+=1; 
          tmpy=0.0;
        }
    }

  // ---- check slabs along x -------------------
  if(NFold(0)>1)  
    {
      int startx=Na_uc*3*NFold(1)*NFold(2);
      int endx= startx + Na_uc*3;
      int nx=1;
      double tmpx=0.0;
      while(nx < NFold(0))
        {
          for(int icx=startx; icx < endx; icx ++ )
            {
              tmpx+=HS.col(icx).sum();
            }
          nx++;
          startx = Na_uc*3*NFold(1)*NFold(2)*nx;
          endx = startx + Na_uc*3;
          /* if all the elements of the sub-matrix are zeros 
             increase the NFold_x by 1 */
          if(tmpx==0.0) nneigh(0,0)+=1; 
          tmpx=0.0;
        }
    }
 
  LM_coord.resize(0,0);
  LM_uc_coord.resize(0,0);

  return nneigh;
}
