

#include"extract_submatrices.hpp"

/* ==========================================================================

                    FUNCTION  :  extract_submatrices

=============================================================================

this function calculates  matrices in OMEN style, 
using the updated structure information from check_slabs function 
and extend uc section 
(new unit cell = Na_uc_new, Lxyz_uc_new, LMs_uc_new)

Note that we are studying transport properties -> 
at least trasport direction (x) must be repeated = NFold_bs(0)=3

input: FCs_cut_sym, LMs_ord, LMs_uc_new, Lxyz,Lxyz_uc_new, NFold_bs

output : Hs_mod Tensor (NFold_bs_x, NFold_bs_y*3*Na_uc_new, NFold_bs_z*3*Na_uc_new) 
         
         it is composed by
         NFold_bs(0) matrices = 3 matrices: H00,H10,H01 
         matrices dimensions: ( NFold_bs_y*3*Na_uc_new, NFold_bs_z*3*Na_uc_new)
         elements : FCs_cut_sym

Note that:
1) 
NFold_bs=[3,3,3] = 3D -> 3 matrices of dim (3*Na_uc_new,3*3*3*Na_uc_new)
NFold_bs=[3,1,3] = 2D -> 3 matrices of dim (3*Na_uc_new,1*3*3*Na_uc_new)
NFold_bs=[3,1,1] = 1D -> 3 matrices of dim (3*Na_uc_new,1*1*3*Na_uc_new)

2)
to fill the matrices, artificially repeated new unit cells are created:
NFold_bs(0) new unit cells along x
NFold_bs(1) new unit cells along y
NFold_bs(2) new unit cells along z
Clearly the numbers of atoms artificially created is > Na
(Remember we only have FCs_cut_sym elements for the original
Na atoms, thus not all the artificially generated atoms 
will be used.)
*/

Eigen::Tensor<double, 3> extract_submatrices(Eigen::MatrixXd FCs_cut_sym,int Na, int Na_uc_new,
                                      Eigen::MatrixXd LMs_ord, Eigen::MatrixXd LMs_uc_new, 
                                      Eigen::MatrixXd Lxyz,Eigen::MatrixXd Lxyz_uc_new,Eigen::MatrixXd NFold_bs)
{

  Eigen::Tensor<double, 3> Hs_mod(NFold_bs(0), 3*Na_uc_new*NFold_bs(1), NFold_bs(2)*3*Na_uc_new);
Hs_mod.setZero();


//----- local variables -----------

 int iLx, iLy, iLz,iz,iy;

/* LM_coord is  LMs_ord matrix, reduced to only coordinates 
   (rid off col 4 = atomic species, col 5 = original order) */    
 Eigen::MatrixXd LM_coord(Na,3); 

for(int i = 0; i < 3 ; i++)
  {
LM_coord.col(i) = LMs_ord.col(i);
  } 
            

 Eigen::MatrixXd xyz1(1,3);

 Eigen::MatrixXd xyz2(1,3);

 xyz1.setZero(); // inizialization
 xyz2.setZero(); // inizialitazion

 double comp_diff = 5*1e-5;

 int IA1, IA2;
//-------------------------------------------------------------------------


for(int ixt =0; ixt < NFold_bs(0); ixt ++)
  { 
     if(ixt==0) iLx = 1;
     else iLx = 0;
     for(int iyt =0; iyt < NFold_bs(1); iyt ++)
       {
         if(NFold_bs(1)==1) 
           {
             iLy =0;
             iyt =1;
             iy =0;
           }
         else
           
           { 
             iy =iyt;
             if(iyt==0) iLy = 1;
             else iLy = 0;
           }
         for(int izt =0; izt < NFold_bs(2); izt ++)
           {
             if(NFold_bs(2)==1) 
               {
                 iLz =0;
                 izt =1;
                 iz =0;
               }
             else
               
               { 
                 iz = izt;
                 if(izt==0) iLz = 1;
                 else iLz = 0;
               }

             

             /*Starting from the atoms in the new uc,
               N  atoms are created 
               (N = NFold_bs(0)*NFold_bs(1)*NFold_bs(2)*Na_uc_new)
               xyz_a(j) = LMs_uc_new(IAnuc1,j)+iLj*Lxyz_uc_new(j) 
               where a = ..,N   and j = x,y,z 
               above them, only those which are within the sc, 
               are searched in the LMs_ord matrix 
               (the corresponding index is denoted by IA1)*/
             for(int IAnuc1 =0; IAnuc1 < Na_uc_new; IAnuc1 ++)
               {
                 xyz1(0,0) = LMs_uc_new(IAnuc1,0)+ iLx*Lxyz_uc_new(0);
                 xyz1(0,1) = LMs_uc_new(IAnuc1,1)+ iLy*Lxyz_uc_new(1);
                 xyz1(0,2) = LMs_uc_new(IAnuc1,2)+ iLz*Lxyz_uc_new(2);
                 // above we created the N atoms
                
                 if(xyz1(0,0)< (Lxyz(0,0)-comp_diff))
                   {
                      if(xyz1(0,1)< (Lxyz(0,1)-comp_diff))
                        {
                           if(xyz1(0,2)< (Lxyz(0,2)-comp_diff))
                             {
                          
                               IA1 = find_vec(xyz1,LM_coord,comp_diff);
                             }
                        }
                   }
         

                 /*Starting from the atoms in the new uc,
                   M  atoms are created 
                   (M =Na_uc_new*(NFold_bs(0)*NFold_bs(1)*NFold_bs(2))^2)
                   xyz2_b(j) = LMs_uc_new(IAnuc1,j)+iLj*Lxyz_uc_new(j) 
                   where b = ..,M   and j = x,y,z 
                   above them, only those which are within the sc, 
                   are searched in the LMs_ord matrix 
                   (the corresponding index is denoted by IA2)*/

                 for(int IAnuc2 =0; IAnuc2 < Na_uc_new; IAnuc2 ++)
                   {
                     xyz2(0,0) = LMs_uc_new(IAnuc2,0)+(ixt-1)*Lxyz_uc_new(0) +iLx*Lxyz_uc_new(0);
                     xyz2(0,1) = LMs_uc_new(IAnuc2,1)+(iyt-1)*Lxyz_uc_new(1)+ iLy*Lxyz_uc_new(1);
                     xyz2(0,2) = LMs_uc_new(IAnuc2,2)+(izt-1)*Lxyz_uc_new(2)+ iLz*Lxyz_uc_new(2);
                     // above we created the M atoms
                     
                     if(xyz2(0,0)< (Lxyz(0,0)-comp_diff))
                       {
                         if(xyz2(0,1)< (Lxyz(0,1)-comp_diff))
                           {
                             if(xyz2(0,2)< (Lxyz(0,2)-comp_diff))
                               {
                                 IA2 = find_vec(xyz2,LM_coord,comp_diff);
                                
                               }
                           }
                       }
                 

                     // fill Hs_mod tensor
                      for(int r=0; r < 3; r++)
                       {
                         for(int c=0; c < 3; c++)
                           {
                              Hs_mod(ixt, 3*IAnuc1+r + 3*Na_uc_new*iy, 3*Na_uc_new*iz + 3*IAnuc2+c)=FCs_cut_sym(3*IA1+r,3*IA2+c);
                           

                           }
                       }


                   }
                 
               }
           }
       }
 
  }




  return Hs_mod;

  Hs_mod.setZero();
}
