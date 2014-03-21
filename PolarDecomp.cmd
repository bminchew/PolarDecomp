
Input command file for PolarDecomp (must have :: for all input lines) 

Image parameters

   Columns (samples) in image                :: 8739
   Starting line                             :: 1
   Ending line (enter number, end, or all)   :: end 
   Starting column                           :: 1
   Ending column (enter number, end, or all) :: end     
           
   Window size (pixels per side)             :: 9 

   Input data folder (no / on the end) :: /home/uavops/Release/Icelnd_02003_09037_006_090610_L090_CX_01
   Output folder (no / on the end)     :: /u/mah-r7/bminchew/Icelnd_Pol/02003_09037_test
   Output file basename                :: snafu.grd
   
   Use UAVSAR file name convention?    :: yes 

   ************************************************************************
   *************  if not using UAVSAR file name convention  ***************
   *** HHHH filename                   :: f
   *** VVVV filename                   ::
   *** HVHV filename                   ::
   *** HHHV filename                   ::
   *** HHVV filename                   ::
   *** HVVV filename                   ::
   ************************************************************************
   ************************************************************************

   Output options ::   enter yes or no for each option

H/A/alpha parameters  

   Output entropy?                              (outputs *.ent)       :: yes
   Output anisotropy?                           (outputs *.ani)       :: yes
   Output alpha (degrees)?                      (outputs *.alp)       :: yes


Pauli decomposition elements (sqrt of diagonal components of T_3)

   Output Pauli_1 (|HH+VV|) value?              (outputs *.t11)       :: yes
   Output Pauli_2 (|HH-VV|) value?              (outputs *.t22)       :: yes
   Output Pauli_3 (2|HV|) value?                (outputs *.t33)       :: yes
   
Eigenvalues and related parameters

   Output largest eigenvalue  (lambda_1)?       (outputs *.ev1)       :: no     
   Output median eigenvalue   (lambda_2)?       (outputs *.ev2)       :: no
   Output smallest eigenvalue (lambda_3)?       (outputs *.ev3)       :: no

   Output pseudo probability for lambda_1?      (outputs *.p1)        :: no
   Output pseudo probability for lambda_2?      (outputs *.p2)        :: no
   Output pseudo probability for lambda_3?      (outputs *.p3)        :: no

   Output averaged intensity?                   (outputs *.lam)       :: no
      (averaged intensity: I = sum(PsuedoPrb_i*lambda_i){i=1-3} [dB])


Shannon entropy parameters

   Output total Shannon entropy?                (outputs *.se)        :: no
   Output Shannon polarimetric component?       (outputs *.sep)       :: no
   Output Shannon intensity component?          (outputs *.sei)       :: no

  
Coherency matrix T_3 (and covariance C_3) outputs

   Output trace(T_3)?                           (outputs *.mat.trt)   :: no
   Output det(T_3)?                             (outputs *.mat.det)   :: no
   Output I_2(T_3)?  (second invariant)         (outputs *.mat.iv2)   :: no


Other eigenvector parameters (weighted mean of each)

   Output Beta (degrees)?                       (outputs *.beta)      :: no
   Output gamma (degrees)?                      (outputs *.gamma)     :: no
   Output delta (degrees)?                      (outputs *.delta)     :: no
   Output phi (degrees)?                        (outputs *.phi)       :: no

Individual eigenvector parameters

   Output alpha_1 (degrees)?                    (outputs *.alp1)      :: no
   Output alpha_2 (degrees)?                    (outputs *.alp2)      :: no
   Output alpha_3 (degrees)?                    (outputs *.alp3)      :: no

   Output Beta_1 (degrees)?                    (outputs *.beta1)      :: no
   Output Beta_2 (degrees)?                    (outputs *.beta2)      :: no
   Output Beta_3 (degrees)?                    (outputs *.beta3)      :: no
   
Other options

   dB multiplier? (usually 10 or 20)                                 :: 10
       
