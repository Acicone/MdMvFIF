
 function [IMFs_s,IMFs_t,stats] = MdMvFIF_v5(f,options,M)

It generates the decomposition of the signal f :

 f = IMFs_s(:, :, :) + IMFs_s(:, :, :) + ... + IMFs_s(:, :, :) +
     IMFs_t(:, :, :) + IMFs_t(:, :, :) + ... + IMFs_t(:, :, :)

where the last row in the matrix IMF is the trend and the other rows
are actual IMFs

                               Inputs

  f         Tensor containing in the first two variables the space variability 
            and in the third one the time

  options    Structure, generated using function Settings_MdMvFIF_v1, containing
             all the parameters needed in the various algorithms

  options_FIF2  Structure, generated using function Settings_FIF2_v2, containing
             all the parameters needed in the various algorithms

  options_MvFIF Structure, generated using function Settings_FIF_v3, containing
             all the parameters needed in the various algorithms 


                              Output

  IMF_s and IMFs_t  Cells containing the IMFs decomposition in space and time

  stats     Statistics regarding the IMFs
              logM      Mask length values used for each IMF
              posF      position of the first minimum in the filter DFT which is forced to become zero
              valF      filter DFT first minimum value before the downward shift

  See also Settings_MdMvFIF_v1, Settings_FIF2_v2, SETTINGS_FIF_V3, mask_v1, Maxmins_v3_8, 
  FIF2_v3, MvFIF_v10.

 Please cite:

 A. Cicone, H. Zhou. 'Numerical Analysis for Iterative Filtering with 
 New Efficient Implementations Based on FFT'. Numerische Mathematik, 147 (1), pages 1-28, 2021. 
 doi: 10.1007/s00211-020-01165-5
 ArXiv http://arxiv.org/abs/1802.01359

 R. Cavassi, A. Cicone, E. Pellegrino, H. Zhou. 'A novel algorithm for the decomposition of non-stationary
 multidimensional and multivariate signals'. Submitted
