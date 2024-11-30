function [IMFs_s,IMFs_t,stats] = MdMvFIF_v5(f,options,options_FIF2,options_MvFIF)
%
%  function [IMFs_s,IMFs_t,stats] = MdMvFIF_v5(f,options,M)
%
% It generates the decomposition of the signal f :
%
%  f = IMFs_s(:, :, :) + IMFs_s(:, :, :) + ... + IMFs_s(:, :, :) +
%      IMFs_t(:, :, :) + IMFs_t(:, :, :) + ... + IMFs_t(:, :, :)
%
% where the last row in the matrix IMF is the trend and the other rows
% are actual IMFs
%
%                                Inputs
%
%   f         Tensor containing in the first two variables the space variability 
%             and in the third one the time
%
%   options    Structure, generated using function Settings_MdMvFIF_v1, containing
%              all the parameters needed in the various algorithms
%
%   options_FIF2  Structure, generated using function Settings_FIF2_v2, containing
%              all the parameters needed in the various algorithms
%
%   options_MvFIF Structure, generated using function Settings_FIF_v3, containing
%              all the parameters needed in the various algorithms 
%
%
%                               Output
%
%   IMF_s and IMFs_t  Cells containing the IMFs decomposition in space and time
%
%   stats     Statistics regarding the IMFs
%               logM      Mask length values used for each IMF
%               posF      position of the first minimum in the filter DFT which is forced to become zero
%               valF      filter DFT first minimum value before the downward shift
%
%   See also Settings_MdMvFIF_v1, Settings_FIF2_v2, SETTINGS_FIF_V3, mask_v1, Maxmins_v3_8, 
%   FIF2_v3, MvFIF_v10.
%
%  Please cite:
%
%  A. Cicone, H. Zhou. 'Numerical Analysis for Iterative Filtering with 
%  New Efficient Implementations Based on FFT'. Numerische Mathematik, 147 (1), pages 1-28, 2021. 
%  doi: 10.1007/s00211-020-01165-5
%  ArXiv http://arxiv.org/abs/1802.01359
%
%  R. Cavassi, A. Cicone, E. Pellegrino, H. Zhou. 'A novel algorithm for the decomposition of non-stationary
%  multidimensional and multivariate signals'. Submitted


if nargin == 0,  help MdMvFIF_v5; IMFs_s=[];IMFs_t=[];stats=[];return; end
if nargin == 1, options = Settings_MdMvFIF_v1; options_FIF2  = Settings_FIF2_v2; options_MvFIF = Settings_FIF_v3;  end
if nargin == 2, options_FIF2  = Settings_FIF2_v2; options_MvFIF = Settings_FIF_v3;  end
if nargin == 3, options_MvFIF = Settings_FIF_v3;  end

% we initiate the outer loop
nIMFs=0;
stats=[];
IMFs_s=cell(1,options.NIMFs+1);
IMFs_t=cell(1,options.NIMFs+1);

while nIMFs < options.NIMFs
    nIMFs=nIMFs+1
    if options.verbose>0
        fprintf('\n IMF # %1.0d\n',nIMFs)
    end
    SD=1;
    % we compute the average space mask length based on h_n
    if length(options_FIF2.alpha)>1
        mask_space = mask_v1(f(:),options_FIF2.alpha(nIMFs),options_FIF2.Xi);
    else
        mask_space = mask_v1(f(:),options_FIF2.alpha,options_FIF2.Xi);
    end
    stats(nIMFs).logM_space=mask_space;
    % we build a temporary signal which contains in a single vector all the
    % column of the original signal
    h_n=reshape(f,size(f,1)*size(f,2),size(f,3));
    h_pp=real(acos(dot(normc(h_n(1:end-1,:)),normc(h_n(2:end,:)))));
    % we compute the average space mask length based on h_pp
    if length(options_MvFIF.alpha)>1
        stats_time= mask_v1(h_pp,options_MvFIF.alpha(nIMFs),options_MvFIF.Xi);
    else
        stats_time= mask_v1(h_pp,options_MvFIF.alpha,options_MvFIF.Xi);
    end
    stats(nIMFs).logM_time=stats_time;
    inStepN=0;
    if nIMFs>1
        if stats(nIMFs).logM_space<=stats(nIMFs-1).logM_space
            if options.verbose>0
                fprintf('Warning mask length along space is decreasing at step %1d. ',nIMFs)
            end
            if options.MonotoneMaskLength==true
                stats(nIMFs).logM_space=ceil(stats(nIMFs-1).logM_space*1.1);
                if options.verbose>0
                    fprintf('The old mask length along space was %1d whereas the new one is forced to be %1d.\n',stats(nIMFs-1).logM_space,stats(nIMFs).logM_space)
                end
            else
                if options.verbose>0
                    fprintf('The old mask length along space was %1d whereas the new one is %1d.\n',stats(nIMFs-1).logM_space,stats(nIMFs).logM_space)
                end
            end
        end
        if stats(nIMFs).logM_time<=stats(nIMFs-1).logM_time
            if options.verbose>0
                fprintf('Warning mask length along time is decreasing at step %1d. ',nIMFs)
            end
            if options.MonotoneMaskLength==true
                stats(nIMFs).logM_time=ceil(stats(nIMFs-1).logM_time*1.1);
                if options.verbose>0
                    fprintf('The old mask length along time was %1d whereas the new one is forced to be %1d.\n',stats(nIMFs-1).logM_time,stats(nIMFs).logM_time)
                end
            else
                if options.verbose>0
                    fprintf('The old mask length along time was %1d whereas the new one is %1d.\n',stats(nIMFs-1).logM_time,stats(nIMFs).logM_time)
                end
            end
        end
    end

    %% IMF calculation
    % we start decomposing over space one step
    IMFs_s{nIMFs}=zeros(size(f));
    IMF_s_trend=zeros(size(f));
    for i=1:size(f,3)
        temp_IMF = FIF2_v3(f(:,:,i),options_FIF2,stats(nIMFs).logM_space);
        if options_FIF2.NIMFs==1
            IMFs_s{nIMFs}(:,:,i)=temp_IMF(:,:,1);
            IMF_s_trend(:,:,i)=temp_IMF(:,:,2);
        else
            disp('Warning: this code can work only when FIF2 extract 1 IMF at the time.')
            disp('Please, set NIMFs to 1 in Settings_FIF2_v2')
        end
    end
    % now we do one step of time decomposition
    temp_IMF=zeros(size(f,3),size(f,1)*size(f,2));
    for i=1:size(f,1)
        for j=1:size(f,2)
            temp_IMF(:,(i-1)*size(f,2)+j)=reshape(IMF_s_trend(i,j,:),[size(f,3),1]);
        end
    end
    temp_IMF2 = MvFIF_v10(temp_IMF,options_MvFIF,stats(nIMFs).logM_time);


    for i=1:size(f,1)
        for j=1:size(f,2)
            IMFs_t{nIMFs}(i,j,:)=temp_IMF2{(i-1)*size(f,2)+j}(:,1);
        end
    end
    f=f-IMFs_s{nIMFs}-IMFs_t{nIMFs}; % trend after computing IMF # nIMFs

end

%% We store the final trend in the space decomposition and clean the IMFs tensors

IMFs_s{nIMFs+1} = f;
IMFs_t{nIMFs+1} = zeros(size(f));
IMFs_t(nIMFs+2:end)=[];
IMFs_s(nIMFs+2:end)=[];

end