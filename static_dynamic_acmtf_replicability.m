% This script shows how to check the replicability of an ACMTF model when jointly analyzing a matrix
% (static) and tensor (dynamic) metabolomics data: 1) randomly leave out 10% of the samples and compute
% the FMS between the best runs, 2) repeat this process ten times and return in total 450 FMS



% We use CMTF Toolbox as well as the Poblano Toolbox to fit the ACMTF model
% When calculating factor match score (FMS), we use score.m from the Tensor toolbox version 3.1
% Parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.

%% add auxilary functions to path
addpath(genpath('./functions'))
% add necessary toolboxes: Tensor Toolbox, CMTF Toolbox, Poblano toolbox and dataset object
% see https://github.com/Lu-source/Joint-analysis-of-static-and-dynamic-metabolomics-data/blob/main/README.md


%% load dataset
load('data.mat','X_T0c','X_T0') % X_T0c and X_T0 are dataset objects

%% data before split----rearrange the data so that Lower BMI subjects appear first
sub_normal   = find(X_T0c.class{1,2}==1 | X_T0c.class{1,2}==4); % Lower BMI subjects
sub_abnormal = find(X_T0c.class{1,2}==2 | X_T0c.class{1,2}==3);% Higher BMI subjects
index_perm   = [sub_normal,sub_abnormal];
X_T0c        = X_T0c(index_perm,:,:);
X_T0         = X_T0(index_perm,:);


%% ACMTF model---set up
options = ncg('defaults');
options.Display ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-8;
options.RelFuncTol   = 1e-8;
P = 2;
R = 2;
beta     = [1e-3 1e-3]; % sparsity on the weihts of rank-one components
nb_runs  = 32;


%% in total, we do N=10 rounds split10; kk--the kk_th round;
N = 10;
Results_all = cell(1,N);
for kk = 1:N
    % randomly shuffle the subjects and save the index in S_perm
    S_perm = [randperm(length(sub_normal)),randperm(length(sub_abnormal))+length(sub_normal)];
    X_split_left_T0c = cell(1,10);
    X_split_left_T0  = cell(1,10);
    Sub_rem_index    = cell(1,10);
    
    % data sets -- randomly remove 1/10 subjects in each data set
    for i = 1:10
        S_perm_rem = S_perm(i:10:end);
        pid_list   = str2num(X_T0c.label{1});
        remove_pid = pid_list(S_perm_rem);
        X_split_left_T0c{i} = removesubject(X_T0c,remove_pid'); % the T0-corrected data after removing 1/10 subjects
        X_split_left_T0{i}  = removesubject(X_T0,remove_pid'); % the T0 data after removing 1/10 subjects
        Sub_rem_index{1,i}   = S_perm_rem; % record the position of the removed subject
    end
    
    % preprocess the data and run ACMTF model for data from each split
    for ii = 1:length(X_split_left_T0c)
        X_correct = X_split_left_T0c{ii}.data;
        X_t0      = X_split_left_T0{ii}.data;
        % preprocess the data: centering and scaling
        Z1 = tensor(preprocess(X_correct));
        Z2 = tensor(preprocess(X_t0));
        % ACMTF model -- data preparation
        W{1} = tensor(~isnan(Z1.data));
        W{2} = tensor(~isnan(Z2.data));
        Z1(find(W{1}.data==0)) = 0;
        Z2(find(W{2}.data==0)) = 0;
        Z.miss{1} = W{1};
        Z.miss{2} = W{2};
        Z.object{1} = Z1;
        Z.object{2} = Z2;
        Z.modes     = {[1 2 3],[1 4]};
        Z.size      = [size(Z1) size(Z2,2)];
        for p=1:2
            norms(p)    = norm(Z.object{p});
            Z.object{p} = Z.object{p}/norms(p);
        end
        % fit ACMTF model
        for i = 1:nb_runs
            [Fac{i}, ~, out{i}]    = acmtf_opt(Z,R,'init','random', 'alg_options', options,  'beta' ,beta, 'alg','ncg');
            goodness_X(i) = out{i}.F; %Stores objective function value
        end
        % ACMTF result -- sorted
        [ff, index] = sort(goodness_X,'ascend');
        Fac_sorted = Fac(index);
        out_sorted = out(index);
        % check the stopping condition for the best run
        if (out{index(1)}.ExitFlag==0) | (out{index(1)}.ExitFlag==3)
            flag_stop = true;
        else
            flag_stop = false;
        end
        % normalize the best factor
        Zhat = Fac_sorted{1};
        for p = 1:P
            temp        = normalize(Zhat{p});
            data.Zhat{p}      = temp;
        end
        % save the sorted factors and function values of the ACMTF model
        data.Z          = Z;
        data.Fac_sorted=Fac_sorted;
        data.func_eval  = ff;%value of the loss function
        % calculate the model fit
        for i=1:2
            data.fit(i) = 100- norm(tensor(times(W{i},full(data.Zhat{i})-Z.object{i})))^2/norm(tensor(times(W{i},Z.object{i})))^2*100;
        end
        
        % record the factorization results from each split
        Results_ii.cov_info=flag_stop;
        Results_ii.info_best=data;
        Results_ii.R=R;
        Results_ii.Sub_rem=Sub_rem_index{1,ii};
        Results_all{1,kk}{1,ii}=Results_ii;
        % clear temporary variables
        clear Results_ii data Z1 Z2 Z W Zhat ...
            l_rec temp Fac_sorted out_sorted ff index
        
    end
    
    clear Sub_rem_index X_split_left_T0c X_split_left_T0
    
end


%% compute FMS_tensor and FMS_matrix: within each split, compute pairwise FMS from unique factorizations
FMS_tensor=[];
FMS_matrix=[];
for kkk=1:length(Results_all) % Loop through all splits, in total 10 splits 
    l=0; % record the number of unique factorization in each split
    for i=1:length(Results_all{1,kkk})
        Fac_aligned = show_spread(R, Results_all{1,kkk}{1,i}.info_best.Fac_sorted, Results_all{1,kkk}{1,i}.info_best.func_eval,1);
        if length(Fac_aligned)<2
            disp (['need more starts in the',num2str(i),'subset and ',num2str(kkk), 'split'])
        else
            ll=0;
            for jj = 1:length(Fac_aligned)
                for kk = jj+1:length(Fac_aligned)
                    ll=ll+1;
                    FMS_score_tensor(ll) = score(Fac_aligned{jj}{1},Fac_aligned{kk}{1},'lambda_penalty',false);
                    FMS_score_matrix(ll)=score(Fac_aligned{jj}{2},Fac_aligned{kk}{2},'lambda_penalty',false);
                end
            end
            if min(FMS_score_tensor)>=0.95
                uni_index_tensorFMS = 1; % factorization with unique factors--tensor
            end
            if min(FMS_score_matrix)>=0.95
                uni_index_matrixFMS = 1; % factorization with unique factors--matrix
            end
            if  uni_index_tensorFMS>0  & uni_index_matrixFMS>0
                l = l+1; % record the number of unique factorization in the kkk_th split
                Results_kkksplit_temp{1,l}=Results_all{1,kkk}{1,i}; % temp variable to store unique factorizations in the kkk_th split
            end
        end
    end
    % compute pairwise FMS in the kkk_th split 
    ll=0;
    for jj=1:length(Results_kkksplit_temp)
        for kk=jj+1:length(Results_kkksplit_temp)
            ll=ll+1;
            Fac_real_temp_jj=Results_kkksplit_temp{1,jj}.info_best.Fac_sorted{1}{1};
            Fac_real_temp_kk=Results_kkksplit_temp{1,kk}.info_best.Fac_sorted{1}{1};
            fms_real_kkk_split(ll,1)=score(ktensor(Fac_real_temp_jj.lambda,Fac_real_temp_jj.U{2},Fac_real_temp_jj.U{3}),...
                ktensor(Fac_real_temp_kk.lambda,Fac_real_temp_kk.U{2},Fac_real_temp_kk.U{3}),'lambda_penalty',false);
            Fac_sim_temp_jj=Results_kkksplit_temp{1,jj}.info_best.Fac_sorted{1}{2};
            Fac_sim_temp_kk=Results_kkksplit_temp{1,kk}.info_best.Fac_sorted{1}{2};
            fms_sim_kkk_split(ll,1)=score(ktensor(Fac_sim_temp_jj.lambda,Fac_sim_temp_jj.U{2}),...
                ktensor(Fac_sim_temp_kk.lambda,Fac_sim_temp_kk.U{2}),'lambda_penalty',false);
        end
    end
    % save all pairwise FMS from all splits in FMS_tensor and FMS_matrix
    FMS_tensor=[FMS_tensor;fms_real_kkk_split];
    FMS_matrix=[FMS_matrix;fms_sim_kkk_split];
    clear    fms_real_kkk_split  fms_sim_kkk_split Results_kkksplit_temp
end











