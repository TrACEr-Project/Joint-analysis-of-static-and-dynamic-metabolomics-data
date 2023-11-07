% This script shows an example to check the replicability of an ACMTF model
% applied to a matrix (static) and a tensor (dynamic)

% We use CMTF Toolbox version 1.1 as well as the Poblano Toolbox to fit the ACMTF model
% When calculating factor match score (FMS), we use score.m from the Tensor toolbox version 3.1
% parts of the scripts may require the dataset object (https://eigenvector.com/software/dataset-object/), publically available.

clear all
clc
close all


%% add auxilary functions to path
addpath(genpath('./functions'))

%% add pckages to path
addpath(genpath('.../tensor_toolbox-v3.1')) %Tensor toolbox is needed;  MATLAB Tensor Toolbox. Copyright 2017, Sandia Corporation, http://www.tensortoolbox.org/
addpath(genpath('.../CMTF_Toolbox_v1_1')) %CMTF toolbox is needed; Available online at http://www.models.life.ku.dk/joda/CMTF Toolbox
addpath(genpath('.../poblano_toolbox_1.1')) % Poblano toolbox is needed; Available online at https://github.com/sandialabs/poblano_toolbox
addpath(genpath('.../dataset')) % dataset object is needed; download here: https://github.com/sandialabs/poblano_toolbox


%%  load dataset
load('data.mat','X_T0c','X_T0') % X_T0c and X_T0 are dataset objects for dynamic and static data, respectively

%% data before split----rearrange the data so that Lower BMI subjects appear first
sub_normal=find(X_T0c.class{1,2}==1 | X_T0c.class{1,2}==4); % Lower BMI subjects
sub_abnormal=find(X_T0c.class{1,2}==2 | X_T0c.class{1,2}==3);% Higher BMI subjects
index_perm=[sub_normal,sub_abnormal];
X_T0c=X_T0c(index_perm,:,:);
X_T0=X_T0(index_perm,:,:);


%%  in total, we do N=10 rounds split10; kk-- the kk_th round;
N=10;
Results_all=cell(1,N*10);
for kk=1:N
    S_perm=[randperm(length(sub_normal)),randperm(length(sub_abnormal))+length(sub_normal)];
    X_split_left_T0c=cell(1,10);
    X_split_left_T0=cell(1,10);
    Results=cell(1,10);
    
    % data with random 1/10 subjects removed
    for i=1:10
        S_perm_rem=S_perm(i:10:end);
        pid_list=str2num(X_T0c.label{1});
        remove_pid=pid_list(S_perm_rem);
        X_split_left_T0c{i}=remove_outliers(X_T0c,remove_pid'); % the T0-corrected data after removing 1/10 part
        X_split_left_T0{i}=remove_outliers(X_T0,remove_pid'); % the T0 data after removing 1/10 part
        Results{i}.Sub_rem=S_perm_rem; % record the position of the removed data
    end
    
    % Preprocess the data and run ACMTF model for data from each split
    for ii=1:length(X_split_left_T0c)
        X_correct=X_split_left_T0c{i}.data;
        X_T0=squeeze(X_split_left_T0{i}.data);
        %%preprocess the data: centering and scaling
        Z1=tensor(preprocess(X_correct));
        Z2=tensor(preprocess(X_T0));
        %%ACMTFmodel---data preparation
        W{1} = tensor(~isnan(Z1.data));
        W{2} = tensor(~isnan(Z2.data));
        Z1(find(W{1}.data==0))=0;
        Z2(find(W{2}.data==0))=0;
        Z.object{1} = Z1;
        Z.object{2} = Z2;
        Z.modes     = {[1 2 3],[1 4]};
        Z.size      = [size(Z1) size(Z2,2)];
        for p=1:2
            norms(p)    = norm(Z.object{p});
            Z.object{p} = Z.object{p}/norms(p);
        end
        %%fit ACMTF-OPT
        options = ncg('defaults');
        options.Display ='final';
        options.MaxFuncEvals = 100000;
        options.MaxIters     = 10000;
        options.StopTol      = 1e-9;
        options.RelFuncTol   = 1e-9;
        R = 2;
        P = 2;
        beta     = [1e-3 1e-3]; % sparsity on the weihts of rank-one components
        nb_runs  = 32;
        goodness_X = zeros(nb_runs,1); %Stores
        for i=1:nb_runs
            if i==1 && R<min(size(Z))
                [Fac{i}, ~, out{i}]    = acmtf_opt(Z,R,'init','nvecs', 'alg_options', options,  'beta' ,beta, 'alg','ncg');
            else
                [Fac{i}, ~, out{i}]    = acmtf_opt(Z,R,'init','random', 'alg_options', options,  'beta' ,beta, 'alg','ncg');
            end
            goodness_X(i) = out{i}.F; %Stores objective function value
        end
        %%ACMTF result
        [ff, index] = sort(goodness_X,'ascend');
        Fac_sorted = Fac(index);
        out_sorted = out(index);
        %check the stopping condition
        if (out{index(1)}.ExitFlag~=2) && (out{index(1)}.ExitFlag~=1)
            flag_stop = true;
        else
            flag_stop = false;
        end
        Zhat = Fac_sorted{1};%To keep consistent with the bar code Evrim made
        l_rec = zeros(P, R);
        for p = 1:P
            temp        = normalize(Zhat{p});
            l_rec(p,:)  = temp.lambda;
        end
        for p=1:P
            data.Zhat{p}      = normalize(Zhat{p});
        end
        
        data.out        = out_sorted{1};
        data.func_eval  = ff;%value of the loss function
        
        data.Z          = Z;
        data.Fac_sorted=Fac_sorted;
        
        
        for i=1:2
            
            data.fit(i) = 100- norm(tensor(times(W{i},full(data.Zhat{i})-Z.object{i})))^2/norm(tensor(times(W{i},Z.object{i})))^2*100;%Z.object{i}
            
        end
        
        %%record the results by ACMTF models from all splits
        Results{1,ii}.cov_info=flag_stop;
        Results{1,ii}.info_best=data;
        Results{1,ii}.R=R;
        Results_all{1,(kk-1)*10+ii}=Results{1,ii};
    end
    
    
end

%% save results by the ACMTF model in each split
% save('Results_all.mat','Results_all')


%% compute FMS_tensor and FMS_matrix based on Results_all
% select the splits with unique factorizations by ACMTF model for further replicability check
% load('Results_all.mat','Results_all')
R=Results_all{1,1}.R;
% need to check uniquess within each split---check uniqueness for lambda,sigma, and factors from static and dynamic data
% store the index of unique splits in uni_index with 1 for unique split
uni_index_lamdas=zeros(length(Results_all),1);
uni_index_tensorFMS=zeros(length(Results_all),1);
uni_index_matrixsimi=zeros(length(Results_all),1);
uni_index=zeros(length(Results_all),1);
for i=1:length(Results_all)
    [Fac_aligned, lamdas, sigmas] = check_spread_only(R, Results_all{1,i}.info_best.Fac_sorted, Results_all{1,i}.info_best.func_eval, 0);
    % check lamdas and sigmas
    a=min(abs(max(lamdas)-min(lamdas))./abs(mean(lamdas)));
    b=min(abs(max(sigmas)-min(sigmas))./abs(mean(sigmas)));
    if length(Fac_aligned)<2
        disp (['need more starts in the',num2str(i),'subset'])
    else
        if max(a,b)<=1e-2
            uni_index_lamdas(i)=1;
        end
        % FMS
        for jj=1:length(Fac_aligned)
            for kk=1:length(Fac_aligned)
                FMS_score_tensor(jj,kk) = score(ktensor(Fac_aligned{jj}{1}.lambda,Fac_aligned{jj}{1}.U{2},Fac_aligned{jj}{1}.U{3}),...
                    ktensor(Fac_aligned{kk}{1}.lambda,Fac_aligned{kk}{1}.U{2},Fac_aligned{kk}{1}.U{3}),'lambda_penalty',false);
                Simi_score_matrix(jj,kk)=score(ktensor(Fac_aligned{jj}{2}.lambda,Fac_aligned{jj}{2}.U{2}),...
                    ktensor(Fac_aligned{kk}{2}.lambda,Fac_aligned{kk}{2}.U{2}),'lambda_penalty',false);
            end
        end
        if min(min(FMS_score_tensor))>=0.95
            uni_index_tensorFMS(i)=1;
        end
        if min(min(Simi_score_matrix))>=0.95
            uni_index_matrixsimi(i)=1;
        end
        if ( uni_index_lamdas(i)>0 & uni_index_tensorFMS(i)>0 ) & uni_index_matrixsimi(i)>0
            uni_index(i)=1;
        end
    end
    
    
end

%%pick only unique splits to check replicability
index_uniq=find(uni_index>0);
for i=1:length(index_uniq)
    Fac_all{i}=Results_all{1,index_uniq(i)}.info_best.Fac_sorted{1};
    ff_all(i)=Results_all{1,index_uniq(i)}.info_best.func_eval(1);
end
%%compte factor match score
[ff_sorted_all, index_ff] = sort(ff_all,'ascend');
for i=1:length(index_ff)
    Fac_sorted_all{i}=Fac_all{index_ff(i)};
end
kk=0;
for ii=1:length(Fac_sorted_all)
    for jj=ii+1:length(Fac_sorted_all)
        kk=kk+1;
        FMS_tensor(kk)=score(ktensor(Fac_sorted_all{ii}{1}.lambda,Fac_sorted_all{ii}{1}.U{2},Fac_sorted_all{ii}{1}.U{3}),...
            ktensor(Fac_sorted_all{jj}{1}.lambda,Fac_sorted_all{jj}{1}.U{2},Fac_sorted_all{jj}{1}.U{3}),'lambda_penalty',false);
        FMS_matrix(kk)=score(ktensor(Fac_sorted_all{ii}{2}.lambda,Fac_sorted_all{ii}{2}.U{2}),...
            ktensor(Fac_sorted_all{jj}{2}.lambda,Fac_sorted_all{jj}{2}.U{2}),'lambda_penalty',false);
    end
end



%% plot FMS_tensor and FMS_matrix
f=figure;
subplot(1,2,1)
boxplot(FMS_tensor); ylabel('FMS_{tensor}'); xlabel(['R=',num2str(R)])
set(gca,'Fontsize',18)
set(findobj(gca, 'type', 'line'), 'LineWidth', 1.5);
index= round(length(FMS_tensor)*0.95);
aaa = sort(FMS_tensor,'descend');
hold on;
plot(r-0.25:0.005:r+0.25, ones(101,1)*aaa(index),'g-','Linewidth',2);
yticks([0.2 0.4 0.6 0.7 0.8 0.9 1.0])
subplot(1,2,2)
boxplot(FMS_matrix); ylabel('FMS_{matrix}'); xlabel(['R=',num2str(R)])
set(findobj(gca, 'type', 'line'), 'LineWidth', 1.5);
index= round(length(FMS_matrix)*0.95);
aaa = sort(FMS_matrix,'descend');
hold on;
plot(r-0.25:0.005:r+0.25, ones(101,1)*aaa(index),'g-','Linewidth',2);
set(gca,'Fontsize',18)
yticks([0.2 0.4 0.6 0.7 0.8 0.9 1.0])








