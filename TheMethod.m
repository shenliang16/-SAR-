
%kps1: 初始匹配对中的特征集合1
%kps2: 初始匹配对中的特征集合2
%Dsc_N： 描述子匹配关系
%Dsc_Dis：描述子距离
%Dsc_Type：描述子类型
%thr_a：阈值参数
function [inliers_g, kp1, kp2, A, info] = TheMethod(kps1, kps2, Dsc_N, Dsc_Dis, Dsc_Type, thr_a, I1, I2)

    % 计算局部变换矩阵
    [ScalMatX, RotMatX,  A, invA, kp1_o, kp2_o, kp1, kp2, kd, cumsum_Kd, idx] = cal_transformation_matrix(kps1, kps2, Dsc_N);


    %% step1:  几何共识校验
    K1 = 8;
    K11 = K1+1;
    Amp = 1.2;
    Repj_Err = 1e5*ones(2, K1);
    N = size(kp1, 1);
    Scaling = kp2(:,3)./kp1(:,3);
    idx = cell(1,2);        idx0 = cell(1,2);
    idx{1} = 1:N; % 设置为1：N， 以保证Nei_sel_S可以直接使用
    idx{2} = 1:N;
    Ratio_test = Dsc_Dis{1}(:, 1)./Dsc_Dis{1}(:, 2);
    idx0{1} = Threshold_With_minNum(1, Ratio_test, 0.7, max(30, N/30), 0);
    idx0{2} = idx{1}(idx0{1}); 
    Thres = [0.5, 1.2];
    switch Dsc_Type
        case 'LAF'
            [neighborX2, neighborY2,   indX, indY,   Nei_sel_S, Nei_sel_L,  DisXt,DisYt] = ... 
                FindNeighbors2_Eclips_SingleSide_SR(kp1', kp2', idx, Scaling, A, invA, K1, Amp, idx0, 1, Thres); 

        case 'SR'
            [neighborX2, neighborY2,   indX, indY,   Nei_sel_S, Nei_sel_L,  DisXt,DisYt] = ... 
                FindNeighbors2_noRec2_SingleSIde_SR(kp1', kp2', idx, Scaling, K1, Amp, idx0, 1, Thres); 
    end 

    info.indX = indX;
    info.indY = indY;
    info.DisXt = DisXt;
    info.DisYt = DisYt;
    info.A = A;
    info.invA = invA;

            %[neighborX, neighborY,   indX, indY,   Nei_sel_S, Nei_sel_L,  DisXt,DisYt,  Xt01,Xt02,  Yt01,Yt02] = FindNeighbors2_Eclips(kp1, kp2, idx, Scaling, A, invA, K1, Amp, idx0);
            %K = 8;
            %[neighborX, DisXt] = knnsearch(kp1, kp1, 'K', K+1); 
            %neighborX = neighborX'; DisXt = DisXt';
    % 增加自身
    Nei_ori = Nei_sel_S;
    Nei_sel_S = [1:N; Nei_sel_S(1:K1,:)];
    [kp1_n, kp2_n, T1, T2] = normalizePoints(kp1(:, 1:2),  kp2(:, 1:2)); 
    Centr_Loc1 = kp1_n(Nei_sel_S(2:end, :), 1:2) - kp1_n(repmat(Nei_sel_S(1, :), K1, 1), 1:2);   
    Centr_Loc2 = kp2_n(Nei_sel_S(2:end, :), 1:2) - kp2_n(repmat(Nei_sel_S(1, :), K1, 1), 1:2);
    Centr_Loc1 = reshape(Centr_Loc1', 2, K1, []);
    Centr_Loc2 = reshape(Centr_Loc2', 2, K1, []);
    lamnda = 250;
    C_ngc = zeros(N,1);
    k=0;
    for i = 1 : size(A,1)
        for j = 1:kd{1}(i)
            k=k+1;
            %T_self = reshape(A(i,j,:), 2,2);
            for m = 1:K1
                T_nei  = reshape(A(Nei_ori(m,k),1,:), 2,2);
                Repj_Err(:, m) = Centr_Loc2(:,m,k) - (T2(1)/T1(1))*T_nei * Centr_Loc1(:,m,k);
            end
            C_ngc(k) = sum(exp(-lamnda*sum(Repj_Err.^2))); %/NumNei
%             if any(i == [74	75	149	338	339	347	348	352	392])
%                 pause(0.0001)
%             end
        end
    end

    %% 
    thr_g = max(min(multithresh(C_ngc,2), 1.8*thr_a), 0.8*thr_a);
    inliers_g = find(C_ngc>thr_g(1));
    %C_ngc([74	75	149	338	339	347	348	352	392])
%     kp1 = kp1(inliers_g1,:); 
%     kp2 = kp2(inliers_g1,:); 
%     kd{1} = kd{1}(inliers_g1,:);  
%     %kd{2} = kd{2}(inliers_g,:); 
%     cumsum_Kd = uint32(cumsum(kd{1}));
%     cumsum_Kd = uint32([[1; cumsum_Kd(1:end-1)+1],    cumsum_Kd]); 
%     A = A(inliers_g1,:,:); 
%     invA = invA(inliers_g1,:,:);
%     Ratio_test = Ratio_test(inliers_g1);
%     inliers_g = inliers_g1;

    %% step1:  几何共识校验
    K1 = 8;
    Amp = 1.2;
    N = size(kp1, 1);
    Scaling = kp2(:,3)./kp1(:,3);
    idx = cell(1,2);        idx0 = cell(1,2);
    idx{1} = 1:N; % 设置为1：N， 以保证Nei_sel_S可以直接使用
    idx{2} = 1:N;
    %idx0{1} = 1:N;
    idx0{1} = inliers_g;
    idx0{2} = idx{1}(idx0{1}); 
    Thres = [0.5, 1.2];
    switch Dsc_Type
        case 'LAF'
            [neighborX2, neighborY2,   indX, indY,   Nei_sel_S, Nei_sel_L] = ... 
                FindNeighbors2_Eclips_SingleSide_SR(kp1', kp2', idx, Scaling, A, invA, K1, Amp, idx0, 1, Thres); 

        case 'SR'
            [neighborX2, neighborY2,   indX, indY,   Nei_sel_S, Nei_sel_L] = ... 
                FindNeighbors2_noRec2_SingleSIde_SR(kp1', kp2', idx, Scaling, K1, Amp, idx0, 1, Thres); 
    end 
    % 增加自身
    Nei_ori = Nei_sel_S;
    K1 = min(K1, size(Nei_sel_S,1));
    Nei_sel_S = [1:N; Nei_sel_S(1:K1,:)];
    [kp1_n, kp2_n, T1, T2] = normalizePoints(kp1(:, 1:2),  kp2(:, 1:2)); 
    Centr_Loc1 = kp1_n(Nei_sel_S(2:end, :), 1:2) - kp1_n(repmat(Nei_sel_S(1, :), K1, 1), 1:2);   
    Centr_Loc2 = kp2_n(Nei_sel_S(2:end, :), 1:2) - kp2_n(repmat(Nei_sel_S(1, :), K1, 1), 1:2);
    Centr_Loc1 = reshape(Centr_Loc1', 2, K1, []);
    Centr_Loc2 = reshape(Centr_Loc2', 2, K1, []);
    lamnda = 250;
    C_ngc = zeros(N,1);
    k=0;
    for i = 1 : size(A,1)
        for j = 1:kd{1}(i)
            k=k+1;
            %T_self = reshape(A(i,j,:), 2,2);
            for m = 1:K1
                T_nei  = reshape(A(Nei_ori(m,k),1,:), 2,2);
                Repj_Err(:, m)  = Centr_Loc2(:,m,k) - (T2(1)/T1(1))*T_nei * Centr_Loc1(:,m,k);
            end
            C_ngc(k) = sum(exp(-lamnda*sum(Repj_Err.^2))); %/NumNei
%             if any(i == [74	75	149	338	339	347	348	352	392])%[5	74	75	149	338	339	352	392	478	480]
%                 pause(0.0001)
%             end
        end
    end

    %% 
    thr_g = max(min(multithresh(C_ngc,2), 2.0*thr_a), 1.5*thr_a);
    inliers_g1 = find(C_ngc>thr_g(1));
    kp1 = kp1(inliers_g1,:); 
    kp2 = kp2(inliers_g1,:); 
    kd{1} = kd{1}(inliers_g1,:);  
    %kd{2} = kd{2}(inliers_g,:); 
    cumsum_Kd = uint32(cumsum(kd{1}));
    cumsum_Kd = uint32([[1; cumsum_Kd(1:end-1)+1],    cumsum_Kd]); 
    A = A(inliers_g1,:,:); 
    invA = invA(inliers_g1,:,:);
    Ratio_test = Ratio_test(inliers_g1);
    inliers_g = inliers_g1;

    
    %% step1:  几何共识校验
    if length(inliers_g)>40
        K1 = 8;
        Amp = 1.2;
        N = size(kp1, 1);
        Scaling = kp2(:,3)./kp1(:,3);
        idx = cell(1,2);        idx0 = cell(1,2);
        idx{1} = 1:N; % 设置为1：N， 以保证Nei_sel_S可以直接使用
        idx{2} = 1:N;
        %idx0{1} = 1:N;
        idx0{1} = 1:N;
        idx0{2} = idx{1}(idx0{1}); 
        Thres = [0.5, 1.2];
        switch Dsc_Type
            case 'LAF'
                [neighborX2, neighborY2,   indX, indY,   Nei_sel_S, Nei_sel_L] = ... 
                    FindNeighbors2_Eclips_SingleSide_SR(kp1', kp2', idx, Scaling, A, invA, K1, Amp, idx0, 1, Thres); 
    
            case 'SR'
                [neighborX2, neighborY2,   indX, indY,   Nei_sel_S, Nei_sel_L] = ... 
                    FindNeighbors2_noRec2_SingleSIde_SR(kp1', kp2', idx, Scaling, K1, Amp, idx0, 1, Thres); 
        end 
        % 增加自身
        Nei_ori = Nei_sel_S;
        K1 = min(K1, size(Nei_sel_S,1));
        Nei_sel_S = [1:N; Nei_sel_S(1:K1,:)];
        [kp1_n, kp2_n, T1, T2] = normalizePoints(kp1(:, 1:2),  kp2(:, 1:2)); 
        Centr_Loc1 = kp1_n(Nei_sel_S(2:end, :), 1:2) - kp1_n(repmat(Nei_sel_S(1, :), K1, 1), 1:2);   
        Centr_Loc2 = kp2_n(Nei_sel_S(2:end, :), 1:2) - kp2_n(repmat(Nei_sel_S(1, :), K1, 1), 1:2);
        Centr_Loc1 = reshape(Centr_Loc1', 2, K1, []);
        Centr_Loc2 = reshape(Centr_Loc2', 2, K1, []);
        lamnda = 250;
        C_ngc = zeros(N,1);
        k=0;
        for i = 1 : size(A,1)
            for j = 1:kd{1}(i)
                k=k+1;
                %T_self = reshape(A(i,j,:), 2,2);
                for m = 1:K1
                    T_nei  = reshape(A(Nei_ori(m,k),1,:), 2,2);
                    Repj_Err(:, m)  = Centr_Loc2(:,m,k) - (T2(1)/T1(1))*T_nei * Centr_Loc1(:,m,k);
                end
                C_ngc(k) = sum(exp(-lamnda*sum(Repj_Err.^2))); %/NumNei
    %             if any(i == [5	149])
    %                 pause(0.1)
    %             end
            end
        end


        %% 
        thr_g = max(min(multithresh(C_ngc,2), 3.0*thr_a), 2.0*thr_a);
        inliers_g2 = find(C_ngc>thr_g(1));
        kp1 = kp1(inliers_g2,:); 
        kp2 = kp2(inliers_g2,:); 
        kd{1} = kd{1}(inliers_g2,:);  
        %kd{2} = kd{2}(inliers_g,:); 
        cumsum_Kd = uint32(cumsum(kd{1}));
        cumsum_Kd = uint32([[1; cumsum_Kd(1:end-1)+1],    cumsum_Kd]); 
        A = A(inliers_g2,:,:); 
        invA = invA(inliers_g2,:,:);
        inliers_g = inliers_g1(inliers_g2);
    end
end


















%%
function [ScalMatX, RotMatX,  TransforMatXY, TransforMatYX, X, Y, X_v, Y_v, Kd, cumsum_Kd, idx] = cal_transformation_matrix(X, Y, Dsc_N, Kdmax, Dsc_Dis, Kd_type, DscThr)
    Dsc_N_ori = Dsc_N;
    [N, ~]=size(X); 
    [N2, ~]=size(Y);
    if nargin<=4
        Kd_type = 'Static';
        Kdmax = 1;
    end
    Kd = cell(2,1);    Kd{1}=ones(N,1);    Kd{2}=ones(N,1);
    switch Kd_type
        case 'Adapti' % 根据欧氏距离，确定备选特征邻域
            % 修改
            Idx_K_Dsc = cell(2,1);
            Idx_K_Dsc{1} = Dsc_Dis{1}<DscThr(1)   &   Dsc_Dis{1}./repmat(Dsc_Dis{1}(:,1), 1, size(Dsc_Dis{1},2)) < DscThr(3);               
            Idx_K_Dsc{1}(:,1) = true;              
            Kd{1} = min(sum(Idx_K_Dsc{1}, 2), Kdmax);         
            Idx_K_Dsc{2} = Dsc_Dis{2}<DscThr(1)   &   Dsc_Dis{2}./repmat(Dsc_Dis{2}(:,1), 1, size(Dsc_Dis{2},2)) < DscThr(3); ;               
            Idx_K_Dsc{2}(:,1) = true;              
            Kd{2} = min(sum(Idx_K_Dsc{2}, 2), Kdmax);
            % 修改
            Dsc_Dis{1}(~Idx_K_Dsc{1}) = nan;
            Dsc_Dis{2}(~Idx_K_Dsc{2}) = nan;  
    
            for n = 1:N
                Dsc_N{1}(n,Kd{1}(n)+1:end) = nan;
            end
            for n = 1:N2
                Dsc_N{2}(n,Kd{2}(n)+1:end) = nan;
            end
        case 'Static' % 固定大小
            Kd{1} = Kdmax*ones(N,1);
            Kd{2} = Kdmax*ones(N2,1);
    end
    [ScalMatX, RotMatX,  TransforMatXY, TransforMatYX, X, Y] = SR_Matrix2...
                                                                (X, Y, Dsc_N, Dsc_N_ori, Kdmax);
    cumsum_Kd = uint32(cumsum(Kd{1}));
    cumsum_Kd = uint32([[1; cumsum_Kd(1:end-1)+1],    cumsum_Kd]);
    N_eff = cumsum_Kd(end, 2);
    X_v = zeros(N_eff, 4);             
    Y_v = zeros(N_eff, 4); 
    for n = 1:N
        idx = Dsc_N{1}(n, 1:Kd{1}(n));
        Y_v(cumsum_Kd(n, 1) : cumsum_Kd(n, 2), :) = Y(idx, :);
        X_v(cumsum_Kd(n, 1) : cumsum_Kd(n, 2), :) = repmat(X(n, :), Kd{1}(n), 1);
    end
end


function [ScalMatX, RotMatX,  TransforMatXY, TransforMatYX, X_ori, Y_ori] = SR_Matrix2(X_ori0, Y_ori0, Dsc_N, Dsc_N2, Kdmax)
    [N, D] = size(X_ori0);
    TransforMatXY = zeros(N, Kdmax, 4);
    TransforMatYX = zeros(N, Kdmax, 4);
    X_ori = Affine2Rigid(X_ori0);
    Y_ori = Affine2Rigid(Y_ori0);
    ScaleX = X_ori(:,3);        ScaXMat1 = repmat(ScaleX, 1, Kdmax);        
    ScaleY = Y_ori(:,3);        ScaYMat1 = ScaleY(Dsc_N2{1}(:,1:Kdmax));     
    RotX = X_ori(:,4);          RotXMat1 = repmat(RotX, 1, Kdmax);             
    RotY = Y_ori(:,4);          RotYMat1 = RotY(Dsc_N2{1}(:,1:Kdmax)); 
    ScalMatX = ScaYMat1./ScaXMat1;  
    RotMatX =  RotYMat1 - RotXMat1;
    NanX = isnan(Dsc_N{1});%% 去除Nan (针对自适应邻域数量，（邻域距离约束下的）)
    ScalMatX(NanX) = nan;  
    RotMatX(NanX) =  nan;
    switch D
        case 4
            %% 4D
            TransforMat_1 = ScalMatX.*cos(RotMatX);
            TransforMat_2 = ScalMatX.*sin(RotMatX);
            TransforMatXY(:, :, 1) = TransforMat_1;
            TransforMatXY(:, :, 4) = TransforMat_1;
            TransforMatXY(:, :, 2) = TransforMat_2;
            TransforMatXY(:, :, 3) = -TransforMat_2;
            TransforMat_1 = cos(RotMatX)./ScalMatX;
            TransforMat_2 = sin(RotMatX)./ScalMatX;
            TransforMatYX(:, :, 1) = TransforMat_1;
            TransforMatYX(:, :, 4) = TransforMat_1;
            TransforMatYX(:, :, 2) = -TransforMat_2;
            TransforMatYX(:, :, 3) = TransforMat_2;
        case 6
            %% 6D
            AnVec_X2Y = zeros(N, Kdmax, 4);  
            AnVec_Y2X = zeros(N, Kdmax, 4);
            for n = 1 : N
                P = reshape(X_ori0(n, 3:6), 2,2);
                for m = 1:Kdmax
                    idx = Dsc_N{1}(n, m);
                    if isnan(idx)
                        continue;
                    end
                    Q = reshape(Y_ori0(idx, 3:6), 2,2);   %取中心点X，Y矩阵
                    temp = Q(1:2,1:2)/P(1:2,1:2);
                    AnVec_X2Y(n, m, :) = temp(:)';
                    temp = P(1:2,1:2)/Q(1:2,1:2);
                    AnVec_Y2X(n, m, :) = temp(:)';
                end
            end
            TransforMatXY = AnVec_X2Y;
            TransforMatYX = AnVec_Y2X;
    end
end
