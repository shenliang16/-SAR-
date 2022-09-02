function [neighborX, neighborY,   indX, indY,   Nei_sel_S, Nei_sel_L,  DisXt_ori,DisYt_ori,  Xt01,Xt02,  Yt01,Yt02,  neighborS_ori, neighborX_ori, neighborY_ori] = FindNeighbors2_Eclips_SingleSide...
                                                           (Xt, Yt, idx, Scaling, TransforMatXY_vec, TransforMatYX_vec, K1, Amp, idx0, Scale_thr_side_selection, Thres)
%   Scaling:    尺度
%   Amp:        邻域放大比例
%   K1：        基础邻域数量
%   idx：       首对应序号
%   idx0：      SIFT_ratio 对应序号

    Rotation = Yt(4,:) - Xt(4,:);
    Xt = Xt(1:2,:);
    Yt = Yt(1:2,:);
    %% 去除重复点（由多对一的匹配引起），这在大缩放下尤为重要
    N = size(Xt,2);
    XtN = Xt(1:2, idx{1}(idx0{1}));
    YtN = Yt(1:2, idx{2}(idx0{2}));
    [~, iaX, icX] = unique(XtN', 'rows');
    [~, iaY, icY] = unique(YtN', 'rows');
    iX_sel = idx{1}(idx0{1}(iaX));
    iY_sel = idx{2}(idx0{2}(iaY)); 
    Xt1 = Xt(:, iX_sel);
    Yt1 = Yt(:, iY_sel);
    s_Xt1 = Scaling(iX_sel);
    s_Yt1 = Scaling(iY_sel);
    r_Xt1 = Rotation(iX_sel);
    r_Yt1 = Rotation(iY_sel);

    K1_X = min([length(iaX)-1, K1]);
    K1_1_X = K1_X+1;
    K1_Y = min([length(iaY)-1, K1]);
    K1_1_Y = K1_Y+1;
    indX = Scaling>Scale_thr_side_selection;
    Ind_left = find(indX);
    indY = ~indX;
    Ind_right = find(indY);                     % 满足以右侧为基准的序号
    %% KNN
        Kmax_X = min([3*K1_X, length(iaX)-1]);%Kmax=round(max(min(100,(N-1)/2),K1));%ff=(ff-min(ff))./(max(ff)-min(ff));
        Kmax_1_X = Kmax_X+1;
        Kmax_Y = min([3*K1_Y, length(iaY)-1]);%Kmax=round(max(min(100,(N-1)/2),K1));%ff=(ff-min(ff))./(max(ff)-min(ff));
        Kmax_1_Y = Kmax_Y+1;

        [neighborX, DisXt] = knnsearch(Xt1', Xt', 'K', Kmax_1_X); neighborX= neighborX'; DisXt = DisXt'; 
        [neighborY, DisYt] = knnsearch(Yt1', Yt', 'K', Kmax_1_Y); neighborY= neighborY'; DisYt = DisYt';

        DisXt_ori = DisXt;
        DisYt_ori = DisYt;
        %% 预筛选
            s_Xt1 = s_Xt1(neighborX); 
            r_Xt1 = r_Xt1(neighborX);
            s_Yt1 = s_Yt1(neighborY);
            r_Yt1 = r_Yt1(neighborY);
            %% 排序
%             sx = abs(log(s_Xt1 ./ s_Xt1(1,:)));
%             rx = abs(wrapToPi(r_Xt1 - r_Xt1(1,:)));
%             sy = abs(log(s_Yt1 ./ s_Yt1(1,:)));
%             ry = abs(wrapToPi(r_Yt1 - r_Yt1(1,:)));
%             sort(sx)
%             sort(sx)
            %% flag
            thr_s = Thres(1);
            thr_r = Thres(2);
            flag_sx = log(s_Xt1 ./ Scaling') < thr_s;
            flag_rx = abs(wrapToPi(r_Xt1 - Rotation)) < thr_r;
            flag_sy = log(s_Yt1 ./ Scaling') < thr_s;
            flag_ry = abs(wrapToPi(r_Yt1 - Rotation)) < thr_r;
            
            flag_x = flag_sx & flag_rx;
            flag_y = flag_sy & flag_ry;
            for loop_i = 1:length(Ind_left)
               i = Ind_left(loop_i);
               idx_i_1 = find(flag_x(:,i));
               idx_i_2 = find(~flag_x(:,i));
               neighborX(:, i) = neighborX([idx_i_1; idx_i_2], i);

               DisXt_i_1 = DisXt(idx_i_1, i);
               DisXt_i_2 = DisXt(idx_i_2, i);
               if ~isempty(DisXt_i_1)
                    DisXt_i_2 = max(DisXt_i_2, DisXt_i_1(end));    %重排序后，距离不递增了
               end
               DisXt(:, i) = [DisXt_i_1; DisXt_i_2];
            end
            for loop_i = 1:length(Ind_right)
               i = Ind_right(loop_i);
               idx_i_1 = find(flag_y(:,i)); 
               idx_i_2 = find(~flag_y(:,i));
               neighborY(:, i) = neighborY([idx_i_1; idx_i_2], i);
               DisYt(:, i) = DisYt([idx_i_1; idx_i_2], i);

               DisYt_i_1 = DisYt(idx_i_1, i);
               DisYt_i_2 = DisYt(idx_i_2, i);
               if ~isempty(DisYt_i_1)
                    DisYt_i_2 = max(DisYt_i_2, DisYt_i_1(end));    %重排序后，距离不递增了
               end
               DisYt(:, i) = [DisYt_i_1; DisYt_i_2];
            end

     %% 邻居的 中心距离/角度 计算
        neighborX2 = neighborX;
        neighborY2 = neighborY;
        Xt0 = [Xt1'; [0,0]];     Xt01=Xt0(:,1);     Xt02=Xt0(:,2);
        Yt0 = [Yt1'; [0,0]];     Yt01=Yt0(:,1);     Yt02=Yt0(:,2);  
        Xt01 = Xt01(neighborX2); 
        Xt02 = Xt02(neighborX2);
        Yt01 = Yt01(neighborY2);
        Yt02 = Yt02(neighborY2);
        Xt01 = Xt01-Xt(1,:);
        Xt02 = Xt02-Xt(2,:);
        Yt01 = Yt01-Yt(1,:);
        Yt02 = Yt02-Yt(2,:);

    %% 恢复序号
        neighborX = iaX(neighborX);
        neighborY = iaY(neighborY);
        neighborX = idx0{1}(neighborX);
        neighborY = idx0{2}(neighborY);

     %% 去掉自身   
        % 修改
            idx_del_X = DisXt(1,:) == 0;
            idx_del_Y = DisYt(1,:) == 0;
            Xt01(1:Kmax_X,idx_del_X)      = Xt01(2:Kmax_1_X,idx_del_X);
            Xt02(1:Kmax_X,idx_del_X)      = Xt02(2:Kmax_1_X,idx_del_X);
            DisXt(1:Kmax_X,idx_del_X)     = DisXt(2:Kmax_1_X,idx_del_X);
            neighborX(1:Kmax_X,idx_del_X) = neighborX(2:Kmax_1_X,idx_del_X);
            
            Yt01(1:Kmax_Y,idx_del_Y)      = Yt01(2:Kmax_1_Y,idx_del_Y);
            Yt02(1:Kmax_Y,idx_del_Y)      = Yt02(2:Kmax_1_Y,idx_del_Y);
            DisYt(1:Kmax_Y,idx_del_Y)     = DisYt(2:Kmax_1_Y,idx_del_Y);
            neighborY(1:Kmax_Y,idx_del_Y) = neighborY(2:Kmax_1_Y,idx_del_Y);

        neighborX_ori = neighborX;
        neighborY_ori = neighborY;

    %% 比例阈值 和 最小数量
        num_min_X = min(K1_X+1, Kmax_X);    
        num_min_Y = min(K1_Y+1, Kmax_Y);

    %% 以左侧Xt为基准
        N_left = length(Ind_left);
        Y_nei = neighborY(2:end, Ind_left);
        Nei_Y_sel = zeros(size(Y_nei)); 
        neighborY(K1_Y+2:end, :) = 0;
            Sel_FlagY = false(size(Nei_Y_sel));
            Y01_ = Yt01(2:end, Ind_left);
            Y02_ = Yt02(2:end, Ind_left);
            ThrDisY = Amp * DisXt(K1_1_X, Ind_left);
            for i = 1 : N_left
                n = Ind_left(i);
                inv_Txy = reshape(TransforMatYX_vec(n, :), 2,2);
                Yi = [Y01_(:,i), Y02_(:,i)];
                for j = 1 : Kmax_Y
                    Dis = inv_Txy * Yi(j, :)';
                    Dis = sum(Dis.^2);
                    Sel_FlagY(j, i) = Dis < ThrDisY(i);
                    if ~Sel_FlagY(j, i)
                        if Dis > 2*ThrDisY(i)
                            break;
                        end
                    end
                end
            end
        Sel_FlagY(1:num_min_Y, :) = true;
            for i = 1:N_left
                temp_i = Y_nei(Sel_FlagY(:, i), i);
                N_i = numel(temp_i);
                Nei_Y_sel(1:N_i, i) = temp_i;

                flags = [true; Sel_FlagY(:, i)];
                n = Ind_left(i);
                neighborY_ori(:, n) = [neighborY_ori(flags, n); neighborY_ori(~flags, n)];
                DisYt(:, n)         = [DisYt(flags, n);         DisYt(~flags, n)];
                Yt01(:, n)          = [Yt01(flags, n);          Yt01(~flags, n)];
                Yt02(:, n)          = [Yt02(flags, n);          Yt02(~flags, n)];
            end
        neighborY(2:end, Ind_left) = Nei_Y_sel;

        
    %% 以又侧Yt为基准
        N_right = length(Ind_right);
        X_nei = neighborX(2:end, Ind_right);
        Nei_X_sel = zeros(size(X_nei));             % 取选择的邻域
        neighborX(K1_X+2:end, :) = 0;                 % 默认只有K1个邻域，其与删除
            Sel_FlagX = false(size(Nei_X_sel));
            X01_ = Xt01(2:end, Ind_right);
            X02_ = Xt02(2:end, Ind_right);
            ThrDisX = Amp * DisYt(K1_1_Y, Ind_right);
            for i = 1 : N_right
                n = Ind_right(i);
                inv_Tyx = reshape(TransforMatXY_vec(n, :), 2,2);
                Xi = [X01_(:,i), X02_(:,i)];
                for j = 1 : Kmax_X
                    Dis = inv_Tyx * Xi(j, :)';
                    Dis = sum(Dis.^2);
                    Sel_FlagX(j, i) = Dis < ThrDisX(i);
                    if ~Sel_FlagX(j, i)
                        if Dis > 2*ThrDisX(i)
                            break;
                        end
                    end
                end
            end
        Sel_FlagX(1:num_min_X, :) = true;
            for i = 1:N_right
                temp_i = X_nei(Sel_FlagX(:, i), i);
                N_i = numel(temp_i);
                Nei_X_sel(1:N_i, i) = temp_i;

                flags = [true; Sel_FlagX(:, i)];
                n = Ind_right(i);
                neighborX_ori(:, n) = [neighborX_ori(flags, n); neighborX_ori(~flags, n)];
                DisXt(:, n)         = [DisXt(flags, n);         DisXt(~flags, n)];
                Xt01(:, n)          = [Xt01(flags, n);          Xt01(~flags, n)];
                Xt02(:, n)          = [Xt02(flags, n);          Xt02(~flags, n)];
            end
        neighborX(2:end, Ind_right) = Nei_X_sel;
        
        % ori
        neighborS_ori = neighborX_ori;
        if Kmax_1_X>Kmax_1_Y
            neighborS_ori(1:Kmax_1_Y, Ind_right) = neighborY_ori(1:Kmax_1_Y, Ind_right);
        else
            neighborS_ori(:, Ind_right) = neighborY_ori(1:Kmax_1_X, Ind_right);
        end
        % 以小的一侧为邻域   
%         Nei_sel_S = neighborX;
%         Nei_sel_L = neighborY;

        Nei_sel_S = zeros(max(K1_1_X, K1_1_Y), N);
        Nei_sel_L = zeros(max(Kmax_1_Y, Kmax_1_Y), N);
        Nei_sel_S(1:K1_1_X, indX) = neighborX(1:K1_1_X, indX);
        Nei_sel_S(1:K1_1_Y, indY) = neighborY(1:K1_1_Y, indY);          
        Nei_sel_L(1:Kmax_1_Y, indX) = neighborY(:, indX);
        Nei_sel_L(1:Kmax_1_X, indY) = neighborX(:, indY);
end


