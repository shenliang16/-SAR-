function [idx0,  threshold] = Threshold_With_minNum(Smaller_Larger, data, Threshold0, minnum, flag)

if nargin<5
    flag = 0;
end

if Smaller_Larger
        N = length(data);
        [dist_in_sort,~]=sort(data);   midInd = round(minnum);
        threshold=max(dist_in_sort(midInd), Threshold0); 
        if flag
            idx0 = find(data <= threshold); 
        else
            idx0 = find(data < threshold); 
        end
%         if length(idx0)>midInd
%             idx0 = idx0(1:midInd);
%         end
else
        N = length(data);
        [dist_in_sort,~]=sort(data, "descend");   midInd = round(minnum);
        threshold=min(dist_in_sort(midInd), Threshold0); 
        if flag
            idx0 = find(data >= threshold); 
        else
            idx0 = find(data > threshold); 
        end
%         if length(idx0)>midInd
%             idx0 = idx0(1:midInd);
%         end
end




%% 性能IMW-10000：SIFT-SURF-KAZE
%% FM-10030
% function [idx0,  threshold] = Threshold_With_minNum(Smaller_Larger, data, Threshold0, minnum, flag)
% 
% if nargin<5
%     flag = 0;
% end
% 
% if Smaller_Larger
%         N = length(data);
%         [dist_in_sort,~]=sort(data);   midInd = round(min(N/2, max(N/10, minnum)));
%         threshold=max(dist_in_sort(midInd), Threshold0); 
%         if flag
%             idx0 = find(data <= threshold); 
%         else
%             idx0 = find(data < threshold); 
%         end
% %         if length(idx0)>midInd
% %             idx0 = idx0(1:midInd);
% %         end
% else
%         N = length(data);
%         [dist_in_sort,~]=sort(data, "descend");   midInd = round(min(N/2, max(N/10, minnum)));
%         threshold=min(dist_in_sort(midInd), Threshold0); 
%         if flag
%             idx0 = find(data >= threshold); 
%         else
%             idx0 = find(data > threshold); 
%         end
% %         if length(idx0)>midInd
% %             idx0 = idx0(1:midInd);
% %         end
% end