function cdf_r = cdf(T,P,topR)
%     clear;
%     topR = 3;
%     T = [1 0 0 1; 0 0 0 1; 0 0 1 1; 0 1 1 0; 0 0 0 0];
%     P = [8 3 1 6; 3 4 3 7; 1 2 5 1; 5 7 2 3; 2 1 1 4];
    [rank,index] = sort(P,'descend');
    numSingl = 0;
    tmp_not_occur = 1;
    for r = 1:topR
         for i = 1:size(T,2);
            ind = find(T(:,i));
            if isempty(ind)
                if (r == 1)
                    numSingl = numSingl + 1;
                end;
                rate(i,r)= 0;
            else
                num = length(ind);
                count = 0;
                    for j = 1:num
                            if(index(r,i) == ind(j))
                                count = count + 1;
                            end
                    end
                    rate(i,r) = count;
            end
         end
         p_cur_occur = sum(rate(:, r)) / (size(rate, 1) - numSingl);
         p_cur_not_occur = 1 - p_cur_occur;
         tmp_not_occur = tmp_not_occur * p_cur_not_occur;
         cdf_r(r) = 1 - tmp_not_occur;
    end
% save result
% disp(['numsingle=',num2str(numSingl)]);
% save cdf_r2.mat cdf_r;
% save rate2.mat rate;

%---------------------------------
%     
% 
%     rate = sum(rate) / (size(rate,1) - numSingl);

%     [rank,index] = sort(P,'descend');
%     for r = 1:topR  
%         numSingl = 0 ;
%          for i = 1:size(T,2);
%             ind = find(T(:,i));
%             if isempty(ind)
%                 numSingl = numSingl + 1;
%                 rate(i,r)= 1;
%             else
%                 num = length(ind);
%                 count = 0;
%                     for j = 1:num
%                         for k = 1:r
%                             if(index(k,i) == ind(j))
%                                 count = count + 1;
%                             end
%                         end
%                     end
%                     rate(i,r) = count / r;
%             end
%          end
%          if(r==1)
%              cdf_r(r) = sum(rate(:,r)) /  (size(rate,1) - numSingl);
%          else
%              cdf_r(r) = sum(rate(:,r)) / (size(rate,1) - numSingl) + cdf_r(r-1);
%           end
%     end

%     rate = sum(rate) / (size(rate,1) - numSingl);
%     
%   for i = 1 : topR
%       if(i==1)
%           cdf_r(i) = rate(i);
%       else
%           cdf_r(i) = rate(i) + rate(i-1);
%       end
%   end

    % rate = mean(rate);

%     x = 1:topR
%     plot(x,cdf_r(1:topR),'r--','linewidth',2);
%     xlabel('Number of genes looked at');
%     ylabel('P(hidden gene among genes looked at)');
%     grid on

    % x =1:4;
    % plot(x,rate(1:4),'c--','linewidth',2);
    % %plot(x,many_rankcount_katz(1:100),'c--',x,many_rankcount_katz_hetersim(1:100),'b:',x,many_rankcount_catapult(1:100),'g-.',x,many_rankcount_catapult_hetersim(1:100),'r-',HumNet_x,HumNet_y,'m+',Pro_x,Pro_y,'k.','linewidth',2);
    % xlabel('Number of genes looked at');
    % ylabel('P(hidden gene among genes looked at)');
    % %legend('katz','HSMP','Catapult','HSSVM','PRINCE','ProDiGe');
    % grid on
                
        
        