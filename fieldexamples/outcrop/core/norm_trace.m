function [ data_out ] = norm_trace( data_in )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
data_out = zeros(size(data_in));
[~,ng] = size(data_in);


for k = 1:ng
%     if abs(data_in(:,k))~=0
%         data_out(:,k) = data_in(:,k)./max(abs(data_in(:,k)));
%     else
%         data_out(:,k) = data_in(:,k);
%     end

        data_out(:,k) = data_in(:,k)./max(abs(data_in(:,k)));
 

    

end


end

