function [rand_i]=randi_probability(probability,number,range)
% random integer with given probability in given range
% Inputs:
% 
% probability: values for probability [p_from integer .... p_to integer];
% number: number of random values
% range: [integers]
% 
% output:
% rand_i: vector (numberx1) with integers in range   
%
% Example usage:
% [rand_i]=randi_probability([0.394 0.34 0.131 0.099 0.036],100,[1 5])


p=probability; 
n=number;

if sum(p)-1>1e-9 
    
    [rand_i]=-1;
    disp('sum p not 1!');
    
elseif length(range)~=length(p)
    
    [rand_i]=-1;
    disp('wrong range');
    
else

rand_i=zeros(n,1);

[ps,idx] = sort(p);
ps = cumsum(ps);
x = rand(n,1);

E = x;
E(x<=ps(1)) = idx(1);

for ii = 2:length(p)
E(x>=ps(ii-1) & x<=ps(ii) ) = idx(ii);
end

ranges=range;

for ii=1:length(p)    
    rand_i(E==ii)=ranges(1,ii);    
end

end

end