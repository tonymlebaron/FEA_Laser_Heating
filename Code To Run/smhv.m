function [ numout ] = smhv(inputarg,varargin)
%Smooth continuous approximation of heaviside function
k=0.5; %Default k of no specified
if length(varargin)>0
k(1)=varargin{1};
end
numout=(1+exp(-2*k*inputarg)).^(-1);


end

