function K = calcK1D(x,u,T2,ks)
%Returns pcp, an nx2 array of all material points that have vaporized
%Format [x1,y1;x2,y2,...xn,yn]
%pcp stands for Phase Change Points
% disp(['time=',num2str(time)])
% x
% y
% u
global pcp %use global variable here to pass back and forth
%K=ks+1e5.*smhv(u-T2);

for i=1:length(x)
    
     K(i)=ks+1e5.*smhv(u(i)-T2);
     if ismember(x(i),pcp)==1 %If pcp already contains (x,y), then material is already vaporized
         %if x(i) ~= 0
            K(i)=ks+1e5; %So we keep K high
            %disp(['keepin it up, x=',num2str(x(i)),'y=',num2str(y(i))])
         %end
%      elseif u(i) >= (T2)-5
%          K(i)=ks+1e5;
%          %disp('Converted!')
%          
%          pcp=vertcat(pcp,x(i));
     end

end


end

