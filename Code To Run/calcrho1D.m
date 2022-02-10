function rho = calcrho1D(x,u,T2,T1,rhos,rhov)
rho=(rhos)*smhv(T2-u)+smhv(u - T1).*rhov;
global pcp
    for i=1:length(x)
        if ismember(x(i),pcp)==1 %If pcp already contains (x,y), then material is already vaporized
            rho(i)=rhov;
        else
            rho(i)=(rhos)*smhv(T2-u)+smhv(u - T1).*rhov;
        end
    end
end

