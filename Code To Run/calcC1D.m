function C = calcC1D(x,u,T2,T1,cps,cpv,delh)
C=(cps)*smhv(T2-u)+smhv(u - T1).*smhv(T2 - u).*delh+cpv.*smhv(u-T2);
global pcp
    for i=1:length(x)
        if ismember(x(i),pcp)==1 %If pcp already contains (x,y), then material is already vaporized
            C(i)=cpv;
        else
            C(i)=(cps)*smhv(T2-u)+smhv(u - T1).*smhv(T2 - u).*delh+cpv.*smhv(u-T2);
        end
    end
end

