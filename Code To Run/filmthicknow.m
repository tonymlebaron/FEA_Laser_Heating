function [thicknow] = filmthicknow(filmthick)
    global pcp
    if length(pcp(:,1))==1
        thicknow=filmthick;
    else
        thicknow=(min(pcp(pcp(:,1)>0,1)));
    end
end

