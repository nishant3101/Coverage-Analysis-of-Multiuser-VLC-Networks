function g = gain(thetaR, fov,n)
    if(thetaR >= 0 && thetaR <= fov)
        g = (n*n)/(sin(fov)*sin(fov));
    else
        g = 0;
    end
    
end