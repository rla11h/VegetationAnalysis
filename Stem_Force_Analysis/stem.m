function thetadot = stem(t,theta,j,b,k)
    a = 5;
    g = -9.81;
    
    thetadot = [theta(2); -b/j*theta(2) - k/j*theta(1)];



end