function R = ddt(Xhat,dt,ddt_order)
% Implentation of several difference approximations to d/dt
% see % see https://en.wikipedia.org/wiki/Finite_difference_coefficient
%
%
if nargin == 2
    ddt_order = 4;
end

R = zeros(size(Xhat));

switch ddt_order
    
    
    case 11   %backward differencing          
        R = ([Xhat(2:end,:)-Xhat(1:end-1,:); Xhat(end,:)-Xhat(end-1,:)])/dt;
        
    case 12  % forward differencing
        R = ([Xhat(2,:)-Xhat(1,:); Xhat(2:end,:)-Xhat(1:end-1,:)])/dt;
        
    case 2 % central differencing
        R(1,:) = (-3*Xhat(1,:) + 4*Xhat(2,:) - Xhat(3,:))/(2*dt);
        R(end,:) = (3*Xhat(end,:) - 4*Xhat(end-1,:) + Xhat(end-2,:))/(2*dt);
        R(2:end-1,:) = (Xhat(3:end,:) - Xhat(1:end-2,:))/(2*dt);
        
    case 3
        error('3rd order implementation is not done')
    
    case 41 % uses mostly central differencing (same as Renee's code); % 4th order 5-point stencil from https://en.wikipedia.org/wiki/Five-point_stencil
        K = length(Xhat);
        R(3:end-2,:) = (1/12*Xhat(1:end-4,:) - 2/3*Xhat(2:end-3,:) + 2/3*Xhat(4:end-1,:) - 1/12*Xhat(5:end,:))/dt;
        R(1:2,:) = (-25/12*Xhat(1:2,:) + 4*Xhat(2:3,:) - 3*Xhat(3:4,:) + 4/3*Xhat(4:5,:) - 1/4*Xhat(5:6,:))/dt;
        R(K-1:K,:) = (25/12*Xhat(K-1:K,:) - 4*Xhat(K-2:K-1,:) +3*Xhat(K-3:K-2,:) -4/3*Xhat(K-4:K-3,:) + 1/4*Xhat(K-5:K-4,:))/dt;
   
    case 42 % uses mostly backwards differencing
        R(5:end,:) = (25/12*Xhat(5:end,:) - 4*Xhat(4:end-1,:) + 3*Xhat(3:end-2,:) - 4/3*Xhat(2:end-3,:) + 1/4*Xhat(1:end-4,:))/dt;
        R(1:4,:) = (-25/12*Xhat(1:4,:) + 4*Xhat(2:5,:) - 3*Xhat(3:6,:) + 4/3*Xhat(4:7,:) - 1/4*Xhat(5:8,:))/dt;
    
    case 43 % uses mostly forwards differencing   
        R(end-3:end,:) = (25/12*Xhat(end-3:end,:) - 4*Xhat(end-4:end-1,:) + 3*Xhat(end-5:end-2,:) - 4/3*Xhat(end-6:end-3,:) + 1/4*Xhat(end-7:end-4,:))/dt;
        R(1:end-4,:) = (-25/12*Xhat(1:end-4,:) + 4*Xhat(2:end-3,:) - 3*Xhat(3:end-2,:) + 4/3*Xhat(4:end-1,:) - 1/4*Xhat(5:end,:))/dt;

end