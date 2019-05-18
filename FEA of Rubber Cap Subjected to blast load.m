
tau=0.2;
rho=0.001;
K0=500;
mu0=5;
lambda0=K0-(2/3)*mu0;
A=2;
b1=0.06;
b2=1.2;
sound_wave=sqrt((K0+(4*mu0)/3)/rho);
for study=1:3
x_elements=25*study;
y_elements=4*study;
nsteps=3500*study;
t=linspace(0,1,nsteps);
dt=t(2)-t(1);

mesh=uniform_mesh(x_elements,y_elements,50,8);
a=zeros(2*length(mesh.X),1);
v=zeros(2*length(mesh.X),1);
u=zeros(2*length(mesh.X),1);

ne=x_elements;
edge_conn=[1:ne; 2:ne+1];
work.dissipated=0;
work.external=0;
%Mass Matrix                                    
M = zeros(mesh.num_nodes,mesh.num_nodes);
for c = mesh.connectivity
        Me = zeros(4,4);
        xe = mesh.X(:,c);
        for q  = mesh.quad_points
            [N,dNdp] = shape(q);
            J=xe*dNdp;
            dNdx=dNdp/J;
            r=xe(1,:)*N;
            Me =Me+2*pi*r*N*N'*rho*mesh.quad_points(end)*det(J);
        end
        M(c,c) = M(c,c) + Me;
end


for b=1:nsteps
    
P=A*(exp(-t(b)/0.2))*(1-(t(b)/0.2));
mesh.x=[u(1:2:end)';u(2:2:end)']+mesh.X;
mesh.v=[v(1:2:end)';v(2:2:end)'];
mesh.u=[u(1:2:end)';u(2:2:end)'];
f=zeros(2*length(mesh.X),1);
fext=zeros(2*length(mesh.X),1);
fint=zeros(2*length(mesh.X),1);
work.internal=0;
wkin1=0;
wkin2=0;
work.kinetic=0;
work.internal=0;
for i=edge_conn
    ie=mesh.x(:,i);
    iv=mesh.v(:,i);
    dx=ie(:,2)-ie(:,1);
    le=ie(1,2)-ie(1,1);
    n=[dx(2), -dx(1)];
    traction=-P*n;
    N1=[0.5;0.5];
    r1=ie(1,:)*N1;
    fx=N1*traction*2*pi*r1;
    fext(2*i-1)=fext(2*i-1)+fx(:,1);
    fext(2*i)=fext(2*i)+fx(:,2);
    work.external=work.external+traction*iv*N1*2*pi*r1*le*dt;
end
plt.ext(b)=work.external;
%Stress Calculations
for g=mesh.connectivity
    Xe=mesh.X(:,g);
    xe=mesh.x(:,g);
    ve=mesh.v(:,g);
    ue=mesh.u(:,g);
    feint=zeros(8,1);
    h=Xe(1,2)-Xe(1,1);
    
    
    for q=mesh.quad_points
        [N,dNdp]=shape(q);
        r=xe(1,:)*N;
        R=Xe(1,:)*N;
        J0=Xe*dNdp;
        J=xe*dNdp;
        dNdx=dNdp/J;
        dNdX=dNdp/J0;
        F=xe*dNdX;
        F(3,3)=r/R;
        
        T=lambda0*log(det(F))*eye(3)+mu0*(F*F'-eye(3));
                sig = [T(1,1);T(2,2);T(3,3);
                       T(1,2);]/det(F);
                   
                B=[dNdx(:,1)',0,0,0,0;
                    0,0,0,0,dNdx(:,2)';
                    N'/r,0,0,0,0;
                    dNdx(:,2)',dNdx(:,1)'];
                
                L=([ve(1,:)'; ve(2,:)']'*B')';
                vol_strain_rate=sum(L(1:3));
                
                L_non_voigt=[L(1),L(4),L(4);L(4),L(2),L(4);L(4),L(4),L(3)];
                D=(L_non_voigt+L_non_voigt')/2;
                D_voigt=[D(1,1);D(2,2);D(3,3);D(1,3)];
                
%                 strain_voigt=([ue(1,:)'; ue(2,:)']'*B')';
%                 vol_strain=strain_voigt(1,1)+strain_voigt(2,1)+strain_voigt(3,1);
%                 dev_strain=strain_voigt-(1/3)*vol_strain*[1;1;1;0];
%                 
                if vol_strain_rate<0
                Q=b1*rho*sound_wave*abs(h)*abs(vol_strain_rate)+b2*rho*h*h*(vol_strain_rate)^2;
                
                else Q=0;
                end
                sig=sig-Q*[1;1;1;0];
                
                wt=2*pi*r*det(J)*q(end);
                
                feint =feint+B'*sig*wt;
                
                %Internal Energy
                work.dissipated=work.dissipated+(Q*[1,1,1,0])*D_voigt*wt;
                
work.internal=work.internal+0.5*lambda0*(log(det(F)))^2-mu0*(log(det(F)))+0.5*mu0*(trace(F'*F)-3)*wt;
                
    end
plt.int(b)=work.internal;
plt.dis(b)=work.dissipated;
    fint(2*g-1)=fint(2*g-1)+feint(1:4);
    fint(2*g)=fint(2*g)+feint(5:8);
end

f=fext-fint;

a(1:2:end)=M\f(1:2:end);
a(2:2:end)=M\f(2:2:end);

%Fixed Boundary condition

%Fixed-Edge Boundary Condition at right end
 right_side = find(mesh.X(1,:) == 50);
     [~, order] = sort(mesh.X(1, right_side));
    right_side = right_side(order);
    right_conn= [right_side(1:end-1); right_side(2:end)];
    
    a(2*right_side-1)=0;
    a(2*right_side)=0;
    
 %X-Values fixed at left end
left_side = find(mesh.X(1,:) == 0);
     [~, order] = sort(mesh.X(1, left_side));
    left_side = left_side(order);
    left_conn= [left_side(1:end-1); left_side(2:end)];
    
    a(2*left_side-1)=0;

if b==1
    alpha=0.5;
else
    alpha=1;
end
    v=v+alpha*dt*a;
    u=u+dt*v;
    umax(b)=max(u);
     
    %Plot
    p.vertices=mesh.x';
    p.faces=mesh.connectivity';
    p.facecolor='w';
    clf
    patch(p)
    xlim([0 50])
    ylim([-30 35])
    pause(0.0001);
    
    %Kinetic Energy
    wkin1=0.5*v(1:2:end)'*M*v(1:2:end);
    wkin2=0.5*v(2:2:end)'*M*v(2:2:end);
    work.kinetic=wkin1+wkin2;
    plt.kin(b)=work.kinetic;
    work.total=-1*(-work.internal-work.kinetic+work.dissipated+work.external);
    plt.total(b)=work.total;
    
end
% figure();
%     plot(t,plt.kin)
%     hold on
%     plot(t,plt.dis)
%     hold on
%     plot(t,plt.int)
%     hold on
%     plot(t,plt.ext)
%     hold on
%     plot(t,plt.total)
%     legend('Kinetic Energy', 'Dissipated Energy', 'Internal Energy', 'External Energy','Total Energy')
    
    figure();
    plot(t,umax)
    hold on
end   
% Creates a unifom mesh of ex by ey elements 
% with length of lx and ly in the x and y directions.
% The origin (lower left) can be optionally specified as x0, y0.
    function [m] = uniform_mesh(ex, ey, lx, ly, x0, y0)
    % if origin is not specified, then set it to zero.
    if nargin < 5, x0 = 0; end
    if nargin < 6, y0 = 0; end       
        
    m.num_nodes = (ex+1)*(ey+1);
    % Nodal reference coordinates.
    m.X    = zeros(2, m.num_nodes);   
    for j=1:ey+1
        for i=1:ex+1
            m.X(:,i+(j-1)*(ex+1)) = [(i-1)*lx/ex + x0; (j-1)*ly/ey + y0];
        end
    end
    m.num_elements = ex*ey;
    m.connectivity = zeros(4, m.num_elements);
    for j=1:ey
        for i=1:ex            
            % first node in element
            n0 = i+(j-1)*(ex+1);            
            m.connectivity(:,i+(j-1)*ex) = [n0; n0+1; n0+ex+2; n0+1+ex];            
        end
    end
        
    % 4 point Gaussian quadrature rule.
    m.quad_points = [-1, 1, 1,-1;
                     -1,-1, 1, 1] / sqrt(3);
    m.quad_points = [m.quad_points; 1,1,1,1];

    % Call as mesh.draw() to draw this mesh.
    m.draw =       @() plot_mesh(m);
    % Call as mesh.draw_nodal(f) to plot mesh colored by nodal value f.
    m.draw_nodal = @(f) plot_nodal(m,f);
end

function [] = plot_mesh(m)
    p.vertices = m.X';
    p.faces = m.connectivity';
    p.facecolor = 'none';
    patch(p);
end

function[] = plot_nodal(m, f)
    p.vertices = m.X';
    p.faces = m.connectivity';
    p.facecolor = 'interp';
    p.facevertexcdata = f;
    patch(p);    
end



function [N, dNdp] = shape(p)
    N = 0.25*[(1-p(1)).*(1-p(2));
              (1+p(1)).*(1-p(2));
              (1+p(1)).*(1+p(2));
              (1-p(1)).*(1+p(2))];
          
    dNdp = 0.25*[-(1-p(2)), -(1-p(1));
                  (1-p(2)), -(1+p(1));
                  (1+p(2)),  (1+p(1));
                 -(1+p(2)),  (1-p(1))];
end
