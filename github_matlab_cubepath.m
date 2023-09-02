clc
for i= 78 %%1:9879;  
    B1 = csvread('/Users/resh/Desktop/RA/CSVfiles/B5_nt_B5_t0.csv',74,0); %%%% CSV data input for first 5000 muons
    B2 = csvread('/Users/resh/Desktop/RA/CSVfiles/B5_nt_B5_t1.csv',74,0); %%%% CSV data input for next 5000 muons
    %size(B1)
    %size(B2)
    B = vertcat(B1,B2); %% Combine both CSV files to get data for 10,000 muons
    %size(B);
    XYZ= B(i,:); 
    i2 = find(XYZ, 1, 'last'); % number of columns
    XYZ=XYZ(5:i2); 
    XYZ=reshape(XYZ',3,[]);
    XYZ=XYZ'; 
    Length= 200; %% Length of cube in mm
    XYZ3=XYZ(:,3)+abs(XYZ(2,3)); %%translating z direction
    XYZ2=XYZ(:,2); %% y-coordinates
    XYZ1=XYZ(:,1); %% x-coordinates
    
    % YZ plane
    exit_disp_y=XYZ2(end); %% y co-ordinate of Detector 4 (D4)
    exit_angle_y=atan((XYZ2(end)-XYZ2(end-1))/(XYZ3(end)-XYZ3(end-1))); %% tan-1[(y of D4 - y of D3)/ (z of D4 - z of D3)] %%exit angle to D3 in y dir
    
    % XZ plane
    exit_disp_x=XYZ1(end); %% x co-ordinate of D4
    exit_angle_x=atan((XYZ1(end)-XYZ1(end-1))/(XYZ3(end)-XYZ3(end-1)));  %%exit angle to D3 in x dir
    
    %%%%%%%%%%%% GMTE algorithm %%%%%%%%%%%%%%%%%%%
    clearvars -except In* Length exit_angle_y exit_disp_y XYZ*;
    const1=13.6^2;beta=1;p_muon=3000;X0=3.2;DE=3; %% X0 (radiation length) is in mm X0_U = 3.2; X0_Al = 88.97; X0_Fe = 17.57; X0_Pb = 5.612; X0_W = 3.504
    xeta1=Length;
    theta=exit_angle_y;
    yeta=XYZ2(end-1);%exit_disp_y
    x1=0;x1_step=1;dt_step=floor((1/x1_step)*xeta1);
    for j=1:1:dt_step
        x1=x1+x1_step;x11(j)=x1;
        yprime=yeta-(xeta1-x1)*sin(theta);
        xeta=x1;
        term=(p_muon^2)-(DE*p_muon*xeta);
        variance_theta=const1*(xeta/X0)*(1/term)*((1+0.038*log(xeta/X0))^2);
        variance_y=(xeta^2)*variance_theta/3;
        variance_th_y=0.5*sqrt(3)*sqrt(variance_theta*variance_y);
        
        % S1 is covariance matrix
        S1=[variance_y variance_th_y;variance_th_y variance_theta];
        R0=[1 x1;0 1];Y0=[0;0];
        xeta=xeta1-x1;
        term=(p_muon^2)-(DE*p_muon*xeta);
        variance_theta=const1*(xeta/X0)*(1/term)*((1+0.038*log(xeta/X0))^2);
        variance_y=(xeta^2)*variance_theta/3;
        variance_th_y=0.5*sqrt(3)*sqrt(variance_theta*variance_y);
        % S2 is covariance matrix
        S2=[variance_y variance_th_y;variance_th_y variance_theta];
        R1=[1 xeta;0 1];Y2=[yeta;theta];
        
        %% y_GMTE
        Y_MLP=inv(inv(S1)+R1'*inv(S2)*R1)*(inv(S1)*R0*Y0+R1'*inv(S2)*Y2);
        y_mlp(:,j)=Y_MLP;
        Error=2.*inv(inv(S1)+R1'*inv(S2)*R1);Err(j)=Error(1,1);    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%% SLP algorithm %%%%%%%%%%%%%%%%%%%
    y_linear=x11.*(yeta/xeta1);
    
    %%%%%%%%%%%% POCA algorithm %%%%%%%%%%%%%%%%%%%
    veta=yeta-tan(theta)*xeta1;
    x_poca=-veta/(tan(theta));
    
    for wq=1:numel(x11)
    %     if x_poca<0
    %         continue
    %     end
        if x11(wq)<x_poca
            y_poca(wq)=0;
        else
            y_poca(wq)=tan(theta)*x11(wq)+veta;
        end
    end
    
    %%%%% 1sigma, 2sigma and 3sigma envelopes %%%%%
    y_mlp_plus1sigma=y_mlp(1,:)'+(Err(:).^0.5);
    y_mlp_minus1sigma=y_mlp(1,:)'-(Err(:).^0.5);
    y_mlp_plus2sigma=y_mlp(1,:)'+2.*(Err(:).^0.5);
    y_mlp_minus2sigma=y_mlp(1,:)'-2.*(Err(:).^0.5);
    y_mlp_plus3sigma=y_mlp(1,:)'+3.*(Err(:).^0.5);
    y_mlp_minus3sigma=y_mlp(1,:)'-3.*(Err(:).^0.5);
    
    %%%%%%% Trajectory of muon taken from GEANT4 %%%%%%%
    y_GEANT4 = XYZ2(1:end-1); 
    xnew = [1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,199]; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% 1st figure for individual muon path %%%%%%%%%%%%
    figure(1)
    hold on
    nexttile
    plot(xnew,y_GEANT4,'Color','k',LineWidth=4);
    hold on
    plot(x11,y_mlp(1,:),'-.','Color','r', LineWidth=2)
    hold on
    plot(x11,y_linear,'-.','Color','b',LineWidth=3);
    hold on
    plot(x11,y_poca,'--','Color','m',LineWidth=3);
    xlabel('X (mm)',"FontSize",15, "FontWeight","bold") ;
    ylabel('Y (mm)',"FontSize",15, "FontWeight","bold") ;
    legend('y-Muon','y-GMTE','y-SLP','y-PoCA','Location','northwest',"FontSize",15, "FontWeight","bold");
    ax = gca; 
    ax.FontSize = 13;
    ax.FontWeight = 'bold';
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%% 2nd figure for sigma cuts %%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    hold on
    nexttile
    plot(xnew,y_GEANT4,'Color','k',LineWidth=3);
    hold on
    plot(x11,y_mlp(1,:),'--','Color','r', LineWidth=3)
    hold on
    plot(x11,y_mlp_plus2sigma,'--',x11,y_mlp_minus2sigma,'--','Color','g',LineWidth=3);
    hold on
    plot(x11,y_mlp_plus3sigma,'-.',x11, y_mlp_minus3sigma,'-.', 'Color','b',LineWidth=3);
    hold on
    legend('y-Muon','y-GMTE', '2 \sigma','','3 \sigma','','y-PoCA','Location','northwest',"FontSize",15, "FontWeight","bold")
    xlabel('X (mm)',"FontSize",15, "FontWeight","bold") 
    ylabel('Y (mm)',"FontSize",15, "FontWeight","bold") 
    str = {'Muon Energy: 3GeV','Uranium cube: 20 cm side'};
    text(2,13,str, "FontSize",15, "FontWeight","bold");
    ax = gca; 
    ax.FontSize = 13;
    ax.FontWeight = 'bold';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




