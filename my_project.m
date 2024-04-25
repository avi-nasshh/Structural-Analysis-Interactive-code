%% SECTION TITLE
% DESCRIPTIVE TEXT
clc
clearvars
%% input data by the user
bay_info =cell(3,2);
prompt="Enter number of bays in the structure in the X direction";
bay_info(1,1)=inputdlg(prompt);
prompt="Enter number of bays in the structure in the Z direction";
bay_info(2,1)=inputdlg(prompt);
prompt="Enter the length of bays in the structure for X direction";
bay_info(1,2)=inputdlg(prompt);
prompt="Enter the length of bays in the structure for Z direction";
bay_info(2,2)=inputdlg(prompt);
prompt="Enter the number of storeys in the structure";
bay_info(3,1)=inputdlg(prompt);
prompt="Enter the height of storeys in the structure";
bay_info(3,2)=inputdlg(prompt);

%% Calculating total number of joints
dum_x=str2double(bay_info{1,1});
dum_y=str2double(bay_info{3,1});
dum_z=str2double(bay_info{2,1});
dum_lx=str2double(bay_info{1,2});
dum_ly=str2double(bay_info{3,2});
dum_lz=str2double(bay_info{2,2});
joint_num= (dum_x+1)*(dum_y+1)*(dum_z+1);
%% To save member co
joint_matrix =cell(joint_num,2);
n=1;
for j=0:dum_y
    for k=0:dum_z
        for i=0:dum_x
            dum_mat=zeros(1,3);
            dum_mat(1,1)=i*dum_lx;
            dum_mat(1,2)=-j*dum_ly;
            dum_mat(1,3)=k*dum_lz;
            joint_matrix{n,1}=dum_mat;
            n=n+1;
        end
    end
end
%% calculating total number of member
%%%% a=number of beams along x-direction
%%%% b=number of columns along y-direction
%%%% c=number of beams along z_direction
%%%% d=total number of beams
a=(dum_x*(dum_z+1)*dum_y);
b=dum_y*(dum_x+1)*(dum_z+1);
c=dum_z*(dum_x+1)*dum_y;
member_num=a+b+c;
%% ASSIGNING MEMBER NUMBER
member_info= cell(member_num,11); %  ["member nodes","type","E","G","Iz","Iy","A","L",J]
n=1;
for j=1:dum_y*2
    if rem(j,2)~=0
        for k=1:dum_z+1
            for i=1:dum_x
                t=(dum_x+1)*(k-1)+i+((j-1)/2)*(dum_z+1)*(dum_x+1);
                member_info{n,1}=[t t+1];
                member_info{n,2}="Bx";
                n=n+1;
            end
        end
        for k=1:dum_z
            for i=1:dum_x+1
                t=(dum_x+1)*(k-1)+i +((j-1)/2)*(dum_z+1)*(dum_x+1);
                member_info{n,1}=[t t+dum_x+1];
                member_info{n,2}="Bz";
                n=n+1;
            end
        end
    else
        for i=1:(dum_z+1)*(dum_x+1)
            t=i+0.5*(j-2)*(dum_z+1)*(dum_x+1);
            member_info{n,1}=[t t+(dum_z+1)*(dum_x+1)];
            member_info{n,2}="C";
            n=n+1;
        end
    end
end
%% PLOTTING OF THE STRUCTURE
for i=1:member_num
        x_1=joint_matrix{member_info{i}(1)}(1);
        x_2=joint_matrix{member_info{i}(2)}(1);
        y_1=joint_matrix{member_info{i}(1)}(2);
        y_2=joint_matrix{member_info{i}(2)}(2);
        z_1=joint_matrix{member_info{i}(1)}(3);
        z_2=joint_matrix{member_info{i}(2)}(3);
        xlabel('X axis')
        ylabel('Y axis')
        zlabel('Z axis')
        line([x_1 x_2],[y_1 y_2],[z_1 z_2])
        
        view(40,40)
        camup([0 1 0])
end
%% ShHOWING JOINT NUMBER
for i=1:joint_num
    cx=joint_matrix{i}(1,1);
    cy=joint_matrix{i}(1,2);
    cz=joint_matrix{i}(1,3);
    text(cx,cy,cz,num2str(i),'Color','red','FontSize',15)
    view(40,40)
end

%% SHOWING MEMBER NUMBER ON THE PLOT
for m=1:member_num
    k=member_info{m};
    dum_1=joint_matrix{k(1,1)};
    dum_2=joint_matrix{k(1,2)};
    mx=(dum_1(1,1)+dum_2(1,1))/2;
    my=(dum_1(1,2)+dum_2(1,2))/2;
    mz=(dum_1(1,3)+dum_2(1,3))/2;
    text(mx,my,mz,num2str(m),'Color','green','FontSize',15)
    movegui('northeast')
end

%% MATERIAL AND GEOMATRICAL PROPERTIES OF BEAM
figure()
rgbImage = imread("beam'.png");
imshow(rgbImage)
movegui("northwest")
prompt ={'Width(mm)','Depth(mm)','E(N/mm^2)','G(N/mm^2)'};
title ='Dimensions and properties of beam';
definput = {'500','500','200000','75000'};
numlines = [1 75;1 75;1 75;1 75];
options.Resize='on';
beam = inputdlg(prompt,title,numlines,definput,options);
bB = str2double(beam{1});
dB = str2double(beam{2});
EB = str2double(beam{3});
GB = str2double(beam{4});

A_B = bB*dB;
Iy_B = dB*((bB)^3)/12;
Iz_B = bB*((dB)^3)/12;
JB = Iy_B + Iz_B;
close(2)
%% MATERIAL AND GEOMATRICAL PROPERTIES OF COLUMN
figure(3)
rgbImage = imread("Column.png");
imshow(rgbImage)
movegui("northwest")
prompt = {'Width(mm)','Depth(mm)','E(N/mm^2)','G(N/mm^2)'};
title = 'Dimensions and Properties of Column';
numlines = [1 75;1 75; 1 75;1 75];
definput = {'300','400','200000','75000'};
options.Resize='on';
column = inputdlg(prompt,title,numlines,definput,options);
bC = str2double(column{1});
dC = str2double(column{2});
EC = str2double(column{3});
GC = str2double(column{4});

A_C = bC*dC;
Iy_C = dC*((bC)^3)/12;
Iz_C = bC*((dC)^3)/12;
JC = Iy_C + Iz_C;

close(3)
%% SAVING MEMBER INFO
for i=1:member_num
    if strcmp(member_info{i,2}(1),'Bx')
        member_info{i,3}=EB;
        member_info{i,4}=GB;
        member_info{i,5}=Iz_B;
        member_info{i,6}=Iy_B;
        member_info{i,7}=bB*dB;
        member_info{i,8}=dum_lx*1000;
        member_info{i,9}=bB*(dB^3)*((1/3)-0.21*(dB/bB)*(1-(dB^4/(12*bB^4))));
        member_info{i,10}=bB;
        member_info{i,11}=dB;
    elseif strcmp(member_info{i,2}(1),'Bz')
        member_info{i,3}=EB;
        member_info{i,4}=GB;
        member_info{i,5}=Iz_B;
        member_info{i,6}=Iy_B;
        member_info{i,7}=bB*dB;
        member_info{i,8}=dum_lz*1000;
        member_info{i,9}=bB*(dB^3)*((1/3)-0.21*(dB/bB)*(1-(dB^4/(12*bB^4))));
        member_info{i,10}=bB;
        member_info{i,11}=dB;
    else
        member_info{i,3}=EC;
        member_info{i,4}=GC; 
        member_info{i,5}=Iz_C;
        member_info{i,6}=Iy_C;
        member_info{i,7}=bC*dC;
        member_info{i,8}=dum_ly*1000;
        member_info{i,9}=dC*(bC^3)*((1/3)-0.21*(bC/dC)*(1-(bC^4/(12*dC^4))));
        member_info{i,10}=bC;
        member_info{i,11}=dC;
    end
end
%% SAVING STRUCTURAL DEGREE OF FREEDOM
for i=1:joint_num
    joint_matrix{i,2}=[6*i-5 6*i-4 6*i-3 6*i-2 6*i-1 6*i];
end

%% SAVING MEMBER PROPERTIES 
member_matrix=cell(member_num,7);
for i=1:member_num
    [E_t,J_t,Iz_t,Iy_t,A_t, L_t,G_t]=member_info{i,3:9};
    member_matrix{i,1}=find_local_k(E_t,J_t,Iz_t,Iy_t,A_t, L_t,G_t);
    t=member_info{i,2};
    k=find_Rm(t);
    member_matrix{i,2}=find_tm(k);
    member_matrix{i,3}=member_matrix{i,2}'*member_matrix{i,1}*member_matrix{i,2}; % Find stifness matrix of members in global level
    t_=member_info{i,1};
    member_matrix{i,4}=[joint_matrix{t_(1),2} joint_matrix{t_(2),2}];
end
num_unrest=(dum_z+1)*(dum_x+1)*dum_y;
%% joint info matrix how many members are actually conneted to this joints
for i=1:1:joint_num
    add=0;
    for j=1:1:member_num
        pp=member_info{j,1};
        if (pp(1)==i) || (pp(2)==i)
            add=add+1;
            dummy(add)=j;
        end
    end
    joint_matrix{i,3}=[add,dummy];
    dummy=zeros;
end
%% half band width
flag1=0;
for i=1:1:num_unrest
    dum_var=joint_matrix{i,3};
    dum_str=joint_matrix{i,2};
    dum_size=size(dum_str);
    min=dum_str(1);

    for j=1:1:dum_size(1,2)
        if dum_str(j)<min
            min=dum_str(j);
        end
    end
    flag=0;
    for k=1:1:dum_var(1)
        dum_var1=member_info{dum_var(1,k+1)};
        for ii=1:1:2
            if dum_var1(1,ii)~=i
                dum_str1=joint_matrix{dum_var1(1,ii),2};
                if flag==0
                    max=dum_str1(1);
                    flag=1;
                end
                dum_size1=size(dum_str1);

                for jj=1:1:dum_size1(1,2)
                    if dum_str1(1,jj)>max && dum_str1(1,jj)<=num_unrest*6 
                        max=dum_str1(1,jj);
                    end
                end
            end
        end
    end
    dum_hbw=max-min+1;
    if flag1==0
        B_W=dum_hbw;
        flag1=1;
    elseif dum_hbw>B_W
        B_W=dum_hbw;
    end
end

%% BAND WIDTH OF Kpp AND ASSIGNING BANDED MATRIX

%B_W=num_unrest*(6)-(num_unrest-(dum_x+1)*(dum_z+1));
K_banded_Global=zeros(num_unrest*6,B_W);
Kxp=zeros((joint_num-num_unrest)*6,num_unrest*6);
Kxx=zeros((joint_num-num_unrest)*6,(joint_num-num_unrest)*6);
for k = 1:1:member_num
    L_aMatrix=member_matrix{k,4};
    for i=1:1:12
        for j= 1:1:12
                ll=L_aMatrix(1,i);
                rr=L_aMatrix(1,j);
                if (ll<=num_unrest*6)&&(rr<=num_unrest*6)
                    if (ll<=rr) && (rr<B_W+ll)
                        dum_mat=R_to_kg(ll,rr);      %dum_mat(1,1),dum_mat(1,2)
                        K_banded_Global(dum_mat(1,1),dum_mat(1,2)) = K_banded_Global(dum_mat(1,1),dum_mat(1,2)) +member_matrix{k,3}(i,j);
                    end
                elseif (ll>num_unrest*6)&&(rr<=num_unrest*6)
                    Kxp(ll-6*num_unrest,rr)=Kxp(ll-6*num_unrest,rr)+ member_matrix{k,3}(i,j);
                elseif (ll>num_unrest*6)&&(rr>num_unrest*6)
                    Kxx(ll-6*num_unrest,rr-6*num_unrest)=Kxx(ll-6*num_unrest,rr-6*num_unrest)+ member_matrix{k,3}(i,j);
                else
                end
        end
    end
end

%% FORCE INPUT
P_force_matrix=zeros(num_unrest*6,1);
prompt="Enter number of forces you wanted to add on the Structure";
title = 'FORCE INPUTS';
numlines = [1 70];
bays_no_X=inputdlg(prompt,title,numlines);
x=str2double(bays_no_X);
figure()
plot3(0,0,0)
a=[-1 1 -1 1 -1 1];
axis(a)
text(-0.06,0,0.2,'.','Color','black','FontSize',40)

x1 = [0 0.4];
y1 = [0 0];
z1=[0 0];
line(x1,y1,z1,'Color','black','LineStyle','-')
x1 = [0 0];
y1 = [0 0.4];
z1=[0 0];
line(x1,y1,z1,'Color','black','LineStyle','-')
x1 = [0 0];
y1 = [0 0];
z1=[0 0.4];
line(x1,y1,z1,'Color','black','LineStyle','-')
text(0.45,0,0,'fx','Color','black','FontSize',15)
text(0,0.6,0,'fz','Color','black','FontSize',15)
text(0,0,0.55,'fy','Color','black','FontSize',15)
text(0.65,0,0,'mx','Color','red','FontSize',15)
text(0,0.9,0,'mz','Color','red','FontSize',15)
text(0,0,0.75,'my','Color','red','FontSize',15)
axis off
movegui("southwest")
for i=1:x
    prompt = {'Enter the node number where you want to add forces','Fx (kN)','Fy (kN)','Fz (kN)','Mx (kNm)','My (kNm)','Mz (kNm)'};
    title = 'FORCE INPUT';
    numlines = [1 75;1 75; 1 75;1 75;1 75;1 75;1 75];
    definput = {'0','0','0','0','0','0','0'};
    options.Resize='on';
    column = inputdlg(prompt,title,numlines,definput,options);
    j = str2double(column{1});
    k=joint_matrix{j,2};
    P_force_matrix(k(1,1)) = 1000*str2double(column{2});
    P_force_matrix(k(1,2)) = 1000*str2double(column{3});
    P_force_matrix(k(1,3)) = str2double(column{4});
    P_force_matrix(k(1,4)) = 10^6*str2double(column{5});
    P_force_matrix(k(1,5)) = 10^6*str2double(column{6});
    P_force_matrix(k(1,6)) = 10^6*str2double(column{7});
end

%% Calculation of all the displacement
u_unres=find_chlsky_band_inv(K_banded_Global,P_force_matrix);
F_rxn=0.001*Kxp*u_unres;
close(2)
%% To calculate member end forces
%member_matrix=cell(member_num,3)
T_disp=zeros(joint_num*6,1);
for i=1:num_unrest*6
    T_disp(i)=u_unres(i);
end
for i=1:member_num
    %mem_asso=zeros(1,12);
    mem_asso=member_matrix{i,4};
    for j=1:12
        member_matrix{i,5}(1,j)=T_disp(mem_asso(1,j));
    end
    member_matrix{i,6}=member_matrix{i,2}*member_matrix{i,5}';
    member_matrix{i,7}=member_matrix{i,1}*member_matrix{i,6};
end

%% writing outputs to text data
fid= fopen("Noutput_file.txt",'w+t');
if fid < 0
    fprintf('error opening file\n');
    return;
end
Title_1="OUTPUT FILE OF STRUCTURAL MODEL CONTAINING MEMBER PROPERTIES AND REACTIONS";
Title_2="MEMBER PROPERTIES";
Label_0= 'S.no';
Label_1= 'TYPE';
Label_2 = 'E';
Label_3 = 'G';
Label_4 = 'Iz';
Label_5 = 'Iy';
Label_6 = 'AREA';
Label_7 = 'LENGTH';
Label_8 = 'J';

fprintf(fid,"\n\n\t\t\t\t\t\t\t\t\t\t%s\n\n\n",Title_1);
fprintf(fid,"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n\n\n",Title_2);
fprintf (fid, '\n') ;
fprintf(fid,  '%4s\t\t%4s\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t\t\t%4s\t\t\t%4s\t\t\t%4s\t\t\t%4s\t\t\t\n \n'  , Label_0,Label_1, Label_2, Label_3,Label_4, Label_5,Label_6,Label_7,Label_8);
for ii = 1:member_num
    fprintf(fid,'%-5d\t:\t',ii);
    fprintf(fid, '%5s\t\t%05e\t\t%05e\t\t%05e\t\t%05e\t\t%05.2f\t\t%05.2f\t\t%05.2f\t\t \n',member_info{ii,2:9});
end
fprintf(fid, '\n' ) ;

Label_0= 'S.no';
Label_1= '1';
Label_2 = '2';
Label_3 = '3';
Label_4 = '4';
Label_5 = '5';
Label_6 = '6';
Label_7 = '7';
Label_8 = '8';
Label_9 = '9';
Label_10 = '10';
Label_11 = '11';
Label_12 = '12';
 fprintf (fid, '\n') ;
  fprintf (fid, '\n') ;

Title_3="MEMBER-END DISPLACEMENTS";
fprintf(fid,"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n\n\n",Title_3);
fprintf (fid, '\n') ;
fprintf(fid,  '%4s\t\t\t\t%4s\t\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t\n \n'  , Label_0,Label_1, Label_2, Label_3,Label_4, Label_5,Label_6,Label_7,Label_8,Label_9,Label_10,Label_11,Label_12);
for ii = 1:member_num
    fprintf(fid,'%-5d\t:\t',ii);
    for jj=1:12
    fprintf(fid, '\t%05e\t',member_matrix{ii,6}(jj));
    end
    fprintf (fid, '\n') ;
end
fprintf (fid, '\n') ;
fprintf (fid, '\n') ;
fprintf (fid, '\n') ;
Title_3="MEMBER-END FORCES";
fprintf(fid,"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n\n\n",Title_3);
fprintf (fid, '\n') ;
fprintf(fid,  '%4s\t\t\t%4s\t\t\t%4s\t\t\t%4s\t\t\t%4s\t\t\t\t\t%4s\t\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t\t%4s\t\t\t\n \n'  , Label_0,Label_1, Label_2, Label_3,Label_4, Label_5,Label_6,Label_7,Label_8,Label_9,Label_10,Label_11,Label_12);
for ii = 1:member_num
    fprintf(fid,'%-5d\t:\t',ii);
    %for jj=1:12
    %fprintf(fid, '\t%05e\t',member_matrix{ii,7}(jj));
    %end
    for jj=1:12
        if fix((jj-1)/3)==0 || fix((jj-1)/3)==2
            fprintf(fid, '%8f\t\t',0.001*member_matrix{ii,7}(jj));
            %yy=yy+1;
        else
            fprintf(fid, '%8f\t\t',0.000001*member_matrix{ii,7}(jj));
            %yy=yy+1;
        end
    end
    fprintf (fid, '\n') ;
end

fprintf (fid, '\n') ;
fprintf (fid, '\n') ;
fprintf (fid, '\n') ;
Title_3="NODAL DISPLACEMENTS";
fprintf(fid,"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n\n\n",Title_3);
fprintf (fid, '\n') ;
fprintf(fid,  '%4s\t\t\t%-4s\t\t%-4s\t\t\t%-4s\t\t\t%-4s\t\t%-4s\t\t%-4s\t\n\n'  , Label_0,Label_1, Label_2, Label_3,Label_4, Label_5,Label_6);
yy=1;
for ii = 1:joint_num
    fprintf(fid,'%-5d\t:\t',ii);
    for jj=1:6
    fprintf(fid, '%5.2f\t\t',T_disp(yy));
    yy=yy+1;
    end
    fprintf (fid, '\n') ;
end
Label_0= 'S.no';
Label_1= 'Fx';
Label_2 = 'Fy';
Label_3 = 'Fz';
Label_4 = 'Mx';
Label_5 = 'My';
Label_6 = 'Mz';
fprintf (fid, '\n') ;
fprintf (fid, '\n') ;
fprintf (fid, '\n') ;
Title_3="REACTIONS";
fprintf(fid,"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n\n\n",Title_3);
fprintf (fid, '\n') ;
fprintf(fid,  '%4s\t\t\t%-4s  \t\t%-4s\t\t\t%-4s\t\t\t%-4s\t\t%-4s\t\t%-4s\t\n\n'  , Label_0,Label_1, Label_2, Label_3,Label_4, Label_5,Label_6);
yy=1;
for ii = 1:(joint_num-num_unrest)
    fprintf(fid,'%-5d\t:\t',num_unrest+ii);
    for jj=1:6
        if fix((jj-1)/3)==0  
            fprintf(fid, '%-5.2f\t\t',F_rxn(yy));
            yy=yy+1;
        else
            fprintf(fid, '%-5.2f\t\t',0.001*F_rxn(yy));
            yy=yy+1;
        end
    end
    fprintf (fid, '\n') ;
end
fprintf(fid, '\n' ) ;
fclose(fid);
%% float Reaction table
for i=1:(joint_num-num_unrest)
    Fx_kN(i,1) = F_rxn(6*i-5,1);
end

for i=1:(joint_num-num_unrest)
    Fy_kN(i,1) = F_rxn(6*i-4,1);
end

for i=1:(joint_num-num_unrest)
    Fz_kN(i,1) = F_rxn(6*i-3,1);
end

for i=1:(joint_num-num_unrest)
    Mx_kNm(i,1) = F_rxn(6*i-2,1)*0.001;
end

for i=1:(joint_num-num_unrest)
    My_kNm(i,1) = F_rxn(6*i-1,1)*0.001;
end

for i=1:(joint_num-num_unrest)
    Mz_kNm(i,1) = F_rxn(6*i,1)*0.001;
end
for i= 1:(joint_num-num_unrest)
    tt(i)=num_unrest+i;
end
Node_Number=tt';



R_react=[Node_Number,Fx_kN,Fy_kN,Fz_kN,Mx_kNm,My_kNm,Mz_kNm];
Reactions = [table(Node_Number),table(Fx_kN),table(Fy_kN),table(Fz_kN),table(Mx_kNm),table(My_kNm),table(Mz_kNm)];
fig_1=uifigure('Position',[50 550 760 150]);
uit1 = uitable(fig_1,"Data",Reactions);
uit1.Position = [10 10 720 120];

%displacements

Ux_mm=zeros(joint_num,1);
Uy_mm=zeros(joint_num,1);
Uz_mm=zeros(joint_num,1);
Rotation_x_rad=zeros(joint_num,1);
Rotation_y_rad=zeros(joint_num,1);
Rotation_z_rad=zeros(joint_num,1);
Joint_No=zeros(joint_num,1);


for i=1:joint_num
    Ux_mm(i,1)=T_disp(6*i-5);
    Uy_mm(i,1)=T_disp(6*i-4);
    Uz_mm(i,1)=T_disp(6*i-3);
    Rotation_x_rad(i,1)=T_disp(6*i-2);
    Rotation_y_rad(i,1)=T_disp(6*i-1);
    Rotation_z_rad(i,1)=T_disp(6*i);
    Joint_No(i,1)=i;
end


nodal_displacements = [table(Joint_No),table(Ux_mm),table(Uy_mm),table(Uz_mm),table(Rotation_x_rad),table(Rotation_y_rad),table(Rotation_z_rad)];

fig_2 = uifigure('Position',[50 150 760 350]);
uit2 = uitable(fig_2,"Data",nodal_displacements);
uit2.Position = [10 10 720 320];

%% FUNCTION TO DEFINE LOCAL STIFFNESS MATRIX
function flk=find_local_k(E_t,J_t,Iz_t,Iy_t,A_t, L_t,G_t) 
            flk = [ E_t*A_t/L_t  0  0 0 0 0 -E_t*A_t/L_t 0 0 0 0 0;
            0 12*E_t*Iz_t/L_t^3 0 0 0 6*E_t*Iz_t/L_t^2 0 -12*E_t*Iz_t/L_t^3 0 0 0 6*E_t*Iz_t/L_t^2;
            0 0 12*E_t*Iy_t/L_t^3 0 -6*E_t*Iy_t/L_t^2 0 0 0 -12*E_t*Iy_t/L_t^3 0 -6*E_t*Iy_t/L_t^2 0; 
            0 0 0 J_t*G_t/L_t 0 0 0 0 0 -J_t*G_t/L_t 0 0;
            0 0 -6*E_t*Iy_t/L_t^2 0 4*E_t*Iy_t/L_t 0 0 0 6*E_t*Iy_t/L_t^2 0 2*E_t*Iy_t/L_t 0 ;
            0 6*E_t*Iz_t/L_t^2 0 0 0 4*E_t*Iz_t/L_t 0 -6*E_t*Iz_t/L_t^2 0 0 0 2*E_t*Iz_t/L_t;
            -E_t*A_t/L_t 0 0 0 0 0 E_t*A_t/L_t 0 0 0 0 0;
            0 -12*E_t*Iz_t/L_t^3 0 0 0 -6*E_t*Iz_t/L_t^2 0 12*E_t*Iz_t/L_t^3 0 0 0 -6*E_t*Iz_t/L_t^2;
            0 0 -12*E_t*Iy_t/L_t^3 0 6*E_t*Iy_t/L_t^2 0 0 0 12*E_t*Iy_t/L_t^3 0 6*E_t*Iy_t/L_t^2 0;
            0 0 0 -J_t*G_t/L_t 0 0 0 0 0 J_t*G_t/L_t 0 0;
            0 0 -6*E_t*Iy_t/L_t^2 0 2*E_t*Iy_t/L_t 0 0 0 6*E_t*Iy_t/L_t^2 0 4*E_t*Iy_t/L_t 0 ;
            0 6*E_t*Iz_t/L_t^2 0 0 0 2*E_t*Iz_t/L_t 0 -6*E_t*Iz_t/L_t^2 0 0 0 4*E_t*Iz_t/L_t
          ];
end

%% Calculate R_transformaton matrix
function Rm=find_Rm(t)
    if strcmp(t,'Bx')
        Rm=[1 0 0;0 1 0;0 0 1];
    elseif strcmp(t,'Bz')
        Rm=[0 0 -1;0 1 0;1 0 0];
    else
        Rm=[0 -1 0;1 0 0;0 0 1];
    end
end
%% Function to get complete Transformation matrix

function tm=find_tm(k)
    Zero=zeros(3);
    tm=[k Zero Zero Zero;Zero k Zero Zero;Zero Zero k Zero;Zero Zero Zero k];
end


%% function to calculate inverse using cholesky
function inv=find_chlsky_band_inv(a,b)
K=size(a);
p=size(b);
w=zeros(p);
u=zeros(p);
upr_tnglr_mat=zeros(K(1),K(1));
for i=1:K(1)
    for j=1:K(1)
        if i==1 && j==1
                crs_dof=R_to_kg(i,j);
                if crs_dof(1,2)>=K(1,2)+crs_dof(1,1)
                    val=0;
                else
                    val=a(crs_dof(1,1),crs_dof(1,2));
                end
                upr_tnglr_mat(i,j)=val^0.5;
                aa=upr_tnglr_mat(1,1);
        elseif i>1 &&j<i
            crs_dof=R_to_kg(1,i);
                if crs_dof(1,2)>=K(1,2)+crs_dof(1,1)
                    val=0;
                else
                    val=a(crs_dof(1,1),crs_dof(1,2));
                end
                upr_tnglr_mat(1,i)=val/aa;
                sum=0;
                for k=1:i-1
                    sum=sum+upr_tnglr_mat(k,i)^2;
                end
                if j>=K(1,2)+i
                    val=0;
                else
                    crs_dof=R_to_kg(i,i);
                    val=a(crs_dof(1,1),crs_dof(1,2));
                end
                upr_tnglr_mat(i,i)= (val-sum)^0.5;
        elseif j>i
                sum=0;
                for k=1:i-1
                    sum=sum+upr_tnglr_mat(k,i)*upr_tnglr_mat(k,j);
                end
                if j>=K(1,2)+i
                    val=0;
                else
                    crs_dof=R_to_kg(i,j);
                    val=a(crs_dof(1,1),crs_dof(1,2));
                end
                upr_tnglr_mat(i,j)= (val-sum)/upr_tnglr_mat(i,i);
        elseif i>j
                upr_tnglr_mat(i,j)=0;
        end
    end
end
ttt=0;
for k=1:i
    ttt=ttt+k;
end
upr_tnglr_mat_t=upr_tnglr_mat';
V=zeros(ttt,1);
i=1;
for j= 1: K(1)
    for k=1:K(1)
        if k<=j
            V(i)=upr_tnglr_mat_t(j,k);
            i=i+1;
        end
    end
end
for i=1:p(1)
    tt=0;
    for k=1:i
        tt=tt+k;
    end
    sum=b(i)/V(tt);
    if rem(i,2)==1
        t_=fix(i/2)*i+1;
        for j=1:i-1
            sum=sum-(V(t_+(j-1))*w(j))/V(tt);
        end
        w(i)=sum;
    else
        t_=0.5*i^2-((i/2)-1);
        for j=1:i-1
            sum=sum-(V(t_+(j-1))*w(j))/V(tt);
        end
        w(i)=sum;
    end
end
VV=zeros(ttt,1);
i=1;
for j= 1: K(1)
    for k=1:K(1)
        if k>=j
            VV(i)=upr_tnglr_mat(j,k);
            i=i+1;
        end
    end
end
kk=ttt-1;
flag=0;
for i=p(1):-1:1
    tt=ttt;
    if i~=p(1)
        flag=1;
    end
    for k=1:p(1)-i
        tt=tt-(k+1);
    end
    sum=w(i)/VV(tt);
    for j=1:p(1)-i
        sum=sum-(VV(kk)*u(p(1)-(j-1)))/VV(tt);
        if i~=p(1)-j
            kk=kk-1;
        end
    end
    u(i)=sum;
    if flag~=0
        kk=kk-2;
    end
end
inv=u;
end

%% tranfer function from banded to global
function t=R_to_kg(i,j)
    if j>=i
        trns_j=j-i+1;
        trns_i=i;
        t=[trns_i,trns_j];
    else
        k=i;
        i=j;
        j=k;
        trns_j=j-i+1;
        trns_i=i;
        t=[trns_i,trns_j];
    end
end