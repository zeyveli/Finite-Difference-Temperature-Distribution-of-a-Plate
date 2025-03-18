%% QUESTION 1 FINITE DIFFERENCE METHOD FOR TEMPERATURE DISTRIBUTION

%properties for the question
side_lenght = 1; %m of the square
delta = 0.1; %m, make dure the grids on each direction are odd numbered to get the middle value
node = side_lenght/delta - 1; %number of nodes on each row
TR = 50; %temperature of the right side
TL = 25; %temperature of the left side
TU = 100; %temperature of the up side
TD = 0; %temperature of the down side

%% Making the coefficient matrix and the resulting matrix

%defining matrices beforehand for convenience
M = zeros(node^2,node^2);
S = zeros(node^2,1);

% corners
M(1,1) = -4; M(1,2) = 1; M(1,node+1) = 1; S(1,1) = -(TD + TL);
M(node,node) = -4; M(node,node-1) = 1; M(node,2*node) = 1; S(node,1) = -(TD + TR);
M((node-1)*node+1,(node-1)*node+1) = -4; M((node-1)*node+1,(node-2)*node+1) = 1; M((node-1)*node+1,(node-1)*node+2) = 1; S((node-1)*node+1,1) = -(TU + TL);
M(node^2,node^2) = -4; M(node^2,node^2-1) = 1; M(node^2,node*(node-1)) = 1; S(node^2,1) = -(TU + TR);

for i = 2:node-1
    % sides
    M(i,i) = -4; M(i,i-1) = 1; M(i,i+1) = 1; M(i,node+i) = 1;  S(i,1) = -TD;
    M((i-1)*node+1,(i-1)*node+1) = -4; M((i-1)*node+1,(i-2)*node+1) = 1; M((i-1)*node+1,i*node+1) = 1; M((i-1)*node+1,(i-1)*node+2) = 1;  S((i-1)*node+1,1) = -TL;
    M(i*node,i*node) = -4; M(i*node,(i-1)*node) = 1; M(i*node,(i+1)*node) = 1; M(i*node,i*node-1) = 1;  S(i*node,1) = -TR;
    M((node-1)*node+i,(node-1)*node+i) = -4; M((node-1)*node+i,(node-1)*node+i-1) = 1; M((node-1)*node+i,(node-1)*node+i+1) = 1; M((node-1)*node+i,(node-2)*node+i) = 1;  S((node-1)*node+i,1) = -TU;

    % interior points
    for j = 2:node-1
        M((i-1)*node+j,(i-1)*node+j) = -4; M((i-1)*node+j,(i-1)*node+j+1) = 1; M((i-1)*node+j,(i-1)*node+j-1) = 1; M((i-1)*node+j,(i-2)*node+j) = 1; M((i-1)*node+j,i*node+j) = 1;
    end
end

%% Solving M*T=S

% solving with partial pivoting
for i = 1:node^2
    % switching rows of the biggest value
    M_piv = abs(M);
    [m,I] = max(M_piv([i:node^2],i));
    P_temp = eye(node^2);
    P_temp([i I+i-1],:) = P_temp([i I+i-1],:);
    P_e_5(:,:,i) = P_temp;
        
    M = P_e_5(:,:,i)*M; %multiplying to find new M
    S = P_e_5(:,:,i)*S; %multiplying to find new S
    
    % row operations
    L_temp = eye(node^2);
    for j = 1:node^2
        if i < j
            L_temp(j,i) = -M(j,i)/M(i,i);
        end
    end
    L_e_5(:,:,i) = L_temp;
    
    M = L_e_5(:,:,i)*M; %multiplying to find new M
    S = L_e_5(:,:,i)*S; %multiplying to find new S
end

% back propogation M*T = S
T = zeros(node^2,1);
for i = node^2:-1:1
    T(i,1) = S(i,1);
    j = node^2;
    while j > i
        T(i,1) = T(i,1) + (-M(i,j))*T(j,1);  
         j = j - 1;
    end
    T(i,1) = T(i,1)/M(i,i);
end

%% Making a table

T_matrix = zeros(node,node);
for i = 1:node
    for j = 1:node
        T_matrix(i,j) = T((node-i)*node+j);
    end
end

% export the matrix T_matrix into an excel file
xlswrite('numerik5_1.xlsx',T_matrix)

%% Plot

% isoterms
iso = zeros(node,node);
for i = 1:node
    iso(i,:) = T_matrix((node-i)+1,:); %I had to make T_max upside down, otherwise it gave a plot upside down
end
contour(delta:1/(node+1):1-delta,delta:1/(node+1):1-delta,iso,node*2) 
colormap hot
colorbar
title('Temperature Distribution','100')
xlabel("0")
colororder({'k','k'})
yyaxis left
ylabel("25")
yyaxis right
ylabel("50")
ylim([delta 1-delta])