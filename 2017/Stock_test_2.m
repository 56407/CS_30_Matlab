%stock price test one

P_Ark = 2.64;
P_Omg  = 8.77;
T_Ark = 18.31;
testC = 0;

T_Kmd = 22.4;
P_Kmd = 2.02;
modX = 0;


%loss percentage of ark
for MOD = 1.15:+0.05:2
    %disp(MOD);
    modx = (MOD);
    testA = P_Ark * MOD;
    disp(testA + ("  Current percentage: ") + (MOD));
    
    
end
display(testA);

mod15 = 1.15;
mod20 = 1.20;



mod25 = 1.25;
mod30 = 1.30;
Lark = P_Ark * lark;
%disp(Lark);

C_Omg = P_Omg * mod15;
%disp(C_Omg);


display(Lark);

%fprintf("This is with the 15 percent modifier\n");

%current price of kmd, can be used in other variables

CP = TKmd * P_Kmd;
%disp(CP);

CPa = P_Ark * T_Ark;
%disp(CPa);

testB = 4.50 * 9.31;
disp(testB);
Total_port = CP + CPa +5.18 + 2.12 + 1.03;
fprintf("This is how much we currently have");


disp(Total_port);


