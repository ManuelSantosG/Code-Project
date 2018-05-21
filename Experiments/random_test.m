clear
rrr=[];
C_u=[];
C_v=[];
delta=[];
for i=1:30;
    [r1,r2,r3,r4,C_uobs1,C_uobs2,C_uobs3,C_uobs4,C_vobs1,C_vobs2,C_vobs3,C_vobs4,delta_obs1,delta_obs2,delta_obs3,delta_obs4]=Bin_exp();
    rrr=[rrr, [r1,r2,r3,r4]'];
    C_u=[C_u, [C_uobs1,C_uobs2,C_uobs3,C_uobs4]'];
    C_v=[C_v, [C_vobs1,C_vobs2,C_vobs3,C_vobs4]'];
    delta=[delta, [delta_obs1,delta_obs2,delta_obs3,delta_obs4]'];
end
