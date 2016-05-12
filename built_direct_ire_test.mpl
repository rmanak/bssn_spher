#################################
# Using direct einstein's tensor
# equation to built IREs
################################
read("eq_ire.mpl"):

res2:=grcomponent(HSR):
tg2:=grcomponent(TG2(up),[x]):

with(DEtools):
addcoords(compact,[t,x],[t,ctfm(x)]);

ire_fields := [res2, tg2, restt, restx, resxx, resthth, ire_rpsi_direct]:
c_ire_fields := [0, 0, 0 , 0 , 0 , 0 , 0]:

# Compactifying Einstein's equation:
for ii from 1 to nops(ire_fields) do
 c_ire_fields[ii] := PDEchangecoords(ire_fields[ii],[t,x],compact,[t,x]);
 c_ire_fields[ii] := subs({diff(ctfm(x),x)=ctfmp(x),diff(ctfm(x),x,x)=ctfmpp(x)},c_ire_fields[ii]); 

 for kk from 1 to nops(grid_functions) do
        fn := grid_functions[kk];
        c_ire_fields[ii] := subs(fn(t,ctfm(x))=fn(t,x),c_ire_fields[ii]);
 end do:

end do:

printf("checking compactification... all residuals should be zero:\n");
for ii from 1 to nops(ire_fields) do
   res := simplify(ire_fields[ii] - eval(c_ire_fields[ii],{ctfm(x) = x, ctfmp(x) = 1, ctfmpp(x) = 0})):
   printf("res = %a\n",res);
end do:


for ii from 1 to nops(ire_fields) do
   c_ire_fields[ii] := collect(collect(eval(expand(c_ire_fields[ii]),{ctfmpp(x)=octfmp(x)*ctfmp(x)^2}),1/ctfmp(x)),1/ctfm(x));
end do:

res2 := c_ire_fields[1]:
tg2 := c_ire_fields[2]:
restt := c_ire_fields[3]:
restx := c_ire_fields[4]:
resxx := c_ire_fields[5]:
resthth := c_ire_fields[6]:
ire_rpsi_direct := c_ire_fields[7]:


Gen_Eval_Code(restx,input="c",proc_name="ire_restx");
Gen_Eval_Code(resxx,input="c",proc_name="ire_resxx");
Gen_Eval_Code(restt,input="c",proc_name="ire_restt");
Gen_Eval_Code(resthth,input="c",proc_name="ire_resthth");
Gen_Eval_Code(tg2-Lamx(t,x),input="c",proc_name="ire_val_lamx");
Gen_Eval_Code(ire_rpsi_direct,input="c",proc_name="ire_rpsi_direct");

pl:=table([ t=[-1,-1],x=[-1,0],y=[-1,-1],z=[-1,-1] ]):
Update_FD_Table(2,pl):

Gen_Eval_Code(res2,input="c",proc_name="ire_hs");
