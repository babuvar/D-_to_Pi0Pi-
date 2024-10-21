int add_acp()
{
double a[100],ea[100], Den=0, Num=0, Avg, Sig;

/*
a[0]=  -0.00789045; 
ea[0]= 0.0280108;
a[1]= 0.00265593; 
ea[1]= 0.0233271;
a[2]= 0.0280664; 
ea[2]= 0.0209975;
a[3]= 0.0299439; 
ea[3]= 0.0169041;
*/

a[0]= -0.00789045;
ea[0]= 0.0280108;
a[1]= 0.0210928;
ea[1]= 0.0116103;



for(int i=0;i<=1;i++){
Num=Num+(a[i]/(ea[i]*ea[i]));
Den=Den+(1/(ea[i]*ea[i]));
}
Avg=Num/Den;
Sig=1/(sqrt(Den));

cout<<"Combined A_raw = "<<Avg<<" +/- "<<Sig<<endl;

return 0;
}
