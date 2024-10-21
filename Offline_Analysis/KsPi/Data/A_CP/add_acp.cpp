int add_acp()
{
double a[100],ea[100], Den=0, Num=0, Avg, Sig;



a[0]= -0.00781403; 
ea[0]= 0.00517793;
a[1]=  -0.00343388;
ea[1]=   0.0104712;



for(int i=0;i<=1;i++){
Num=Num+(a[i]/(ea[i]*ea[i]));
Den=Den+(1/(ea[i]*ea[i]));
}
Avg=Num/Den;
Sig=1/(sqrt(Den));

cout<<"Combined A_raw = "<<Avg<<" +/- "<<Sig<<endl;

return 0;
}
