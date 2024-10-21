int add_acp()
{
double a[100],ea[100], Den=0, Num=0, Avg, Sig;
//a[0]= -0.010884;
//a[1]= 0.0953276;
//a[2]= 0.0397003;

//Category
//a[0]= -0.0195158;
//a[1]= 0.0162854;

//Merged fit of 2 and 3 bins
//a[1]= 0.0784584;//*
//a[1]= 0.0663729;

//ea[0]= 0.0207101;
//ea[1]= 0.0571001;
//ea[2]= 0.0865508;

//Category
//ea[0]= 0.0234895;
//ea[1]= 0.0425554;
 
//Merged fit of 2 and 3 bins
//ea[1]=  0.0476622;//*
//ea[1]= 0.0510919;

a[0]= -0.00793263;
a[1]= -0.0012647;
ea[0]= 0.00518099;
ea[1]= 0.0111542;


for(int i=0;i<=1;i++){
Num=Num+(a[i]/(ea[i]*ea[i]));
Den=Den+(1/(ea[i]*ea[i]));
}
Avg=Num/Den;
Sig=1/(sqrt(Den));

cout<<"Combined A_raw = "<<Avg<<" +/- "<<Sig<<endl;

return 0;
}
