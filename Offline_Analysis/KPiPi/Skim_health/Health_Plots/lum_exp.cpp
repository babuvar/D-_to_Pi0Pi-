void lum_exp()
{
int count=-1, a[1000],atemp=7,bcumul,A[50],COUNT=0;
float b[1000],c,B[50]={0};
ifstream fin;
//fin.open("info_onres.txt");
//fin.open("info_5s.txt");
fin.open("info_continuum.txt");
while(!fin.eof()){
count++;
fin>>a[count]>>c>>c>>b[count];
}
fin.close();

for (int i=0;i<count;i++){
//cout<<a[i]<<"\t"<<b[i]<<endl;
if(a[i]==atemp){ A[COUNT]=a[i]; B[COUNT]=B[COUNT]+b[i]; }
else{atemp=a[i]; COUNT++;A[COUNT]=a[i]; B[COUNT]=B[COUNT]+b[i]; }


}

for (int i=0;i<=COUNT;i++){
cout<<A[i]<<"\t"<<B[i]<<endl;
}

}
