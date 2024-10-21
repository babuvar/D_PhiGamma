#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cstdlib>
#include<cstring>
using namespace std;
int main()
{
char c[10];
string s,s1,s2;
s2=".txt";
int i,j=0,k,A[1000],B[1000],C[1000]={0},A1[100],B1[100];
ifstream fin;
stringstream ss;



for(i=1;i<=104;i++){//for loop
ss.str("");
ss<<i;
s1=ss.str();
s=s1+s2;
strcpy(c,s.c_str());
fin.open(c);
while(!fin.eof())
{
j++;
fin>>A[j];
fin>>B[j];
//cout<<A[j]<<"\t"<<B[j]<<endl;
if(fin.eof()) {j--;C[j]=1; break;}
}
//cout<<"----------------------------------------------------------------------"<<endl;
fin.close();
}//for loop

//cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;

//printing the stuff
/*
i=0;
for(k=1;k<=j;k++)
{
cout<<A[k]<<"\t"<<B[k]<<endl;
if(C[k]==1){i++; cout<<"---------------------------------------------------------------------- "<<i<<endl;}
}

cout<<"**************************"<<endl;*/

int n=0,ii,st;
for(k=1;k<=j;k++){
st=0;
for(ii=1;ii<=n;ii++){
if(A[k]==A1[ii]){
st=1;
B1[ii]=B1[ii]+B[k];}
}

if(st==0){
n++;
A1[n]=A[k];
B1[n]=B[k];
}
}// for loop

ofstream fout;
fout.open("TotalPID.txt");

//printing A1 & B1
for(int ij=1;ij<=n;ij++){
fout<<A1[ij]<<"\t"<<B1[ij]<<endl;
}
fout.close();

return 0;
}
