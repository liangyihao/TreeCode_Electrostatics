/*
这段代码集中了几乎所有树操作的代码，外用函数
void init_tree()///external
void move_tree(int I,double dx,double dy,double dz)///移动离子并更新树
double Energy_tree(int I,tnode*p,double Dx,double Dy.double Dz)//计算以p->的内容为根的树中离子在第I个离子处的能量，返回能量量纲为1/[L]，注意换算
bool insert_particle(int I,tnode*p)///将第I个粒子的信息插入根节点为p的树
bool delete_particle(int I,tnode*p)///将第I个粒子的信息从根节点为p的树中删除
*/
#include "Tree_Def.h"
#include <cmath>
#include <iostream>
using namespace std;
int MAXLEVEL=4;///MAXLEVEL不要大于5
int MAXP=2;///表示最高三级展开
double BAC=0.5;
double err=1.0/1000;///err的量纲1/[L]
const int N0=2;///对于非叶结点,当所存粒子数小于等于N0时,激活这个结点的particle数组，存储这个结点含有的粒子
///树信息
tnode root;//root 的level为0，节点的level是从0~maxlevel
double tMax_x,tMax_y,tMax_z,tMin_x,tMin_y,tMin_z;///树计算所用的盒子范围(真实电荷及镜像电荷的范围)
int tposition[100000][10];///电荷在tree上面的索引
int tNpar;///树处理的离子数-1，即0~tNpar
int index_p[100000][10];///The index of particle in each node's of tree

int co[8][3]={{0,0,0},
              {0,0,1},
              {0,1,0},
              {0,1,1},
              {1,0,0},
              {1,0,1},
              {1,1,0},
              {1,1,1},
              };///在生成结点的最大最小坐标时使用


void expand(tnode* r)
{
for(int i=0;i<=MAXP;i++)for(int j=0;j<=MAXP;j++)for(int k=0;k<=MAXP;k++)r->moment[i][j][k]=0;

if(r->level<MAXLEVEL){
r->particle=new int[N0+1];
for(int i=0;i<=N0;i++)r->particle[i]=-1;
}
if(r->level==MAXLEVEL){
r->particle=new int[MAX_Npar_incell];
for(int i=0;i<=MAX_Npar_incell-1;i++)r->particle[i]=-1;
}
r->numpar=0;///先清理一下多极矩信息和离子信息
if(r->level>=MAXLEVEL)return;///如果到达最细的剖分那么不继续扩展
for(int i=0;i<=7;i++)
{
r->child[i]=new tnode;
if(co[i][0]==0)
  {
    r->child[i]->x_max=r->x_mid;
    r->child[i]->x_min=r->x_min;
  }else
  {
    r->child[i]->x_max=r->x_max;
    r->child[i]->x_min=r->x_mid;
  }

if(co[i][1]==0)
  {
    r->child[i]->y_max=r->y_mid;
    r->child[i]->y_min=r->y_min;
  }else
  {
    r->child[i]->y_max=r->y_max;
    r->child[i]->y_min=r->y_mid;
  }
if(co[i][2]==0)
  {
    r->child[i]->z_max=r->z_mid;
    r->child[i]->z_min=r->z_min;
  }else
  {
    r->child[i]->z_max=r->z_max;
    r->child[i]->z_min=r->z_mid;
  }
r->child[i]->x_mid=(r->child[i]->x_max+r->child[i]->x_min)/2;
r->child[i]->y_mid=(r->child[i]->y_max+r->child[i]->y_min)/2;
r->child[i]->z_mid=(r->child[i]->z_max+r->child[i]->z_min)/2;
r->child[i]->sqradius=((r->child[i]->x_max-r->child[i]->x_min)*(r->child[i]->x_max-r->child[i]->x_min)+(r->child[i]->y_max-r->child[i]->y_min)*(r->child[i]->y_max-r->child[i]->y_min)+(r->child[i]->z_max-r->child[i]->z_min)*(r->child[i]->z_max-r->child[i]->z_min))/4;
r->child[i]->radius=sqrt(r->child[i]->sqradius);
r->child[i]->level=r->level+1;

expand(r->child[i]);
}///creat children
///!////////////////////////
}

extern double x[70000],y[70000],z[70000],q[70000];
extern int tposition[100000][10];

bool insert_particle(int I,tnode*p)///将第I个粒子的信息插入根节点为p的树,!!!高重用率代码段
{
if((x[I]>=p->x_max)||(x[I]<p->x_min)||(y[I]>=p->y_max)||(y[I]<p->y_min)||(z[I]>=p->z_max)||(z[I]<p->z_min))return false;
p->numpar++;
if((p->level==MAXLEVEL)&&(p->numpar>MAX_Npar_incell))cout<<"wrong in box x:["<<p->x_min<<','<<p->x_max<<")  y:["<<p->y_min<<','<<p->y_max<<")  z:["<<p->z_min<<','<<p->z_max<<")"<<endl;
for(int i=0;i<=MAXP;i++)
 for(int j=0;j<=MAXP-i;j++)
  for(int k=0;k<=MAXP-i-j;k++)
   p->moment[i][j][k]+=q[I]*pow(x[I]-p->x_mid,i)*pow(y[I]-p->y_mid,j)*pow(z[I]-p->z_mid,k);

if(p->level==MAXLEVEL)
 {
  p->particle[p->numpar-1]=I;
  index_p[I][p->level]=p->numpar-1;
  return true;
 }
if(p->numpar<=N0)
 {
  p->particle[p->numpar-1]=I;
  index_p[I][p->level]=p->numpar-1;
 }
if(p->level<MAXLEVEL)
 for(int i=0;i<=7;i++)
 if(insert_particle(I,p->child[i]))
 {
 tposition[I][p->level]=i;
 return true;
 }
return false;
}

void delete_particle(int I,tnode*p)
{
tnode*p1=p;
tnode*p_temp[10];int temp=-1;
bool mark=true;
while(mark)
 {
  p1->numpar--;
  if((p1->numpar==N0)&&(p1->level!=MAXLEVEL)){temp++;p_temp[temp]=p1;}
  if((p1->numpar<N0)||(p1->level==MAXLEVEL))
   {
	   p1->particle[index_p[I][p1->level]]=p1->particle[p1->numpar];
	   index_p[p1->particle[p1->numpar]][p1->level]=index_p[I][p1->level];
   }
  for(int i=0;i<=MAXP;i++)for(int j=0;j<=MAXP-i;j++)for(int k=0;k<=MAXP-i-j;k++)p1->moment[i][j][k]-=q[I]*pow(x[I]-p1->x_mid,i)*pow(y[I]-p1->y_mid,j)*pow(z[I]-p1->z_mid,k);
  if(p1->level==MAXLEVEL)mark=false;
  p1=p1->child[tposition[I][p1->level]];
 }///多极矩贡献删掉

int n_t=-1;tnode*p2;
if(temp==-1)return;
for(int i=0;i<=7;i++)
 {
  p2=p_temp[temp]->child[i];
  for(int j=0;j<=p2->numpar-1;j++)
    {
      n_t++;
      p_temp[temp]->particle[n_t]=p2->particle[j];
      for(int k=temp;k>=0;k--)index_p[p_temp[k]->particle[n_t]][p_temp[k]->level]=n_t;
    }
 }
for(int i=temp-1;i>=0;i--)for(int j=0;j<=N0-1;j++)
 {
     p_temp[i]->particle[j]=p_temp[temp]->particle[j];
 }
}
//!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void init_tree()///external
{
for(int i=0;i<=9999;i++)for(int j=0;j<=9;j++)tposition[i][j]=-1;
root.x_max=tMax_x;root.y_max=tMax_y;root.z_max=tMax_z;
root.x_min=tMin_x;root.y_min=tMin_y;root.z_min=tMin_z;
root.x_mid=(tMax_x+tMin_x)/2;root.y_mid=(tMax_y+tMin_y)/2;root.z_mid=(tMax_z+tMin_z)/2;
root.sqradius=((tMax_x-tMin_x)*(tMax_x-tMin_x)+(tMax_y-tMin_y)*(tMax_y-tMin_y)+(tMax_z-tMin_z)*(tMax_z-tMin_z))/4;
root.radius=sqrt(root.sqradius);
root.level=0;
expand(&root);
for(int i=0;i<=tNpar;i++)insert_particle(i,&root);
}
//!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//extern double T_x,T_y,Max_x,Min_x,Max_y,Min_y;
void move_tree(int I,double dx,double dy,double dz)///移动离子并更新树
{
double x_temp,y_temp,z_temp;
x_temp=x[I]+dx;y_temp=y[I]+dy;z_temp=z[I]+dz;
//if((abs(x_temp-Max_x)<1e-8)||(x_temp>Max_x))x_temp-=T_x;if(x_temp<Min_x)x_temp+=T_x;
//if((abs(y_temp-Max_y)<1e-8)||(y_temp>Max_y))y_temp-=T_y;if(y_temp<Min_y)y_temp+=T_y;
//if(abs(x_temp-Max_x)<1e-8)x_temp=Min_x;if(x_temp>Max_x)x_temp-=T_x;if(x_temp<Min_x)x_temp+=T_x;
//if(abs(y_temp-Max_y)<1e-8)y_temp=Min_y;if(y_temp>Max_y)y_temp-=T_y;if(y_temp<Min_y)y_temp+=T_y;
tnode*p=&root;tnode*oldp=p;
while((p->x_min<=x_temp)&&(x_temp<p->x_max)&&(p->y_min<=y_temp)&&(y_temp<p->y_max)&&(p->z_min<z_temp)&&(z_temp<p->z_max))
{
    for(int i=0;i<=MAXP;i++)
     for(int j=0;j<=MAXP-i;j++)
      for(int k=0;k<=MAXP-i-j;k++)
       p->moment[i][j][k]+=-q[I]*pow(x[I]-p->x_mid,i)*pow(y[I]-p->y_mid,j)*pow(z[I]-p->z_mid,k)+q[I]*pow(x_temp-p->x_mid,i)*pow(y_temp-p->y_mid,j)*pow(z_temp-p->z_mid,k);

    if(p->level==MAXLEVEL){x[I]=x_temp;y[I]=y_temp;z[I]=z_temp;return;}
    oldp=p;
    p=p->child[tposition[I][p->level]];
}
p=oldp;///这样p就指向了修改前的正确的最深节点

if(oldp->level<MAXLEVEL)
{
    oldp=oldp->child[tposition[I][p->level]];
    delete_particle(I,oldp);
}
 ///将错误的结点进行修改
 x[I]=x_temp;y[I]=y_temp;z[I]=z_temp;
 for(int i=0;i<=7;i++)if(insert_particle(I,p->child[i]))tposition[I][p->level]=i;///将离子插入正确的子树

}
//!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Energy_tree(int I,tnode*p,double Dx,double Dy,double Dz)//计算以p->的内容为根的树中离子在第I个离子处的能量
{
if(p->numpar==0)return 0;
double ans=0;
double R,Rx,Ry,Rz,r;
Rx=p->x_mid-x[I]+Dx;
Ry=p->y_mid-y[I]+Dy;
Rz=p->z_mid-z[I]+Dz;
R=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
r=p->radius/R;
double T[10][10][10];
double t11,t12,t13,t21,t22,t23;
double dx,dy,dz,q_temp;
///if((pow(r,MAXP+1)/(1-r)<=err*R)&&(r<1))//满足级数近似条件
double inv_R2=1.0/(R*R),inv_R3=1.0/(R*R*R),inv_R5=1.0/(R*R*R*R*R);
if((r<=BAC)&&(p->numpar>N0))
 {
    T[0][0][0]=1.0/R;
    ans+=T[0][0][0]*(p->moment[0][0][0]);
    if(MAXP>=1)
    {
    T[1][0][0]=-Rx*inv_R3;
    T[0][1][0]=-Ry*inv_R3;
    T[0][0][1]=-Rz*inv_R3;
    ans+=(T[1][0][0]*(p->moment[1][0][0])+T[0][1][0]*(p->moment[0][1][0])+T[0][0][1]*(p->moment[0][0][1]));
    }
    if(MAXP>=2)
    {
        T[2][0][0]=(3*Rx*Rx*inv_R2-1)*0.5*inv_R3;
        T[0][2][0]=(3*Ry*Ry*inv_R2-1)*0.5*inv_R3;
        T[0][0][2]=(3*Rz*Rz*inv_R2-1)*0.5*inv_R3;
        T[1][1][0]=(3*Rx*Ry)*inv_R5;
        T[1][0][1]=(3*Rx*Rz)*inv_R5;
        T[0][1][1]=(3*Ry*Rz)*inv_R5;
        ans+=(T[2][0][0]*(p->moment[2][0][0])+T[0][2][0]*(p->moment[0][2][0])+T[0][0][2]*(p->moment[0][0][2])+T[1][1][0]*(p->moment[1][1][0])+T[1][0][1]*(p->moment[1][0][1])+T[0][1][1]*(p->moment[0][1][1]));
    }
    for(int P=3;P<=MAXP;P++)
     for(int i=0;i<=P;i++)
      for(int j=0;j<=P-i;j++)
      {
          int k=P-i-j;
          t11=(i-1<0)?0:T[i-1][j][k];
          t12=(j-1<0)?0:T[i][j-1][k];
          t13=(k-1<0)?0:T[i][j][k-1];
          t21=(i-2<0)?0:T[i-2][j][k];
          t22=(j-2<0)?0:T[i][j-2][k];
          t23=(k-2<0)?0:T[i][j][k-2];
          T[i][j][k]=-1.0/(P*R*R)*((2*P-1)*(Rx*t11+Ry*t12+Rz*t13)+(P-1)*(t21+t22+t23));
          ans+=T[i][j][k]*(p->moment[i][j][k]);
      }
 }
else if((p->level==MAXLEVEL)||(p->numpar<=N0))//不满足级数近似条件但已经达到叶子,或者该节点所含粒子数小于等于N0
 {
     for(int i=0;i<=p->numpar-1;i++)
      if((I!=p->particle[i])||(Dx!=0)||(Dy!=0)||(Dz!=0))
       {
        dx=x[p->particle[i]]-x[I]+Dx;
        dy=y[p->particle[i]]-y[I]+Dy;
        dz=z[p->particle[i]]-z[I]+Dz;
        q_temp=q[p->particle[i]];
        ans+=q_temp/sqrt(dx*dx+dy*dy+dz*dz);
       }
 }
else//不满足近似条件，继续往下分
 {
     for(int i=0;i<=7;i++)ans+=Energy_tree(I,p->child[i],Dx,Dy,Dz);
 }
return ans;
}
