#define MAX_Npar_incell 80
struct tnode
{
	int numpar;///numpar�����ӣ���ôparticle�±��Ǵ�0~numpar-1
	double x_min,y_min,z_min;
	double x_max,y_max,z_max;
	double x_mid,y_mid,z_mid;
	double radius,sqradius;
	int level;
	struct tnode *child[8];
    //struct tnode *father;
    double moment[3][3][3];
    int *particle;//[MAX_Npar_incell];///�������������������㹻������particle[i]��ʾ�ڼ������
};

