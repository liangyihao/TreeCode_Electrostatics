#define MAX_Npar_incell 80
struct tnode
{
	int numpar;///numpar个离子，那么particle下标是从0~numpar-1
	double x_min,y_min,z_min;
	double x_max,y_max,z_max;
	double x_mid,y_mid,z_mid;
	double radius,sqradius;
	int level;
	struct tnode *child[8];
    //struct tnode *father;
    double moment[3][3][3];
    int *particle;//[MAX_Npar_incell];///这个变量仅在树的深度足够后开启，particle[i]表示第几个电荷
};

